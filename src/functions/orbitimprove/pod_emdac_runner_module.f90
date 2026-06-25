!> @file pod_emdac_runner_module.f90
!> @brief EMDAC-N 轨道改进集成封装层
module pod_emdac_runner_module
    use pod_global, only: DP, MAX_STRING_LEN
    use pod_config, only: config
    use pod_uq_gmm_state_module, only: uq_gmm_state_type
    use pod_filter_emdac_module, only: emdac_filter
    use pod_obs_io_module, only: obs_record, preload_observations, &
                                  station_record, preload_stations, find_station_by_id, &
                                  ref_orbit_record, preload_reference_orbits, find_reference_by_et
    use pod_measurement_base_module, only: observation_station
    use pod_basicmath_module, only: PI
    use pod_data_format_module, only: load_initial_opm, write_json_opm, &
                                       write_residual_line, write_error_line, &
                                       write_gmm_diag_lines, write_gmm_snapshots
    use pod_error_analysis_module, only: compute_orbit_error, &
                                       compute_covariance_condition_number

    implicit none
    private
    
    ! 仅对外暴露这个唯一的极简接口
    public :: run_emdac_orbit_determination

contains

    !> ======================================================================
    !> 核心集成接口：执行完整的 EMDAC 轨道定轨流程
    !> ======================================================================
    subroutine run_emdac_orbit_determination(obs_file, site_json_file, ref_orbit_file, &
                                             initial_json_file, output_opm_file, &
                                             output_residual_file, output_error_file, &
                                             gmm_in_switch, n_components, max_da_order, &
                                             opt_particles, opt_em_max_iter, opt_em_tol, &
                                             output_gmm_file)

        character(len=*), intent(in) :: obs_file             ! 观测数据文件 (.obs)
        character(len=*), intent(in) :: site_json_file       ! 测站配置文件 (.json)
        character(len=*), intent(in) :: initial_json_file    ! 初始先验状态文件 (.json)
        character(len=*), intent(in) :: output_opm_file      ! 输出定轨结果文件
        character(len=*), intent(in) :: output_residual_file ! 输出残差文件
        character(len=*), intent(in), optional :: ref_orbit_file    ! 参考轨道文件 (.ref)
        character(len=*), intent(in), optional :: output_error_file ! 输出误差文件 (.err)
        character(len=*), intent(in), optional :: output_gmm_file   ! GMM 快照输出文件 (.gmms.json)
        logical, intent(in) :: gmm_in_switch
        integer,  intent(in) :: n_components
        integer,  intent(in) :: max_da_order

        logical :: has_gmm_loaded
        integer,  intent(in), optional :: opt_particles      ! 粒子总数
        integer,  intent(in), optional :: opt_em_max_iter    ! EM 算法最大迭代次数
        real(DP), intent(in), optional :: opt_em_tol         ! EM 算法收敛容差
    
        ! 局部变量：核心对象
        type(emdac_filter) :: my_filter
        type(observation_station) :: current_station
        type(obs_record), allocatable :: obs_list(:)
        type(station_record), allocatable :: station_list(:)
        type(ref_orbit_record), allocatable :: ref_list(:)
        
        ! 状态与时间变量
        real(DP) :: initial_mean(6), final_mean(6)
        real(DP) :: initial_cov(6,6), final_cov(6,6), noise_Q(6,6)
        type(uq_gmm_state_type) :: initial_gmm

        real(DP) :: y_meas(2), noise_R(2,2)
        real(DP) :: et_current, et_obs, dt
        integer :: obs_count, i
        
        ! 【新增】用于自适应阶数的内部变量，绝不污染外层接口
        integer :: current_order
        logical :: is_first_step

        real(DP) :: step_res(6)    ! 存储最近一次测量更新的 6 列残差
        real(DP) :: step_comp(2)   ! 存储最近一次计算出的预测观测值 (Lon, Lat)
        real(DP) :: prior_pos_err(3), prior_vel_err(3)
        real(DP) :: posterior_pos_err(3), posterior_vel_err(3)
        real(DP) :: prior_pos_rms, prior_vel_rms, prior_mahalanobis_d
        real(DP) :: posterior_pos_rms, posterior_vel_rms, posterior_mahalanobis_d
        real(DP) :: prior_cond_p, prior_lambda_min_p, prior_lambda_max_p
        real(DP) :: posterior_cond_p, posterior_lambda_min_p, posterior_lambda_max_p
        real(DP) :: ref_state_at_obs(6)
        integer :: prior_cov_info, posterior_cov_info
        logical :: ref_found

        ! GMM 快照收集
        type(uq_gmm_state_type), allocatable :: gmm_snapshots(:)
        real(DP), allocatable :: gmm_epochs(:)
        logical :: collect_gmms
        type(uq_gmm_state_type) :: step_gmm
        logical :: is_first_gmm_diag
        
        ! 1. 测量噪声协方差设置 (例如光学赤经赤纬，0.1角秒精度)
        noise_R = 0.0_DP
        noise_R(1,1) = (config%measurement_noise_ra_arcsec * PI / 180.0_DP / 3600.0_DP)**2
        noise_R(2,2) = (config%measurement_noise_dec_arcsec * PI / 180.0_DP / 3600.0_DP)**2

        ! 2. 初始状态加载 (从文件读取先验值)
        if (gmm_in_switch) then
            call load_initial_opm(initial_json_file, et_current, initial_mean, initial_cov, initial_gmm, has_gmm_loaded)
            call my_filter%init(et_current, initial_mean,initial_cov, n_components, max_da_order,initial_gmm = initial_gmm)
        else
            ! 从 JSON 文件加载单一高斯先验状态
            call load_initial_opm(initial_json_file, et_current, initial_mean, initial_cov)
            call my_filter%init(et_current, initial_mean,initial_cov, n_components, max_da_order)
        end if

        write(*,*) '  [Runner] 滤波器初始化完成'
        write(*,*) '    初始均值: ', initial_mean
        write(*,*) '    初始历元: ', et_current
        
        ! 3. 滤波器装配与初始化
        if (present(opt_particles)) then
            call my_filter%set_n_particles(opt_particles)
        end if
        if (present(opt_em_tol)) then
            call my_filter%set_em_parameters(tol = opt_em_tol)
        end if
        if (present(opt_em_max_iter)) then
            call my_filter%set_em_parameters(max_iter = opt_em_max_iter)
        end if
        
        ! 初始化阶数控制逻辑标志
        is_first_step = .true.
        is_first_gmm_diag = .true.
        
        ! 4. 预加载全部观测、测站与参考轨道
        call preload_observations(obs_file, obs_list)
        call preload_stations(site_json_file, station_list)
        if (present(ref_orbit_file)) then
            call preload_reference_orbits(ref_orbit_file, ref_list)
        end if

        if (present(ref_orbit_file)) then
            write(*,*) '  [Runner] 预加载完成：观测 ', size(obs_list), &
                       ' 条，测站 ', size(station_list), ' 个，参考轨道 ', size(ref_list), ' 条'
        else
            write(*,*) '  [Runner] 预加载完成：观测 ', size(obs_list), &
                       ' 条，测站 ', size(station_list), ' 个'
        end if

        ! GMM 快照收集初始化
        collect_gmms = present(output_gmm_file)
        if (collect_gmms) then
            allocate(gmm_snapshots(size(obs_list)))
            allocate(gmm_epochs(size(obs_list)))
            write(*,*) '  [Runner] GMM 快照输出启用，将记录 ', size(obs_list), ' 步'
        end if

        ! 5. 核心数据同化流 (Time Update + Measurement Update)
        do obs_count = 1, size(obs_list)
            et_obs = obs_list(obs_count)%et
            y_meas(1) = obs_list(obs_count)%ra
            y_meas(2) = obs_list(obs_count)%dec
            current_station = find_station_by_id(obs_list(obs_count)%station_id, station_list)

            write(*,*) 'Station ECEF:', current_station%ecef_position
            write(*,*) 'Station geodetic:', current_station%latitude, current_station%longitude, current_station%altitude
            
            call my_filter%get_current_epoch(et_current)
            write(*,*) '  传播前时间为: ', et_current
            call my_filter%get_current_state(final_mean)
            write(*,*) '  传播前位置为: ',  final_mean

            ! 计算积分步长
            dt = et_obs - et_current

            ! 给出noise_Q
            noise_Q = 0.0_DP
            if (config%use_process_noise) then
                do i = 1, 3
                    noise_Q(i,i)     = (dt**4 / 4.0_DP) * config%process_noise_sigma_a**2
                    noise_Q(i+3,i+3) = dt**2 * config%process_noise_sigma_a**2
                end do
            end if
            ! ==========================================================
            ! 智能 DA 阶数调整逻辑 (完全基于步长时间判定)
            ! ==========================================================
           if (is_first_step) then
                ! 1. 如果是第一步，且用户指定了 opt_da_order，无条件遵从用户输入
                current_order = max_da_order
            else
                ! 2. 后续步骤 (或用户没指定的第一步)，走基于步长 dt 的智能选择逻辑
                if (abs(dt) > 3.0_DP * 86400.0_DP) then      ! 步长大于3天 -> 最大阶数
                    current_order = max_da_order
                else if (abs(dt) > 86400.0_DP) then          ! 步长在1天到3天之间 -> 3阶
                    current_order = 3
                else if (abs(dt) < 0.5*3600.0_DP) then           ! 步长小于0.5小时 -> 1阶
                    current_order = 1
                else                                         ! 步长在1小时到1天之间 -> 2阶
                    current_order = 2
                end if
            end if
            ! 应用最新计算出的阶数
            call my_filter%set_da_order(current_order)
            is_first_step = .false. ! 第一步已走完，切断强制覆盖机制
            ! ==========================================================
            
            write(*,'(A,I0,A,F10.2,A,I1)') '  [Runner] 处理观测 #', obs_count, &
                  ' dt:', dt, 's, DA阶数:', current_order
            
            if (config%use_process_noise) then
                call my_filter%time_update(et_obs, noise_Q)
            else
                call my_filter%time_update(et_obs)
            end if
            call my_filter%get_current_epoch(et_current)
            write(*,*) '  传播后时间为: ', et_current
            call my_filter%get_current_state(final_mean)
            write(*,*) '  时间更新后结果为: ',  final_mean
            call my_filter%get_current_cov(final_cov)

            if (present(ref_orbit_file)) then
                call find_reference_by_et(et_obs, ref_list, ref_state_at_obs, ref_found, tolerance=1.0_DP)
                if (ref_found) then
                    call compute_orbit_error(final_mean, final_cov, ref_state_at_obs, &
                                              prior_pos_err, prior_vel_err, &
                                              prior_pos_rms, prior_vel_rms, &
                                              prior_mahalanobis_d)
                    prior_cond_p = -1.0_DP
                    prior_lambda_min_p = -1.0_DP
                    prior_lambda_max_p = -1.0_DP
                    call compute_covariance_condition_number(final_cov, prior_cond_p, prior_cov_info, &
                                                             lambda_min=prior_lambda_min_p, &
                                                             lambda_max=prior_lambda_max_p)
                    if (present(output_error_file)) then
                        call my_filter%get_current_gmm(step_gmm)
                        call write_gmm_diag_lines(output_error_file, et_obs, "PRIOR", &
                                                  step_gmm, ref_state_at_obs, is_first_gmm_diag)
                        is_first_gmm_diag = .false.
                    end if
                else
                    write(*,*) '[WARN] No reference orbit found within 1 second of observation ET: ', et_obs
                end if
            end if

            call my_filter%measurement_update(y_meas, noise_R, et_obs, current_station)
            call my_filter%get_current_epoch(et_current)
            write(*,*) '  测量更新后时间为: ', et_current
            call my_filter%get_current_state(final_mean)
            write(*,*) '  测量更新后结果为: ',  final_mean

            call my_filter%get_current_cov(final_cov)

            if (present(ref_orbit_file)) then
                if (ref_found) then
                    call compute_orbit_error(final_mean, final_cov, ref_state_at_obs, &
                                              posterior_pos_err, posterior_vel_err, &
                                              posterior_pos_rms, posterior_vel_rms, &
                                              posterior_mahalanobis_d)
                    posterior_cond_p = -1.0_DP
                    posterior_lambda_min_p = -1.0_DP
                    posterior_lambda_max_p = -1.0_DP
                    call compute_covariance_condition_number(final_cov, posterior_cond_p, &
                                                             posterior_cov_info, &
                                                             lambda_min=posterior_lambda_min_p, &
                                                             lambda_max=posterior_lambda_max_p)

                    if (present(output_error_file)) then
                        call write_error_line(output_error_file, et_obs, prior_pos_err, prior_vel_err, &
                                               prior_pos_rms, prior_vel_rms, prior_mahalanobis_d, &
                                               (obs_count == 1), &
                                               posterior_mahalanobis_d=posterior_mahalanobis_d, &
                                               prior_cond_p=prior_cond_p, &
                                               prior_lambda_min_p=prior_lambda_min_p, &
                                               prior_lambda_max_p=prior_lambda_max_p, &
                                               posterior_cond_p=posterior_cond_p, &
                                               posterior_lambda_min_p=posterior_lambda_min_p, &
                                               posterior_lambda_max_p=posterior_lambda_max_p)
                    end if

                    if (present(output_error_file)) then
                        call my_filter%get_current_gmm(step_gmm)
                        call write_gmm_diag_lines(output_error_file, et_obs, "POSTERIOR", &
                                                  step_gmm, ref_state_at_obs, is_first_gmm_diag)
                        is_first_gmm_diag = .false.
                    end if
                end if
            end if

            call my_filter%get_last_residual(step_res, step_comp)

            ! 调用写入函数
            ! obs_count == 1 时为 true，创建新文件；之后为 false，追加写入
            call write_residual_line(output_residual_file, et_obs, y_meas, step_comp, &
                                    step_res, trim(current_station%name), (obs_count == 1))

            ! 收集测量更新后的 GMM 快照
            if (collect_gmms) then
                call my_filter%get_current_gmm(gmm_snapshots(obs_count))
                gmm_epochs(obs_count) = et_obs
            end if

        end do
        
        write(*,*) '  [Runner] 滤波结束，有效观测数: ', obs_count - 1

        ! 写入 GMM 快照文件
        if (collect_gmms) then
            call write_gmm_snapshots(output_gmm_file, gmm_epochs, gmm_snapshots, size(obs_list))
            deallocate(gmm_snapshots, gmm_epochs)
        end if

        ! 5. 结果提取与落盘
        call my_filter%get_current_epoch(et_current)
        call my_filter%get_current_state(final_mean)
        call my_filter%get_current_cov(final_cov)
        
        call write_json_opm(output_opm_file, final_mean, final_cov, my_filter%gmm_state, 0.0_DP, "DRO", et_current)
        
    end subroutine run_emdac_orbit_determination

end module pod_emdac_runner_module