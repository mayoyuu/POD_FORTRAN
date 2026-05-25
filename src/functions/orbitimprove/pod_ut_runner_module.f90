!> @file pod_ut_runner_module.f90
!> @brief 无迹变换（UT）轨道确定集成封装层
module pod_ut_runner_module
    use pod_global, only: DP, MAX_STRING_LEN
    use pod_filter_ut_module, only: ut_filter
    use pod_obs_io_module, only: obs_record, preload_observations, &
                                  station_record, preload_stations, find_station_by_id, &
                                  ref_orbit_record, preload_reference_orbits
    use pod_measurement_base_module, only: observation_station
    use pod_basicmath_module, only: PI
    use pod_data_format_module, only: load_initial_opm, write_json_opm, &
                                       write_residual_line, write_error_line
    use pod_error_analysis_module, only: compute_orbit_error

    implicit none
    private

    public :: run_ut_orbit_determination

contains

    !> ======================================================================
    !> 执行基于无迹卡尔曼滤波（UKF）的轨道确定流程
    !> ======================================================================
    subroutine run_ut_orbit_determination(obs_file, site_json_file, ref_orbit_file, &
                                          initial_json_file, output_opm_file, &
                                          output_residual_file, output_error_file)
        character(len=*), intent(in) :: obs_file
        character(len=*), intent(in) :: site_json_file
        character(len=*), intent(in) :: initial_json_file
        character(len=*), intent(in) :: output_opm_file
        character(len=*), intent(in) :: output_residual_file
        character(len=*), intent(in), optional :: ref_orbit_file
        character(len=*), intent(in), optional :: output_error_file

        ! 核心对象
        type(ut_filter)       :: my_filter
        type(observation_station) :: current_station
        type(obs_record), allocatable :: obs_list(:)
        type(station_record), allocatable :: station_list(:)
        type(ref_orbit_record), allocatable :: ref_list(:)

        ! 状态与时间
        real(DP) :: initial_mean(6), final_mean(6)
        real(DP) :: initial_cov(6,6), final_cov(6,6), noise_Q(6,6)
        real(DP) :: y_meas(2), noise_R(2,2)
        real(DP) :: et_current, et_obs, dt
        integer  :: obs_count, i
        real(DP), parameter :: sigma_a = 1.0e-11_DP   ! km/s²

        real(DP) :: step_res(6)    ! 存储最近一次测量更新的 6 列残差
        real(DP) :: step_comp(2)   ! 存储最近一次计算出的预测观测值 (Lon, Lat)
        real(DP) :: pos_err(3), vel_err(3)
        real(DP) :: pos_rms, vel_rms, mahalanobis_d
        

        ! ---- 1. 测量噪声协方差（光学赤经赤纬，0.1角秒精度）----
        noise_R = 0.0_DP
        noise_R(1,1) = (0.1_DP * PI / 180.0_DP / 3600.0_DP)**2
        noise_R(2,2) = noise_R(1,1)

        ! ---- 2. 从 JSON 加载初始状态（仅均值和协方差）----
        call load_initial_opm(initial_json_file, et_current, initial_mean, initial_cov)

        ! ---- 3. 初始化 UT 滤波器 ----
        call my_filter%filter_init(et_current, initial_mean, initial_cov)

        write(*,*) '  [UT Runner] 滤波器初始化完成'
        write(*,*) '    初始历元：', et_current
        write(*,*) '    初始均值：', initial_mean(1:3)

        ! ---- 4. 预加载全部观测、测站与参考轨道 ----
        call preload_observations(obs_file, obs_list)
        call preload_stations(site_json_file, station_list)
        if (present(ref_orbit_file)) then
            call preload_reference_orbits(ref_orbit_file, ref_list)
        end if

        if (present(ref_orbit_file)) then
            write(*,*) '  [UT Runner] 预加载完成：观测 ', size(obs_list), &
                       ' 条，测站 ', size(station_list), ' 个，参考轨道 ', size(ref_list), ' 条'
        else
            write(*,*) '  [UT Runner] 预加载完成：观测 ', size(obs_list), &
                       ' 条，测站 ', size(station_list), ' 个'
        end if

        ! ---- 5. 主循环：逐观测进行时间更新 + 测量更新 ----
        do obs_count = 1, size(obs_list)
            et_obs = obs_list(obs_count)%et
            y_meas(1) = obs_list(obs_count)%ra
            y_meas(2) = obs_list(obs_count)%dec
            current_station = find_station_by_id(obs_list(obs_count)%station_id, station_list)

            ! 当前滤波器时刻
            call my_filter%get_current_epoch(et_current)
            dt = et_obs - et_current

            ! 构造与时间步长相关的过程噪声 Q
            noise_Q = 0.0_DP
            do i = 1, 3
                noise_Q(i,i)       = (dt**4 / 4.0_DP) * sigma_a**2
                noise_Q(i+3,i+3)   = dt**2 * sigma_a**2
            end do

            write(*,*) '  [UT Runner] 处理观测 #', obs_count, &
                    '  观测时刻 et_obs = ', et_obs, ' 秒', &
                    '  时间步长 dt =', dt, ' 秒'

            ! call my_filter%time_update(et_obs, noise_Q)
            call my_filter%time_update(et_obs)
            call my_filter%get_current_epoch(et_current)
            call my_filter%get_current_state(final_mean)
            call my_filter%get_current_cov(final_cov)

            if (present(ref_orbit_file)) then
                call compute_orbit_error(final_mean, final_cov, ref_list(obs_count)%state, &
                                          pos_err, vel_err, pos_rms, vel_rms, mahalanobis_d)

                if (present(output_error_file)) then
                    call write_error_line(output_error_file, et_obs, pos_err, vel_err, &
                                           pos_rms, vel_rms, mahalanobis_d, (obs_count == 1))
                end if
            end if

            call my_filter%measurement_update(y_meas, noise_R, et_obs, current_station)
            call my_filter%get_current_epoch(et_current)
            call my_filter%get_current_state(final_mean)

            call my_filter%get_last_residual(step_res, step_comp)

            call write_residual_line(output_residual_file, et_obs, y_meas, step_comp, &
                                    step_res, trim(current_station%name), (obs_count == 1))

        end do

        write(*,*) '  [UT Runner] 滤波结束，有效观测数：', obs_count - 1

        ! ---- 5. 提取最终结果并写入 JSON 文件 ----
        call my_filter%get_current_epoch(et_current)
        call my_filter%get_current_state(final_mean)
        call my_filter%get_current_cov(final_cov)

        ! 使用不带 GMM 状态参数的 write_json_opm 重载版本（若库中提供）
        ! 若库中仅有带 gmm_state 的版本，可改为：
        !   call write_json_opm(output_file_name, final_mean, final_cov, et_current)
        call write_json_opm(output_opm_file, final_mean, final_cov, &
                    rms = 0.0_DP, obj_id = "DRO", et_last = et_current)

        write(*,*) '  [UT Runner] 结果已写入：', trim(output_opm_file)

    end subroutine run_ut_orbit_determination

end module pod_ut_runner_module