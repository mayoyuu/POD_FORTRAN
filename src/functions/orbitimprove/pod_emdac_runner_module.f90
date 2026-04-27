!> @file pod_emdac_runner_module.f90
!> @brief EMDAC-N 轨道改进集成封装层
module pod_emdac_runner_module
    use pod_global, only: DP, MAX_STRING_LEN
    use pod_uq_gmm_state_module, only: uq_gmm_state_type
    use pod_filter_emdac_module, only: emdac_filter
    use pod_obs_io_module, only: load_single_observation
    use pod_measurement_base_module, only: observation_station
    use pod_basicmath_module, only: PI
    ! 假设引入了 JSON 读写模块
    use pod_data_format_module, only: load_initial_opm, write_json_opm

    implicit none
    private
    
    ! 仅对外暴露这个唯一的极简接口
    public :: run_emdac_orbit_determination

contains

    !> ======================================================================
    !> 核心集成接口：执行完整的 EMDAC 轨道定轨流程
    !> ======================================================================
    subroutine run_emdac_orbit_determination(obs_file, site_json_file, gmm_in_switch, &
                                             initial_json_file, output_json_file, n_components,&
                                             max_da_order,opt_particles, &
                                             opt_em_max_iter, opt_em_tol)
        
        character(len=*), intent(in) :: obs_file           ! 观测文件路径 (.obs)
        character(len=*), intent(in) :: site_json_file     ! 测站配置文件路径 (.json)
        character(len=*), intent(in) :: initial_json_file  ! 初始先验状态文件路径 (.opm/.json)
        character(len=*), intent(in) :: output_json_file   ! 输出定轨结果文件路径 (.opm/.json)
        logical, intent(in) :: gmm_in_switch               ! GMM 初始化开关
        integer,  intent(in) :: n_components       ! GMM 分量数量
        integer,  intent(in) :: max_da_order       ! DA 阶数
    
        logical :: has_gmm_loaded
        integer,  intent(in), optional :: opt_particles      ! 粒子总数
        integer,  intent(in), optional :: opt_em_max_iter    ! EM 算法最大迭代次数
        real(DP), intent(in), optional :: opt_em_tol         ! EM 算法收敛容差
    
        ! 局部变量：核心对象
        type(emdac_filter) :: my_filter
        type(observation_station) :: current_station
        
        ! 状态与时间变量
        real(DP) :: initial_mean(6), final_mean(6)
        real(DP) :: initial_cov(6,6), final_cov(6,6)
        type(uq_gmm_state_type) :: initial_gmm

        real(DP) :: y_meas(2), noise_R(2,2)
        real(DP) :: et_current, et_obs, dt
        logical :: is_eof
        integer :: obs_count, i
        
        ! 【新增】用于自适应阶数的内部变量，绝不污染外层接口
        integer :: current_order
        logical :: is_first_step
        
        ! 1. 测量噪声协方差设置 (例如光学赤经赤纬，0.1角秒精度)
        noise_R = 0.0_DP
        noise_R(1,1) = (0.1_DP * PI / 180.0_DP / 3600.0_DP)**2 
        noise_R(2,2) = noise_R(1,1)

        ! 2. 初始状态加载 (从文件读取先验值)
        if (gmm_in_switch) then
            call load_initial_opm(initial_json_file, et_current, initial_mean, initial_cov, initial_gmm, has_gmm_loaded)
            call my_filter%init(et_current, initial_mean,initial_cov, n_components, max_da_order,initial_gmm = initial_gmm)
        else
            ! 从 JSON 文件加载单一高斯先验状态
            call load_initial_opm(initial_json_file, et_current, initial_mean, initial_cov)
            call my_filter%init(et_current, initial_mean,initial_cov, n_components, max_da_order)
        end if
        
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
        
        ! 4. 核心数据同化流 (Time Update + Measurement Update)
        obs_count = 1
        do
            call load_single_observation(obs_file, site_json_file, obs_count, &
                                         et_obs, y_meas(1), y_meas(2), current_station, is_eof)
            if (is_eof) exit
            
            ! 计算积分步长
            dt = et_obs - et_current
            
            ! ==========================================================
            ! 智能 DA 阶数调整逻辑 (完全基于步长时间判定)
            ! ==========================================================
            if (is_first_step) then
                ! 1. 如果是第一步，且用户指定了 opt_da_order，无条件遵从用户输入
                current_order = max_da_order
            else
                ! 2. 后续步骤 (或用户没指定的第一步)，走基于步长 dt 的智能选择逻辑
                if (abs(dt) > 86400.0_DP) then
                    current_order = 3    ! 步长超过1天 -> 3阶
                else if (abs(dt) < 3600.0_DP) then
                    current_order = 1    ! 步长小于1小时 -> 1阶
                else
                    current_order = 2    ! 步长在1小时到1天之间 -> 2阶
                end if
            end if
            
            ! 应用最新计算出的阶数
            call my_filter%set_da_order(current_order)
            is_first_step = .false. ! 第一步已走完，切断强制覆盖机制
            ! ==========================================================
            
            write(*,'(A,I0,A,F10.2,A,I1)') '  [Runner] 处理观测 #', obs_count, &
                  ' dt:', dt, 's, DA阶数:', current_order
            
            call my_filter%time_update(et_obs)
            call my_filter%measurement_update(y_meas, noise_R, et_obs, current_station)

            et_current = et_obs
            obs_count = obs_count + 1
        end do
        
        write(*,*) '  [Runner] 滤波结束，有效观测数: ', obs_count - 1
        
        ! 5. 结果提取与落盘
        call my_filter%get_current_epoch(et_current)
        call my_filter%get_current_state(final_mean)
        call my_filter%get_current_cov(final_cov)
        
        call write_json_opm(output_json_file, final_mean, final_cov, my_filter%gmm_state, 0.0_DP, "CAT_TARGET", et_current)
        
    end subroutine run_emdac_orbit_determination

end module pod_emdac_runner_module