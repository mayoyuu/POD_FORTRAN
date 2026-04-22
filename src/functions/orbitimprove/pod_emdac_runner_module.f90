!> @file pod_emdac_runner_module.f90
!> @brief EMDAC-N 轨道改进集成封装层
module pod_emdac_runner_module
    use pod_global, only: DP, MAX_STRING_LEN
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
    subroutine run_emdac_orbit_determination(obs_file, site_json_file, &
                                             initial_json_file, output_json_file, &
                                             opt_particles, opt_da_order, &
                                             opt_em_max_iter, opt_em_tol, n_components)
        
        character(len=*), intent(in) :: obs_file           ! 观测文件路径 (.obs)
        character(len=*), intent(in) :: site_json_file     ! 测站配置文件路径 (.json)
        character(len=*), intent(in) :: initial_json_file  ! 初始先验状态文件路径 (.opm/.json)
        character(len=*), intent(in) :: output_json_file   ! 输出定轨结果文件路径 (.opm/.json)
        
        integer,  intent(in) :: opt_particles      ! 粒子总数
        integer,  intent(in) :: opt_da_order       ! DA 阶数
        integer,  intent(in) :: opt_em_max_iter    ! EM 算法最大迭代次数
        real(DP), intent(in) :: opt_em_tol         ! EM 算法收敛容差
        integer,  intent(in) :: n_components       ! GMM 分量数量
        
        ! 局部变量：核心对象
        type(emdac_filter) :: my_filter
        type(observation_station) :: current_station
        
        ! 状态与时间变量
        real(DP) :: initial_mean(6), final_mean(6)
        real(DP) :: initial_cov(6,6), final_cov(6,6)
        real(DP) :: y_meas(2), noise_R(2,2)
        real(DP) :: et_current, et_obs, dt
        logical :: is_eof
        integer :: obs_count, i
        
        ! 1. 测量噪声协方差设置 (例如光学赤经赤纬，0.1角秒精度)
        noise_R = 0.0_DP
        noise_R(1,1) = (0.1_DP * PI / 180.0_DP / 3600.0_DP)**2 
        noise_R(2,2) = noise_R(1,1)

        ! 2. 初始状态加载 (从文件读取先验值)
        ! call load_initial_opm(initial_json_file, et_current, initial_mean, initial_cov)
        et_current = 0.0_DP 
        initial_mean = 0.0_DP
        initial_cov = 0.0_DP; do i=1,6; initial_cov(i,i) = 1.0e-6_DP; end do
        
        ! 3. 滤波器装配与初始化
        my_filter%n_particles = opt_particles
        my_filter%em_tol = opt_em_tol
        my_filter%em_max_iter = opt_em_max_iter
        
        call my_filter%set_da_order(opt_da_order)
        call my_filter%init(initial_mean, initial_cov, n_components, opt_particles)
        
        ! 4. 核心数据同化流 (Time Update + Measurement Update)
        obs_count = 1
        do
            call load_single_observation(obs_file, site_json_file, obs_count, &
                                         et_obs, y_meas(1), y_meas(2), current_station, is_eof)
            if (is_eof) exit
            
            write(*,'(A,I0,A,F14.3)') '  [Runner] 处理观测 #', obs_count, ' 历元: ', et_obs
            
            dt = et_obs - et_current
            if (abs(dt) > 1.0e-6_DP) then
                call my_filter%time_update(0.0_DP, dt)
            end if
            
            call my_filter%measurement_update(y_meas, noise_R, et_obs, current_station)
            
            et_current = et_obs
            obs_count = obs_count + 1
        end do
        
        write(*,*) '  [Runner] 滤波结束，有效观测数: ', obs_count - 1
        
        ! 5. 结果提取与落盘
        final_mean = my_filter%state_mean
        ! call my_filter%gmm_state%compute_global_covariance(final_cov) 
        
        call write_json_opm(output_json_file, final_mean, final_cov, 0.0_DP, "CAT_TARGET", et_current)
        
    end subroutine run_emdac_orbit_determination

end module pod_emdac_runner_module