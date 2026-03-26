program test_uq_api
    use pod_global, only: DP
    use pod_spice, only: spice_init
    use pod_force_model_module, only: set_perturbation_switches, set_gravity_degrees, set_active_planets, init_gravity_network
    
    ! 引入我们新改版的 API 和必需的常量开关
    use pod_uq_propagation, only: run_uq_propagation, METHOD_MC, METHOD_DA
    use pod_uq_base_module, only: INTEG_RKF78, INTEG_RK4
    use pod_uq_state_module, only: uq_state_type
    
    implicit none

    ! 定义测试所需的输入参数
    real(DP) :: nominal_orbit(6)
    real(DP) :: initial_covariance(6,6)
    real(DP) :: t0, tf, dt
    
    ! 声明用于接收结果的容器
    type(uq_state_type) :: initial_dist, final_dist_da, final_dist_mc

    ! ==========================================
    ! 0. 初始化物理环境
    ! ==========================================
    call spice_init()
    call set_perturbation_switches(.true., .true., .true., .false., .false., .false.)
    call set_gravity_degrees(4, 0)
    call set_active_planets([3, 10]) 
    call init_gravity_network()

    ! ==========================================
    ! 1. 准备初始状态 (硬编码或从文件读取)
    ! ==========================================
    nominal_orbit = [100000.0_DP, 50000.0_DP, 20000.0_DP, 1.5_DP, 2.5_DP, 0.5_DP]
    
    initial_covariance = 0.0_DP
    initial_covariance(1,1) = 1.0_DP**2
    initial_covariance(2,2) = 1.0_DP**2
    initial_covariance(3,3) = 1.0_DP**2
    initial_covariance(4,4) = 1.0D-3**2
    initial_covariance(5,5) = 1.0D-3**2
    initial_covariance(6,6) = 1.0D-3**2

    t0 = 0.0_DP
    tf = 3600.0_DP
    dt = 60.0_DP

    ! ==========================================
    ! 2. 静默调用 API 执行误差传播
    ! ==========================================
    write(*,*) '=== 测试 1: 使用 DA 算法传播 ==='
    call run_uq_propagation( &
        nominal_state = nominal_orbit, &
        initial_cov   = initial_covariance, &
        t_start       = t0, &
        t_end         = tf, &
        dt_init       = dt, &
        method_switch = METHOD_DA, &         ! <--- 开关：使用 DA
        integrator_switch = INTEG_RKF78, &   ! <--- 开关：使用 RKF78
        n_particles   = 10000, &             ! <--- 粒子数量
        save_results_to_file = .false., &    ! <--- 开关：不保存文件
        initial_state_out = initial_dist, &
        final_state_out   = final_dist_da &
    )

    write(*,*) '=== 测试 2: 使用 MC 算法对比 ==='
    call run_uq_propagation( &
        nominal_state = nominal_orbit, &
        initial_cov   = initial_covariance, &
        t_start       = t0, &
        t_end         = tf, &
        dt_init       = dt, &
        method_switch = METHOD_MC, &         ! <--- 开关：切换为 MC
        integrator_switch = INTEG_RKF78, &
        n_particles   = 10000, &
        save_results_to_file = .true., &     ! <--- 开关：保存文件以便分析
        initial_state_out = initial_dist, &
        final_state_out   = final_dist_mc &
    )

    ! 测试结束，可以在外层清理内存
    if (allocated(initial_dist%samples)) deallocate(initial_dist%samples)
    if (allocated(final_dist_da%samples)) deallocate(final_dist_da%samples)
    if (allocated(final_dist_mc%samples)) deallocate(final_dist_mc%samples)

end program test_uq_api