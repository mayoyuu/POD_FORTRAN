program test_uq_api
    use pod_global, only: DP
    use pod_engine_module, only: pod_engine_init
    use pod_spice, only: str2et
    use pod_force_model_module
    use pod_dace_classes
    
    ! 引入我们新改版的 API 和必需的常量开关
    use pod_uq_propagation, only: run_uq_propagation, METHOD_MC, METHOD_DA
    use pod_uq_state_module, only: uq_state_type
    
    implicit none

    ! 所需参数
    real(DP) :: nominal_orbit(6)
    real(DP) :: initial_covariance(6,6)
    real(DP) :: t0, tf, dt, epoch_start
    type(uq_state_type) :: initial_dist, final_dist_da, final_dist_mc

    ! ==========================================
    ! 0. 初始化物理环境
    ! ==========================================
    call pod_engine_init('dummy_test_config.txt')

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
    tf = 86400.0_DP
    dt = 60.0_DP
    call str2et('2024-03-09T12:00:00', epoch_start) 

   ! ==========================================
    ! 2. 极致清爽的 API 调用
    ! ==========================================
    call dace_initialize(2, 6)
    write(*,*) '=== 测试 1: 使用 DA 算法传播 ==='
    call run_uq_propagation( &
        nominal_state = nominal_orbit, &
        initial_cov   = initial_covariance, &
        epoch0        = epoch_start, &       ! <--- 传入物理基准历元
        t_start       = t0, &
        t_end         = tf, &
        method_switch = METHOD_DA, &         ! <--- 只需告诉 API：用 DA 算法
        n_particles   = 10000, &
        save_results_to_file = .false., &
        initial_state_out = initial_dist, &
        final_state_out   = final_dist_da, &
        file_prefix   = './output/da_run_orbit_1'&! <--- 传入文件前缀 (可选参数)
    )
    ! 注意：integrator_switch 参数彻底消失了！API 会自动用最高精度的 RKF78 处理。

    write(*,*) '=== 测试 2: 使用 MC 算法对比 ==='
    call run_uq_propagation( &
        nominal_state = nominal_orbit, &
        initial_cov   = initial_covariance, &
        epoch0        = epoch_start, &
        t_start       = t0, &
        t_end         = tf, &
        method_switch = METHOD_MC, &         ! <--- 只需告诉 API：切换为 MC 算法
        n_particles   = 10000, &
        save_results_to_file = .true., &
        initial_state_out = initial_dist, &
        final_state_out   = final_dist_mc, &
        file_prefix   = './output/mc_run_orbit_1'&! <--- 传入文件前缀 (可选参数)
    )

    ! 清理内存
    if (allocated(initial_dist%samples)) deallocate(initial_dist%samples)
    if (allocated(final_dist_da%samples)) deallocate(final_dist_da%samples)
    if (allocated(final_dist_mc%samples)) deallocate(final_dist_mc%samples)

end program test_uq_api