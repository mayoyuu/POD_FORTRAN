program test_da_force_model
    use pod_global, only: DP
    use pod_spice, only: spice_init, str2et
    use pod_config, only: config, load_config, print_config, resolve_config_dependencies
    
    use pod_da_force_model_module, only: da_compute_acceleration,init_gravity_network, set_propagation_epoch
    use pod_da_integrator_module, only: da_adaptive_step_integrate, METHOD_RKF45, METHOD_RKF78
    use pod_dace_classes
    
    implicit none

    real(DP) :: tdb_epoch
    real(DP), dimension(6) :: initial_state_real_dp
    
    ! DA 状态变量 (物理真实量纲)
    type(AlgebraicVector) :: final_state_45_real, final_state_78_real
    
    ! DA 状态变量 (无量纲化积分器专用)
    type(AlgebraicVector) :: initial_state_nondim
    real(DP) :: t_start_nondim, t_end_nondim
    
    ! 积分器控制与统计
    real(DP) :: cpu_start, cpu_end, time_45, time_78
    integer :: n_steps_45, n_steps_78
    real(DP), allocatable, dimension(:) :: times
    type(AlgebraicVector), allocatable, dimension(:) :: states
    integer :: i
    
    ! Eval 测试变量
    real(DP) :: eval_inputs(6)
    real(DP), allocatable :: eval_results(:)

    write(*,*) '=================================================='
    write(*,*) '      🚀 POD DA 动力学链路全量集成测试 (无量纲版)      '
    write(*,*) '=================================================='

    ! =========================================================
    ! 阶段 0: 环境与物理模型配置 (与 DP 测试完全对齐)
    ! =========================================================
    call spice_init()
    
    call load_config('dummy_test_config.txt')
    
    ! --- 物理模型摄动配置 ---
    config%use_earth_nspheric = .true.
    config%use_moon_nspheric  = .true.
    config%use_third_body     = .true.
    config%use_srp            = .false.
    config%earth_degree       = 10
    config%moon_degree        = 10
    config%use_planet         = .false.
    config%use_planet(3)      = .true.
    config%use_planet(10)     = .true.
    config%use_planet(11)     = .true.

    ! --- 注入地月空间归一化单位 ---
    config%LU   = 384400.0_DP                  ! km
    config%TU   = 375189.032672856_DP          ! s
    config%MU   = 6.045626292270688E24_DP      ! kg
    config%VU   = 1.024550212739266_DP         ! km/s
    config%AccU = 2.730757367307436E-6_DP      ! km/s^2

    ! --- 注入 WRMS 容差 ---
    config%rkf45_rel_tol = 1.0e-12_DP
    config%rkf45_abs_tol = 1.0e-12_DP
    config%rkf78_rel_tol = 1.0e-12_DP
    config%rkf78_abs_tol = 1.0e-14_DP
    config%max_propagation_steps = 500000
    
    call resolve_config_dependencies()
    call init_gravity_network()
    call print_config()

    ! =========================================================
    ! 阶段 1: 历元注入与 DA 状态向量初始化 (空间不确定性去量纲化)
    ! =========================================================
    write(*,*) '[INFO] 正在初始化 DA 独立变量并去量纲化...'
    
    call str2et('2024-03-09T12:00:00', tdb_epoch)
    call set_propagation_epoch(tdb_epoch)
    
    initial_state_real_dp = [100000.0_DP, 50000.0_DP, 20000.0_DP, 1.5_DP, 2.5_DP, 0.5_DP]
    
    ! 初始化 DACE 引擎 (2阶，6个变量)
    call dace_initialize(2, 6)
    
    call initial_state_nondim%init(6)
    
    ! 【核心魔法】赋予 DA 物理误差，并在底层重载中自动完成整体无量纲化
    do i = 1, 3
        ! 位置无量纲化: (Pos_real + d_Pos_real) / LU
        initial_state_nondim%elements(i) = (initial_state_real_dp(i) + da_var(i)) / config%LU
        ! 速度无量纲化: (Vel_real + d_Vel_real) / VU
        initial_state_nondim%elements(i+3) = (initial_state_real_dp(i+3) + da_var(i+3)) / config%VU
    end do
    
    ! 时间无量纲化
    t_start_nondim = 0.0_DP
    t_end_nondim   = 86400.0_DP / config%TU 
    
    ! =========================================================
    ! 阶段 2: 测试 RKF45 变步长积分器与 STM 提取
    ! =========================================================
    write(*,*) '[INFO] 正在进行 RKF45 变步长积分 (DA 版)...'

    call cpu_time(cpu_start)
    ! 注意：调用的是我上一版为你更新过接口参数的 da_adaptive_step_integrate
    call da_adaptive_step_integrate( &
        state = initial_state_nondim, &
        t_start = t_start_nondim, &
        t_end = t_end_nondim, &
        integrator_method = METHOD_RKF45, &
        times = times, &
        states = states, &
        n_steps = n_steps_45 &
    )
    call cpu_time(cpu_end)
    time_45 = cpu_end - cpu_start
    
    ! 还原物理量纲 (重新乘回 LU 和 VU)
    call final_state_45_real%init(6)
    do i = 1, 3
        final_state_45_real%elements(i) = states(n_steps_45)%elements(i) * config%LU
        final_state_45_real%elements(i+3) = states(n_steps_45)%elements(i+3) * config%VU
    end do
    
    write(*,"(A, F8.4, A)") " [RKF45-DA] 耗时       : ", time_45, " 秒"
    write(*,"(A, I8)")      " [RKF45-DA] 总计算步数 : ", n_steps_45
    write(*,"(A, 3(E20.12))") " 最终位置速度常数项 (km) : ", final_state_45_real%cons()
    
    write(*,*) '--------------------------------------------------'
    write(*,*) '✨ RKF45 最终物理状态的 STM 左上角 3x3 ✨'
    write(*,*) '--------------------------------------------------'
    do i = 1, 3
        write(*, "(3(ES15.6, 2X))") &
            final_state_45_real%elements(i)%get_deriv_value(1), &
            final_state_45_real%elements(i)%get_deriv_value(2), &
            final_state_45_real%elements(i)%get_deriv_value(3)
    end do
    
    ! 必须释放内存以供 RKF78 使用
    if (allocated(times)) deallocate(times)
    if (allocated(states)) deallocate(states)

    ! =========================================================
    ! 阶段 3: 测试 RKF78 变步长积分器与 STM 提取
    ! =========================================================
    write(*,*) ''
    write(*,*) '[INFO] 正在进行 RKF78 变步长积分 (DA 版)...'
    
    call cpu_time(cpu_start)
    call da_adaptive_step_integrate( &
        state = initial_state_nondim, &
        t_start = t_start_nondim, &
        t_end = t_end_nondim, &
        integrator_method = METHOD_RKF78, &
        times = times, &
        states = states, &
        n_steps = n_steps_78 &
    )
    call cpu_time(cpu_end)
    time_78 = cpu_end - cpu_start
    
    ! 还原物理量纲
    call final_state_78_real%init(6)
    do i = 1, 3
        final_state_78_real%elements(i) = states(n_steps_78)%elements(i) * config%LU
        final_state_78_real%elements(i+3) = states(n_steps_78)%elements(i+3) * config%VU
    end do
    
    write(*,"(A, F8.4, A)") " [RKF78-DA] 耗时       : ", time_78, " 秒"
    write(*,"(A, I8)")      " [RKF78-DA] 总计算步数 : ", n_steps_78
    write(*,"(A, 3(E20.12))") " 最终位置速度常数项 (km) : ", final_state_78_real%cons()
    
    write(*,*) '--------------------------------------------------'
    write(*,*) '✨ RKF78 最终物理状态的 STM 左上角 3x3 ✨'
    write(*,*) '--------------------------------------------------'
    do i = 1, 3
        write(*, "(3(ES15.6, 2X))") &
            final_state_78_real%elements(i)%get_deriv_value(1), &
            final_state_78_real%elements(i)%get_deriv_value(2), &
            final_state_78_real%elements(i)%get_deriv_value(3)
    end do

    ! =========================================================
    ! 阶段 4: 结果交叉比对与求值测试
    ! =========================================================
    write(*,*) ""
    write(*,*) "=================================================="
    write(*,*) "               DA 积分结果一致性校验              "
    write(*,*) "=================================================="
    write(*,"(A, 3(E20.12))") " 绝对结果差异 (km)   : ", abs(final_state_45_real%cons() - final_state_78_real%cons())

    write(*,*) '🎉 求值 .eval 测试开始 🎉'
    eval_inputs = [0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP]
    ! 测试: 对物理化之后的 78 阶 DA 向量求值，应与 cons() 一致
    eval_results = final_state_78_real%eval(eval_inputs) 
    write(*,*) 'RKF78_real 零点求值结果 (数值): ', eval_results(1:3)
    
    if (allocated(times)) deallocate(times)
    if (allocated(states)) deallocate(states)
    
end program test_da_force_model