program test_force_model
    use pod_global, only: DP
    use pod_spice, only: spice_init, str2et  
    use pod_config, only: config, load_config, print_config, resolve_config_dependencies
    use pod_force_model_module, only: init_gravity_network, set_propagation_epoch
    use pod_integrator_module, only: adaptive_step_integrate, METHOD_RKF45, METHOD_RKF78
    
    implicit none
    
    real(DP) :: tdb_epoch
    real(DP), dimension(3) :: position, velocity
    
    ! 物理真实量纲变量
    real(DP), dimension(6) :: initial_state_real, final_state_45_real, final_state_78_real
    
    ! 无量纲化(积分器专用)变量
    real(DP), dimension(6) :: initial_state_nondim
    real(DP) :: t_start_nondim, t_end_nondim
    
    ! 统计与输出变量
    real(DP) :: cpu_start, cpu_end, time_45, time_78
    integer :: n_steps_45, n_steps_78
    real(DP), allocatable, dimension(:) :: times
    real(DP), allocatable, dimension(:,:) :: states
    
    ! 1. SPICE 初始化
    call spice_init()
    
    write(*,*) "=================================================="
    write(*,*) "   CAT POD System - 地月空间高精度积分测试 (无量纲化) "
    write(*,*) "=================================================="
    
    ! 2. 核心配置：加载或覆盖配置
    call load_config('dummy_test_config.txt')
    
    ! --- 物理模型配置 ---
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

    ! --- 极度关键：注入地月空间归一化单位 ---
    config%LU   = 384400.0_DP                  ! 长度单位 (km)
    config%TU   = 375189.032672856_DP          ! 时间单位 (s)
    config%MU   = 6.045626292270688E24_DP      ! 质量单位 (kg)
    config%VU   = 1.024550212739266_DP         ! 速度单位 (km/s)
    config%AccU = 2.730757367307436E-6_DP      ! 加速度单位 (km/s^2)

    ! --- 极度关键：注入分离后的 WRMS 容差 ---
    ! RKF45: 中阶求解器，容差相对较松
    config%rkf45_rel_tol = 1.0e-12_DP
    config%rkf45_abs_tol = 1.0e-12_DP
    ! RKF78: 高阶求解器，必须使用极其严格的容差激发其大步长潜力
    config%rkf78_rel_tol = 1.0e-12_DP
    config%rkf78_abs_tol = 1.0e-14_DP
    
    config%max_propagation_steps = 500000
    
    call resolve_config_dependencies()
    call init_gravity_network()
    call print_config()
    
    ! 3. 设定物理初始状态 (2024-03-09T12:00:00)
    call str2et('2024-03-09T12:00:00', tdb_epoch)
    write(*,*) ">>> 传播历元 (TDB) : ", tdb_epoch, " 秒"
    position = [100000.0_DP, 50000.0_DP, 20000.0_DP]  ! km
    velocity = [1.5_DP, 2.5_DP, 0.5_DP]               ! km/s
    initial_state_real = [position, velocity]
    
    ! =========================================================
    ! 4. 【核心架构操作】历元注入与状态去量纲化
    ! =========================================================
    
    ! A. 注入全局物理历元基准到力模型闭包中
    call set_propagation_epoch(tdb_epoch)
    
    ! B. 去量纲化：转换为积分器所需的纯数学相对状态
    initial_state_nondim(1:3) = initial_state_real(1:3) / config%LU
    initial_state_nondim(4:6) = initial_state_real(4:6) / config%VU
    
    t_start_nondim = 0.0_DP
    t_end_nondim   = 86400.0_DP / config%TU  ! 传播 1 天的无量纲时间
    
    ! =========================================================
    ! 5. 极简积分器调用测试 
    ! =========================================================
    
    write(*,*) ""
    write(*,*) ">>> 正在运行 RKF 4(5) 自适应积分..."
    call cpu_time(cpu_start)
    ! 积分器只看到 0.0 ~ 0.230(无量纲天) 以及微小的状态量
    call adaptive_step_integrate( &
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
    
    ! 还原物理量纲 (Dimensionalize)
    final_state_45_real(1:3) = states(n_steps_45, 1:3) * config%LU
    final_state_45_real(4:6) = states(n_steps_45, 4:6) * config%VU
    
    write(*,"(A, F8.4, A)") " [RKF45] 耗时       : ", time_45, " 秒"
    write(*,"(A, I8)")      " [RKF45] 总计算步数 : ", n_steps_45
    write(*,"(A, 3(E20.12))") " 最终位置 (km) : ", final_state_45_real(1:3)

    ! 【防溢出】必须释放内存，准备下一次积分测试！
    if (allocated(times)) deallocate(times)
    if (allocated(states)) deallocate(states)

    write(*,*) ""
    write(*,*) ">>> 正在运行 RKF 7(8) 自适应积分..."
    call cpu_time(cpu_start)
    call adaptive_step_integrate( &
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

    ! 还原物理量纲 (Dimensionalize)
    final_state_78_real(1:3) = states(n_steps_78, 1:3) * config%LU
    final_state_78_real(4:6) = states(n_steps_78, 4:6) * config%VU

    write(*,"(A, F8.4, A)") " [RKF78] 耗时       : ", time_78, " 秒"
    write(*,"(A, I8)")      " [RKF78] 总计算步数 : ", n_steps_78
    write(*,"(A, 3(E20.12))") " 最终位置 (km) : ", final_state_78_real(1:3)

    ! =========================================================
    ! 6. 结果交叉比对
    ! =========================================================
    write(*,*) ""
    write(*,*) "=================================================="
    write(*,*) "               积分结果一致性校验                 "
    write(*,*) "=================================================="
    write(*,"(A, 3(E20.12))") " RKF45 最终位置 (km) : ", final_state_45_real(1:3)
    write(*,"(A, 3(E20.12))") " RKF78 最终位置 (km) : ", final_state_78_real(1:3)
    write(*,"(A, 3(E20.12))") " 绝对位置差异 (km)   : ", abs(final_state_45_real(1:3) - final_state_78_real(1:3))
    
    ! 退出前养成良好习惯，清理内存
    if (allocated(times)) deallocate(times)
    if (allocated(states)) deallocate(states)
    
end program test_force_model