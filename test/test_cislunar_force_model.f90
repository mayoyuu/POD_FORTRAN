program test_cislunar_force_model
    use pod_global, only: DP
    use pod_spice, only: spice_init, str2et  
    use pod_config, only: config, load_config  ! <-- 引入全局配置
    use pod_force_model_module
    use pod_integrator_module, only: adaptive_step_integrate, METHOD_RKF45, METHOD_RKF78
    
    implicit none
    
    real(DP) :: tdb_epoch
    real(DP), dimension(3) :: position, velocity, acceleration
    
    real(DP), dimension(6) :: initial_state, final_state_45, final_state_78
    real(DP) :: t_start, t_end, cpu_start, cpu_end, time_45, time_78
    integer :: n_steps_45, n_steps_78
    real(DP), allocatable, dimension(:) :: times
    real(DP), allocatable, dimension(:,:) :: states
    
    ! 1. SPICE 初始化
    call spice_init()
    
    write(*,*) "=================================================="
    write(*,*) "   CAT POD System - 地月空间力学模型与积分测试    "
    write(*,*) "=================================================="
    
    ! 2. 核心：通过修改 config 对象来配置这次任务的物理环境
    ! (如果未来写好了配置文件，这里只需一句 call load_config('test.ini') 即可)
    config%use_earth_nspheric = .true.
    config%use_moon_nspheric  = .true.
    config%use_third_body     = .true.
    config%use_srp            = .false.
    config%earth_degree       = 10
    config%moon_degree        = 10
    
    ! 激活 地球(3), 月球(10), 太阳(11)
    config%use_planet = .false. 
    config%use_planet(3)  = .true.
    config%use_planet(10) = .true.
    config%use_planet(11) = .true.

    ! 覆盖积分器默认配置 (仅供本次严格测试使用)
    config%rkf45_tolerance = 1.0D-7
    config%max_propagation_steps = 500000
    
    ! 3. 初始化物理引擎 (它会自己去读上面的 config)
    call init_gravity_network()
    call print_force_model_status()
    
    ! 4. 设定初始状态
    call str2et('2024-03-09T12:00:00', tdb_epoch)
    position = [100000.0_DP, 50000.0_DP, 20000.0_DP]
    velocity = [1.5_DP, 2.5_DP, 0.5_DP]
    initial_state = [position, velocity]
    
    t_start = tdb_epoch
    t_end = tdb_epoch + 86400.0_DP 
    
    ! 5. 极简的积分器调用：
    ! 参数极其干净！积分器内部会自动从 config 中提取容差、最小步长等
    write(*,*) ">>> 正在运行 RKF 4(5) 自适应积分..."
    call cpu_time(cpu_start)
    call adaptive_step_integrate( &
        state = initial_state, &
        t_start = t_start, &
        t_end = t_end, &
        integrator_method = METHOD_RKF45, &
        times = times, &
        states = states, &
        n_steps = n_steps_45 &
    )
    call cpu_time(cpu_end)
    time_45 = cpu_end - cpu_start
    
    final_state_45 = states(n_steps_45, :)
    write(*,"(A, F8.4, A)") " [RKF45] 耗时       : ", time_45, " 秒"
    write(*,"(A, I8)")      " [RKF45] 总计算步数 : ", n_steps_45
    write(*,"(A, 3(E20.12))") " 最终位置 (km) : ", final_state_45(1:3)

    write(*,*) "?>> 正在运行 RKF 7(8) 自适应积分..."
    call cpu_time(cpu_start)
    call adaptive_step_integrate( &
        state = initial_state, &
        t_start = t_start, &
        t_end = t_end, &
        integrator_method = METHOD_RKF78, &
        times = times, &
        states = states, &
        n_steps = n_steps_78 &
    )
    call cpu_time(cpu_end)
    time_78 = cpu_end - cpu_start

    final_state_78 = states(n_steps_78, :)
    write(*,"(A, F8.4, A)") " [RKF78] 耗时       : ", time_78, " 秒"
    write(*,"(A, I8)")      " [RKF78] 总计算步数 : ", n_steps_78
    write(*,"(A, 3(E20.12))") " 最终位置 (km) : ", final_state_78(1:3)

        ! 结果对比
    write(*,*) ""
    write(*,*) "=================================================="
    write(*,*) "               积分结果一致性校验                 "
    write(*,*) "=================================================="
    write(*,"(A, 3(E20.12))") " RKF45 最终位置 (km) : ", final_state_45(1:3)
    write(*,"(A, 3(E20.12))") " RKF78 最终位置 (km) : ", final_state_78(1:3)
    write(*,"(A, 3(E20.12))") " 位置差异 (km)       : ", abs(final_state_45(1:3) - final_state_78(1:3))
    
end program test_cislunar_force_model

