program test_cislunar_force_model
    use pod_global, only: DP
    use pod_spice, only: spice_init, str2et  
    use pod_force_model_module
    use pod_integrator_module, only: adaptive_step_integrate, METHOD_RKF45, METHOD_RKF78
    
    implicit none
    
    real(DP) :: tdb_epoch
    real(DP), dimension(3) :: position, velocity, acceleration
    
    ! --- 积分器相关变量 ---
    real(DP), dimension(6) :: initial_state, final_state_45, final_state_78
    real(DP) :: t_start, t_end, tolerance
    integer :: max_steps, n_steps_45, n_steps_78
    real(DP), allocatable, dimension(:) :: times
    real(DP), allocatable, dimension(:,:) :: states
    real(DP) :: cpu_start, cpu_end, time_45, time_78
    
    ! 1. 全局只允许在这里进行唯一的一次 SPICE 初始化！
    call spice_init()
    
    write(*,*) "=================================================="
    write(*,*) "   CAT POD System - 地月空间力学模型与积分测试    "
    write(*,*) "=================================================="
    
    ! 2. 任务级力模型配置 
    call set_perturbation_switches(.true., .true., .true., .false., .false., .false.)
    call set_gravity_degrees(10, 10)
    call set_active_planets([3, 10, 11])
    
    ! 3. 初始化多体网并打印面板
    call init_gravity_network()
    call print_force_model_status()
    
    ! 4. 设定测试初始状态
    call str2et('2024-03-09T12:00:00', tdb_epoch)
    
    position = [100000.0_DP, 50000.0_DP, 20000.0_DP]
    velocity = [1.5_DP, 2.5_DP, 0.5_DP]
    
    ! 组装 6D 状态向量
    initial_state(1:3) = position
    initial_state(4:6) = velocity
    
    write(*,*) ">>> 测试历元 (TDB 秒数) : ", tdb_epoch
    write(*,*) ">>> 卫星初始位置 (km): ", position
    write(*,*) ">>> 卫星初始速度 (km/s): ", velocity
    
    ! 5. 单步力模型基准测试
    write(*,*) ">>> 正在计算初始历元加速度..."
    call compute_acceleration(position, velocity, tdb_epoch, acceleration)
    write(*,"(A, 3(E20.12))") " 初始总加速度 (km/s^2): ", acceleration
    
    ! ==================================================
    ! 6. 积分器性能与精度测试
    ! ==================================================
    write(*,*) ""
    write(*,*) "=================================================="
    write(*,*) "               轨道积分器效能评估                 "
    write(*,*) "=================================================="
    
    ! 设置积分任务参数 (例如：向前推演 1 天)
    t_start = tdb_epoch
    t_end = tdb_epoch + 86400.0_DP 
    tolerance = 1.0D-7  ! 设置一个极其严格的容差
    max_steps = 500000
    
    write(*,*) "积分时长: 86400.0 秒 (1 天)"
    write(*,*) "局部截断误差容差: ", tolerance
    write(*,*) ""
    
    ! --- 测试 1: RKF45 ---
    write(*,*) ">>> 正在运行 RKF 4(5) 自适应积分..."
    call cpu_time(cpu_start)
    call adaptive_step_integrate(initial_state, t_start, t_end, max_steps, tolerance, &
                                 METHOD_RKF45, times, states, n_steps_45)
    call cpu_time(cpu_end)
    time_45 = cpu_end - cpu_start
    final_state_45 = states(n_steps_45, :)
    
    write(*,"(A, I8)")      " [RKF45] 总计算步数 : ", n_steps_45
    write(*,"(A, F8.4, A)") " [RKF45] 耗时       : ", time_45, " 秒"
    
    ! 释放内存，为下一次测试做准备
    deallocate(times, states)
    
    ! --- 测试 2: RKF78 ---
    write(*,*) ">>> 正在运行 RKF 7(8) 自适应积分..."
    call cpu_time(cpu_start)
    call adaptive_step_integrate(initial_state, t_start, t_end, max_steps, tolerance, &
                                 METHOD_RKF78, times, states, n_steps_78)
    call cpu_time(cpu_end)
    time_78 = cpu_end - cpu_start
    final_state_78 = states(n_steps_78, :)
    
    write(*,"(A, I8)")      " [RKF78] 总计算步数 : ", n_steps_78
    write(*,"(A, F8.4, A)") " [RKF78] 耗时       : ", time_78, " 秒"
    
    deallocate(times, states)
    
    ! --- 结果交叉对比 ---
    write(*,*) ""
    write(*,*) "=================================================="
    write(*,*) "               积分结果一致性校验                 "
    write(*,*) "=================================================="
    write(*,"(A, 3(E20.12))") " RKF45 最终位置 (km) : ", final_state_45(1:3)
    write(*,"(A, 3(E20.12))") " RKF78 最终位置 (km) : ", final_state_78(1:3)
    write(*,"(A, 3(E20.12))") " 位置差异 (km)       : ", abs(final_state_45(1:3) - final_state_78(1:3))
    
end program test_cislunar_force_model