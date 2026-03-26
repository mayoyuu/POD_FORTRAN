program test_da_pipeline
    use pod_global, only: DP
    use pod_spice, only: spice_init, str2et
    use pod_dace_classes
    use pod_da_force_model_module
    use pod_da_integrator_module
    ! 如果你有底层的 dace 初始化库，确保引入。如果没有封装则忽略。
    ! use dace_core
    implicit none

    type(AlgebraicVector) :: state_da, deriv_da, state_next
    real(DP), dimension(6) :: r0
    real(DP) :: tdb_epoch
    integer :: i
     
    ! --- 积分器相关变量 ---
    type(AlgebraicVector) :: final_state_45, final_state_78
    real(DP) :: t_start, t_end, tolerance,dt
    integer :: max_steps, n_steps_45, n_steps_78
    real(DP), allocatable, dimension(:) :: times
    type(AlgebraicVector), allocatable, dimension(:) :: states
    real(DP) :: cpu_start, cpu_end

    real(DP) :: eval_inputs(6)
    real(DP), allocatable :: eval_results(:)


    write(*,*) '=================================================='
    write(*,*) '        🚀 POD DA 动力学链路全量集成测试        '
    write(*,*) '=================================================='

    ! =========================================================
    ! 阶段 0: 环境与模型初始化
    ! =========================================================
    ! 1. 全局只允许在这里进行唯一的一次 SPICE 初始化！
    call spice_init()
    write(*,*) '[INFO] 正在配置力学模型参数...'
    ! 开启地球高阶引力场(比如4阶)，关闭其他摄动以保持测试纯粹
    call set_perturbation_switches(.true., .true., .true., .false., .false., .false.)
    call set_gravity_degrees(4, 0) 
    call set_active_planets([3,10,11]) ! 只激活地球 (ID=3)

    ! 3. 初始化多体网并打印面板
    call init_gravity_network()
    call print_force_model_status()
    
    ! 4. 设定测试初始状态
    call str2et('2024-03-09T12:00:00', tdb_epoch)
    
    ! 打印我们写好的监控面板
    call print_force_model_status()

    ! =========================================================
    ! 阶段 1: DA 状态向量初始化 (注入空间不确定性)
    ! =========================================================
    write(*,*) '[INFO] 正在初始化 DA 独立变量...'
    
    ! 设定一个典型的 LEO 轨道初值 (位置 km, 速度 km/s)
    r0 = [100000.0_DP, 50000.0_DP, 20000.0_DP, 1.5_DP, 2.5_DP, 0.5_DP]
    call dace_initialize(2, 6) ! 初始化 DACE 引擎，阶数=2（位置和速度），变量数=6（状态向量维度）

    call state_da%init(6)
    do i = 1, 6
        state_da%elements(i) = r0(i) + da_var(i) ! 这里 da_var(i) 是一个 DA 变量，代表第 i 个状态分量的微小偏移
    end do
    
    write(*,*) '初始状态: ', state_da%elements(1), state_da%elements(2), state_da%elements(3), &
                     state_da%elements(4), state_da%elements(5), state_da%elements(6)

    ! =========================================================
    ! 阶段 2: 测试 Force Model 与右函数
    ! =========================================================
    write(*,*) '[INFO] 正在调用 da_compute_derivatives...'
    ! 设置积分任务参数 (例如：向前推演 1 天)
    t_start = tdb_epoch
    t_end = tdb_epoch + 86400.0_DP 
    tolerance = 1.0D-7  ! 设置一个极其严格的容差
    max_steps = 500000
    dt = 60.0_DP  ! RK4 固定步长测试用

    call deriv_da%init(6)
    ! 这将穿透你的 compute_acceleration，并调用引力场里的 f_zonal_da 和 f_tesseral_da
    call da_compute_derivatives(state_da, t_start, deriv_da)
    
    write(*,*) '计算成功！当前状态导数 (速度与加速度):'
    write(*,*) '常数项: ', deriv_da%cons()
    
    ! 验证偏导数是否成功传递：X方向的加速度，对初始X位置应该有非零的偏导数
    write(*,*) 'd(Acc_x) / d(Pos_x) = ', deriv_da%elements(4)%deriv(1)

    ! =========================================================
    ! 阶段 3: 测试 RK4 积分器与状态转移矩阵 (STM) 提取
    ! =========================================================
    write(*,*) '[INFO] 正在进行 RK4 单步积分 (dt = ', dt, ' s)...'
    call state_next%init(6)
    
    ! 自动调用 RK4 算法，并测量 CPU 时间
    call cpu_time(cpu_start)
    call da_rk4_integrate(state_da, dt, t_start, state_next)
    call cpu_time(cpu_end)
    write(*,*) 'RK4 单步积分完成！耗时: ', cpu_end - cpu_start, ' 秒'
    write(*,*) '60秒后的状态 (常数项): ', state_next%cons()
    write(*,*) '--------------------------------------------------'
    write(*,*) '✨ 状态转移矩阵 (STM) 的左上角 3x3 ✨'
    write(*,*) '--------------------------------------------------'
   do i = 1, 3
        write(*, "(3(ES15.6, 2X))") &
            state_next%elements(i)%get_deriv_value(1), &
            state_next%elements(i)%get_deriv_value(2), &
            state_next%elements(i)%get_deriv_value(3)
    end do
    write(*,*) '=================================================='


    ! 调用 RKF45 变步长积分器进行更长时间的积分，并提取 STM
    write(*,*) '[INFO] 正在进行 RKF45 变步长积分...'

    call cpu_time(cpu_start)
    call da_adaptive_step_integrate(state_da, t_start, t_end, max_steps, tolerance, &
                                 METHOD_RKF45, times, states, n_steps_45)
    call cpu_time(cpu_end)
    final_state_45 = states(n_steps_45)
    write(*,*) 'RKF45 积分完成！耗时: ', cpu_end - cpu_start, ' 秒'
    write(*,*) '最终的状态: ', final_state_45%cons()
    write(*,*) '--------------------------------------------------'
    write(*,*) '✨ RKF45 最终状态的 STM 左上角 3x3 ✨'
    write(*,*) '--------------------------------------------------'
    do i = 1, 3
        write(*, "(3(ES15.6, 2X))") &
            final_state_45%elements(i)%get_deriv_value(1), &
            final_state_45%elements(i)%get_deriv_value(2), &
            final_state_45%elements(i)%get_deriv_value(3)
    end do
    write(*,*) '=================================================='

    ! 调用 RKF78 变步长积分器进行更长时间的积分，并提取 STM
    write(*,*) '[INFO] 正在进行 RKF78 变步长积分...'
    call cpu_time(cpu_start)
    call da_adaptive_step_integrate(state_da, t_start, t_end, max_steps, tolerance, &
                                 METHOD_RKF78, times, states, n_steps_78)
    call cpu_time(cpu_end)
    write(*,*) 'RKF78 积分完成！耗时: ', cpu_end - cpu_start, ' 秒'
    final_state_78 = states(n_steps_78)
    write(*,*) '最终的状态: ', final_state_78%cons()
    write(*,*) '--------------------------------------------------'
    ! write(*,*) '✨ RKF78 最终状态 (DA) ✨'
    ! write(*,*) '--------------------------------------------------'
    ! call final_state_78%print()
    write(*,*) '--------------------------------------------------'
    write(*,*) '✨ RKF78 最终状态的 STM 左上角 3x3 ✨'
    write(*,*) '--------------------------------------------------'
    do i = 1, 3
        write(*, "(3(ES15.6, 2X))") &
            final_state_78%elements(i)%get_deriv_value(1), &
            final_state_78%elements(i)%get_deriv_value(2), &
            final_state_78%elements(i)%get_deriv_value(3)
    end do
    write(*,*) '=================================================='

    write(*,*) '🎉 求值.eval测试开始 🎉'
    eval_inputs = [0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP]
    eval_results = final_state_78%eval(eval_inputs) 
    ! 4. 打印结果
    write(*,*) 'RKF78 最终状态求值结果 (数值): ', eval_results
end program test_da_pipeline