module pod_integrator_module
    use pod_global, only: DP
    use pod_config, only: config
    use pod_force_model_module, only: compute_acceleration
    
    implicit none

    ! --- 积分器类型常量定义 (暴露给外部用户使用) ---
    integer, parameter, public :: METHOD_RKF45 = 1
    integer, parameter, public :: METHOD_RKF78 = 2
    
contains

    subroutine compute_derivatives(state, time, derivatives)
        real(DP), dimension(6), intent(in) :: state
        real(DP), intent(in) :: time
        real(DP), dimension(6), intent(out) :: derivatives
        
        real(DP), dimension(3) :: position, velocity, acceleration
        
        ! 提取位置和速度
        position = state(1:3)
        velocity = state(4:6)
        
        ! 计算加速度
        call compute_acceleration(position, velocity, time, acceleration)
        
        ! 构建导数向量 [dx/dt, dy/dt, dz/dt, dvx/dt, dvy/dt, dvz/dt]
        derivatives(1:3) = velocity
        derivatives(4:6) = acceleration
    end subroutine compute_derivatives
    

    subroutine rk4_integrate(state, dt, time, new_state)
        real(DP), dimension(6), intent(in) :: state
        real(DP), intent(in) :: dt, time
        real(DP), dimension(6), intent(out) :: new_state
        
        real(DP), dimension(6) :: k1, k2, k3, k4
        real(DP), dimension(6) :: temp_state
        
        ! RK4 积分方法
        ! k1 = f(t, y)
        call compute_derivatives(state, time, k1)
        
        ! k2 = f(t + dt/2, y + dt*k1/2)
        temp_state = state + 0.5_DP * dt * k1
        call compute_derivatives(temp_state, time + 0.5_DP * dt, k2)
        
        ! k3 = f(t + dt/2, y + dt*k2/2)
        temp_state = state + 0.5_DP * dt * k2
        call compute_derivatives(temp_state, time + 0.5_DP * dt, k3)
        
        ! k4 = f(t + dt, y + dt*k3)
        temp_state = state + dt * k3
        call compute_derivatives(temp_state, time + dt, k4)
        
        ! y(t + dt) = y(t) + dt*(k1 + 2*k2 + 2*k3 + k4)/6
        new_state = state + dt * (k1 + 2.0_DP * k2 + 2.0_DP * k3 + k4) / 6.0_DP
    end subroutine rk4_integrate
    
    subroutine rkf45_integrate(state, dt, time, new_state)
        real(DP), dimension(6), intent(in) :: state
        real(DP), intent(in) :: dt, time
        real(DP), dimension(6), intent(out) :: new_state
        
        real(DP), dimension(6) :: temp_4th, temp_5th
        real(DP) :: err_est
        
        ! 1. 直接调用核心计算引擎
        call rkf45_step(state, dt, time, temp_4th, temp_5th, err_est)
        
        ! 2. 丢弃4阶解和误差，只把最高精度的5阶解作为结果返回
        new_state = temp_5th
    end subroutine rkf45_integrate

    subroutine rkf78_integrate(state, dt, time, new_state)
        real(DP), dimension(6), intent(in) :: state
        real(DP), intent(in) :: dt, time
        real(DP), dimension(6), intent(out) :: new_state
        
        real(DP), dimension(6) :: state_7th, state_8th
        real(DP) :: error_estimate, tolerance
        
        tolerance = config%propagation_tolerance
        
        call rkf78_step(state, dt, time, state_7th, state_8th, error_estimate)
        
        ! 返回高阶解
        new_state = state_8th
        
        ! 简单的误差检查逻辑
        if (error_estimate > tolerance) then
            ! 这里可以添加警告，但在自适应循环中，这个逻辑会被上层接管
        end if
    end subroutine rkf78_integrate
    
    subroutine adaptive_step_integrate(state, t_start, t_end, integrator_method, &
                                       times, states, n_steps, &
                                       max_steps_in, tolerance_in, dt_min_in, dt_max_in)
        ! 核心物理参数
        real(DP), dimension(6), intent(in) :: state
        real(DP), intent(in) :: t_start, t_end
        integer, intent(in) :: integrator_method
        
        ! 输出参数
        real(DP), allocatable, dimension(:), intent(out) :: times
        real(DP), allocatable, dimension(:,:), intent(out) :: states
        integer, intent(out) :: n_steps
        
        ! 可选的控制参数 (如果外部不传，就用 config 里的默认值)
        integer, intent(in), optional :: max_steps_in
        real(DP), intent(in), optional :: tolerance_in, dt_min_in, dt_max_in
        
        ! 内部工作变量
        integer :: max_steps
        real(DP) :: tolerance, dt_min, dt_max, safety_factor
        real(DP) :: current_time, dt
        real(DP), dimension(6) :: current_state
        real(DP), dimension(6) :: next_state_4th, next_state_5th, next_state_7th, next_state_8th
        real(DP), dimension(6) :: next_state_high ! 用于统一暂存不同算法的高阶输出
        real(DP) :: error_estimate,  exp_power
        integer :: i
        
       ! ==========================================
        ! 1. 智能参数绑定 (结合 optional 与全局 config)
        ! ==========================================
        max_steps = config%max_propagation_steps
        if (present(max_steps_in)) max_steps = max_steps_in
        
        ! 根据不同的积分器，从 config 中提取专属的默认参数
        if (integrator_method == METHOD_RKF45) then
            tolerance = config%rkf45_tolerance
            dt_min = config%rkf45_min_step
            dt_max = config%rkf45_max_step
        else if (integrator_method == METHOD_RKF78) then
            tolerance = config%rkf78_tolerance
            dt_min = config%rkf78_min_step
            dt_max = config%rkf78_max_step
        else
            print *, "Error: Unsupported integrator method!"
            stop
        end if
        
        ! 如果用户显式传入了参数，则覆盖 config 的默认设定
        if (present(tolerance_in)) tolerance = tolerance_in
        if (present(dt_min_in)) dt_min = dt_min_in
        if (present(dt_max_in)) dt_max = dt_max_in
        
        safety_factor = 0.9_DP ! 这个也可以未来加进 config 里
        
        ! ==========================================
        ! 2. 内存由积分器全权管理
        ! ==========================================
        ! 确保外部传入的 times 和 states 是未分配状态，积分器自己负责分配
        if (allocated(times)) deallocate(times)
        if (allocated(states)) deallocate(states)
        allocate(times(max_steps))
        allocate(states(max_steps, 6))
        
        ! 初始状态
        current_state = state
        current_time = t_start
        times(1) = current_time
        states(1, :) = current_state
        
        ! 初始步长
        dt = min(dt_max, (t_end - t_start) / 100.0_DP)
        
        n_steps = 1
        
        do i = 2, max_steps
            ! 检查是否到达结束时间
            if (current_time >= t_end) exit

           ! ==========================================
            ! 1. 积分器路由分发 (使用易读的常量)
            ! ==========================================
            if (integrator_method == METHOD_RKF78) then
                call rkf78_step(current_state, dt, current_time, next_state_7th, next_state_8th, error_estimate)
                next_state_high = next_state_8th
                exp_power = 0.125_DP              ! 1/8
            else if (integrator_method == METHOD_RKF45) then
                call rkf45_step(current_state, dt, current_time, next_state_4th, next_state_5th, error_estimate)
                next_state_high = next_state_5th
                exp_power = 0.25_DP               ! 1/4
            else
                ! 防御性编程：用户传了不支持的常量
                print *, "Error: Unsupported integrator method!"
                stop
            end if

            ! =========================================================
            ! 2. 检查误差并决定是否推进
            ! =========================================================
            if (error_estimate <= tolerance) then
                ! 误差达标：正式接受这一步
                current_state = next_state_high
                current_time = current_time + dt
                n_steps = n_steps + 1
                
                times(n_steps) = current_time
                states(n_steps, :) = current_state
                
                ! 调整（尝试放大）步长，使用动态的 exp_power
                if (error_estimate > 0.0_DP) then
                    dt = safety_factor * dt * (tolerance / error_estimate)**exp_power
                end if
                dt = max(dt_min, min(dt_max, dt))
                
            else
                ! 误差超标：拒绝这一步（注意：不更新 current_state 和 current_time）
                
                ! 检查是否已经触底
                if (dt <= dt_min) then
                    ! 【防抱死补丁】已经缩小到极小步长依然无法满足容差，强行接受步长并前进！
                    current_state = next_state_high
                    current_time = current_time + dt
                    n_steps = n_steps + 1
                    
                    times(n_steps) = current_time
                    states(n_steps, :) = current_state
                    
                    ! 保持 dt 为 dt_min 继续尝试下一步
                else
                    ! 正常拒绝这一步，减小步长，重新计算当前时刻。使用动态的 exp_power
                    dt = safety_factor * dt * (tolerance / error_estimate)**exp_power
                    dt = max(dt_min, dt)
                end if
            end if
            
            ! 确保不超过结束时间
            if (current_time + dt > t_end) then
                dt = t_end - current_time
            end if
        end do
        
       ! 调整数组大小 
        if (n_steps < max_steps) then
            block
                real(DP), allocatable, dimension(:) :: temp_times
                real(DP), allocatable, dimension(:,:) :: temp_states
                
                allocate(temp_times(n_steps))
                allocate(temp_states(n_steps, 6))
                
                temp_times = times(1:n_steps)
                temp_states = states(1:n_steps, :)
                
                call move_alloc(temp_times, times)
                call move_alloc(temp_states, states)
            end block
        end if
    end subroutine adaptive_step_integrate

    subroutine rkf45_step(state, dt, time, state_4th, state_5th, error_estimate)
        real(DP), dimension(6), intent(in) :: state
        real(DP), intent(in) :: dt, time
        real(DP), dimension(6), intent(out) :: state_4th, state_5th
        real(DP), intent(out) :: error_estimate
        
        real(DP), dimension(6) :: k1, k2, k3, k4, k5, k6
        real(DP), dimension(6) :: temp_state
        
        ! RKF45 系数
        real(DP), parameter :: a2 = 1.0_DP / 5.0_DP
        real(DP), parameter :: a3 = 3.0_DP / 10.0_DP
        real(DP), parameter :: a4 = 4.0_DP / 5.0_DP
        real(DP), parameter :: a5 = 8.0_DP / 9.0_DP
        
        real(DP), parameter :: b21 = 1.0_DP / 5.0_DP
        real(DP), parameter :: b31 = 3.0_DP / 40.0_DP
        real(DP), parameter :: b32 = 9.0_DP / 40.0_DP
        real(DP), parameter :: b41 = 44.0_DP / 45.0_DP
        real(DP), parameter :: b42 = -56.0_DP / 15.0_DP
        real(DP), parameter :: b43 = 32.0_DP / 9.0_DP
        real(DP), parameter :: b51 = 19372.0_DP / 6561.0_DP
        real(DP), parameter :: b52 = -25360.0_DP / 2187.0_DP
        real(DP), parameter :: b53 = 64448.0_DP / 6561.0_DP
        real(DP), parameter :: b54 = -212.0_DP / 729.0_DP
        real(DP), parameter :: b61 = 9017.0_DP / 3168.0_DP
        real(DP), parameter :: b62 = -355.0_DP / 33.0_DP
        real(DP), parameter :: b63 = 46732.0_DP / 5247.0_DP
        real(DP), parameter :: b64 = 49.0_DP / 176.0_DP
        real(DP), parameter :: b65 = -5103.0_DP / 18656.0_DP
        
        real(DP), parameter :: c1 = 35.0_DP / 384.0_DP
        real(DP), parameter :: c3 = 500.0_DP / 1113.0_DP
        real(DP), parameter :: c4 = 125.0_DP / 192.0_DP
        real(DP), parameter :: c5 = -2187.0_DP / 6784.0_DP
        real(DP), parameter :: c6 = 11.0_DP / 84.0_DP
        
        real(DP), parameter :: d1 = 5179.0_DP / 57600.0_DP
        real(DP), parameter :: d3 = 7571.0_DP / 16695.0_DP
        real(DP), parameter :: d4 = 393.0_DP / 640.0_DP
        real(DP), parameter :: d5 = -92097.0_DP / 339200.0_DP
        real(DP), parameter :: d6 = 187.0_DP / 2100.0_DP
        real(DP), parameter :: d7 = 1.0_DP / 40.0_DP
        
        ! 计算k1到k6
        call compute_derivatives(state, time, k1)
        
        temp_state = state + dt * b21 * k1
        call compute_derivatives(temp_state, time + dt * b21, k2)
        
        temp_state = state + dt * (b31 * k1 + b32 * k2)
        call compute_derivatives(temp_state, time + dt * (b31 + b32), k3)
        
        temp_state = state + dt * (b41 * k1 + b42 * k2 + b43 * k3)
        call compute_derivatives(temp_state, time + dt * (b41 + b42 + b43), k4)
        
        temp_state = state + dt * (b51 * k1 + b52 * k2 + b53 * k3 + b54 * k4)
        call compute_derivatives(temp_state, time + dt * (b51 + b52 + b53 + b54), k5)
        
        temp_state = state + dt * (b61 * k1 + b62 * k2 + b63 * k3 + b64 * k4 + b65 * k5)
        call compute_derivatives(temp_state, time + dt * (b61 + b62 + b63 + b64 + b65), k6)
        
        ! 4阶解
        state_4th = state + dt * (c1 * k1 + c3 * k3 + c4 * k4 + c5 * k5 + c6 * k6)
        
        ! 5阶解
        state_5th = state + dt * (d1 * k1 + d3 * k3 + d4 * k4 + d5 * k5 + d6 * k6 + d7 * k6)
        
        ! 误差估计
        error_estimate = sqrt(sum((state_5th - state_4th)**2))
    end subroutine rkf45_step

    ! ======================================================================
    ! RKF 7(8) 单步推进核心
    ! 计算高精度的新状态，并返回截断误差估计向量
    ! ======================================================================
    subroutine rkf78_step(state, dt, time, state_7th, state_8th, error_estimate)
        real(DP), dimension(6), intent(in) :: state
        real(DP), intent(in) :: dt, time
        real(DP), dimension(6), intent(out) :: state_7th, state_8th
        real(DP), intent(out) :: error_estimate
        
        real(DP), dimension(6) :: f0, f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12
        
        ! --- RKF78 系数 (Fehlberg 7(8)) ---
        real(DP), parameter :: a1 = 2.0_DP/27.0_DP,   a2 = 1.0_DP/9.0_DP,   a3 = 1.0_DP/6.0_DP
        real(DP), parameter :: a4 = 5.0_DP/12.0_DP,   a5 = 1.0_DP/2.0_DP,   a6 = 5.0_DP/6.0_DP
        real(DP), parameter :: a7 = 1.0_DP/6.0_DP,    a8 = 2.0_DP/3.0_DP,   a9 = 1.0_DP/3.0_DP
        
        real(DP), parameter :: b10 = 2.0_DP/27.0_DP
        real(DP), parameter :: b20 = 1.0_DP/36.0_DP,    b21 = 1.0_DP/12.0_DP
        real(DP), parameter :: b30 = 1.0_DP/24.0_DP,    b32 = 1.0_DP/8.0_DP
        real(DP), parameter :: b40 = 5.0_DP/12.0_DP,    b42 = -25.0_DP/16.0_DP, b43 = 25.0_DP/16.0_DP
        real(DP), parameter :: b50 = 1.0_DP/20.0_DP,    b53 = 1.0_DP/4.0_DP,    b54 = 1.0_DP/5.0_DP
        real(DP), parameter :: b60 = -25.0_DP/108.0_DP, b63 = 125.0_DP/108.0_DP, &
                               b64 = -65.0_DP/27.0_DP,  b65 = 125.0_DP/54.0_DP
        real(DP), parameter :: b70 = 31.0_DP/300.0_DP,  b74 = 61.0_DP/225.0_DP,  &
                               b75 = -2.0_DP/9.0_DP,    b76 = 13.0_DP/900.0_DP
                               
        ! 从这里开始进行安全折行处理
        real(DP), parameter :: b80 = 2.0_DP,            b83 = -53.0_DP/6.0_DP,   &
                               b84 = 704.0_DP/45.0_DP,  b85 = -107.0_DP/9.0_DP,  &
                               b86 = 67.0_DP/90.0_DP,   b87 = 3.0_DP
                               
        real(DP), parameter :: b90 = -91.0_DP/108.0_DP, b93 = 23.0_DP/108.0_DP,  &
                               b94 = -976.0_DP/135.0_DP, b95 = 311.0_DP/54.0_DP, &
                               b96 = -19.0_DP/60.0_DP,  b97 = 17.0_DP/6.0_DP,    &
                               b98 = -1.0_DP/12.0_DP
                               
        real(DP), parameter :: b100 = 2383.0_DP/4100.0_DP, b103 = -341.0_DP/164.0_DP, &
                               b104 = 4496.0_DP/1025.0_DP, b105 = -301.0_DP/82.0_DP,  &
                               b106 = 2133.0_DP/4100.0_DP, b107 = 45.0_DP/82.0_DP,    &
                               b108 = 45.0_DP/164.0_DP,    b109 = 18.0_DP/41.0_DP
                               
        real(DP), parameter :: b110 = 3.0_DP/205.0_DP,   b115 = -6.0_DP/41.0_DP,    &
                               b116 = -3.0_DP/205.0_DP,  b117 = -3.0_DP/41.0_DP,    &
                               b118 = 3.0_DP/41.0_DP,    b119 = 6.0_DP/41.0_DP
                               
        real(DP), parameter :: b120 = -1777.0_DP/4100.0_DP, b123 = -341.0_DP/164.0_DP, &
                               b124 = 4496.0_DP/1025.0_DP,  b125 = -289.0_DP/82.0_DP,  &
                               b126 = 2193.0_DP/4100.0_DP,  b127 = 51.0_DP/82.0_DP,    &
                               b128 = 33.0_DP/164.0_DP,     b129 = 12.0_DP/41.0_DP
        
        real(DP), parameter :: c5 = 34.0_DP/105.0_DP,    c6 = 9.0_DP/35.0_DP,       &
                               c7 = 9.0_DP/35.0_DP,      c8 = 9.0_DP/280.0_DP,      &
                               c9 = 9.0_DP/280.0_DP,     c11 = 41.0_DP/840.0_DP,    &
                               c12 = 41.0_DP/840.0_DP
                               
        real(DP), parameter :: err_factor = 41.0_DP/840.0_DP
        ! ==========================================
        ! 修正点 1：清理了 dt == 0.0_DP 时的历史变量名
        ! ==========================================
        if (dt == 0.0_DP) then
            state_7th = state
            state_8th = state
            error_estimate = 0.0_DP
            return
        end if

        ! 13 阶段导数计算
        call compute_derivatives(state, time, f0)
        
        call compute_derivatives(state + dt*(f0*b10), time + dt*a1, f1)
        
        call compute_derivatives(state + dt*(f0*b20 + f1*b21), time + dt*a2, f2)
        
        call compute_derivatives(state + dt*(f0*b30 + f2*b32), time + dt*a3, f3)
        
        call compute_derivatives(state + dt*(f0*b40 + f2*b42 + f3*b43), time + dt*a4, f4)
        
        call compute_derivatives(state + dt*(f0*b50 + f3*b53 + f4*b54), time + dt*a5, f5)
        
        call compute_derivatives(state + dt*(f0*b60 + f3*b63 + f4*b64 + f5*b65), &
                                 time + dt*a6, f6)
        
        call compute_derivatives(state + dt*(f0*b70 + f4*b74 + f5*b75 + f6*b76), &
                                 time + dt*a7, f7)
        
        call compute_derivatives(state + dt*(f0*b80 + f3*b83 + f4*b84 + f5*b85 + &
                                 f6*b86 + f7*b87), time + dt*a8, f8)
        
        call compute_derivatives(state + dt*(f0*b90 + f3*b93 + f4*b94 + f5*b95 + &
                                 f6*b96 + f7*b97 + f8*b98), time + dt*a9, f9)
        
        call compute_derivatives(state + dt*(f0*b100 + f3*b103 + f4*b104 + f5*b105 + &
                                 f6*b106 + f7*b107 + f8*b108 + f9*b109), &
                                 time + dt, f10)
        
        call compute_derivatives(state + dt*(f0*b110 + f5*b115 + f6*b116 + f7*b117 + &
                                 f8*b118 + f9*b119), time, f11)
        
        call compute_derivatives(state + dt*(f0*b120 + f3*b123 + f4*b124 + f5*b125 + &
                                 f6*b126 + f7*b127 + f8*b128 + f9*b129 + f11), &
                                 time + dt, f12)

        ! 8阶解 (作为主结果)
        state_8th = state + dt*(f5*c5 + f6*c6 + f7*c7 + f8*c8 + f9*c9 + &
                                f11*c11 + f12*c12)
        
        ! 7阶解 (用于误差评估)
        state_7th = state_8th - (err_factor * (f0 + f10 - f11 - f12) * dt)
        
        ! ==========================================
        ! 修正点 2：修复标量化计算时丢失变量声明的问题
        ! ==========================================
        error_estimate = sqrt(sum((state_8th - state_7th)**2))
        
    end subroutine rkf78_step

end module pod_integrator_module
