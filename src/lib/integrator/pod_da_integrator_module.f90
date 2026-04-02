module pod_da_integrator_module
    use pod_global, only: DP
    use pod_config, only: config
    use pod_da_force_model_module, only: da_compute_acceleration,current_epoch0
    use pod_dace_classes

    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    
    implicit none

    ! --- 积分器类型常量定义 (暴露给外部用户使用) ---
    integer, parameter, public :: METHOD_RKF45 = 1
    integer, parameter, public :: METHOD_RKF78 = 2
    
contains

    ! =========================================================
    ! 接口极致纯粹，积分器完全不知道 epoch0 的存在 (DA版)
    ! =========================================================
    subroutine da_compute_derivatives(state_unit, time_unit, derivatives_unit)
        type(AlgebraicVector), intent(in) :: state_unit
        real(DP), intent(in) :: time_unit   ! 积分器传来的无量纲相对时间
        type(AlgebraicVector), intent(inout) :: derivatives_unit
        
        real(DP) :: real_time
        type(AlgebraicVector) :: position, velocity, acceleration
        integer :: i
        
        call position%init(3)
        call velocity%init(3)
        call acceleration%init(3)
        if (derivatives_unit%size /= 6) call derivatives_unit%init(6)
        
        ! 1. 还原物理状态 (空间的量纲还原)
        do i = 1, 3
            position%elements(i) = state_unit%elements(i) * config%LU
            velocity%elements(i) = state_unit%elements(i+3) * config%VU
        end do
        
        ! 2. 调用物理引擎 (还原时间并计算加速度)
        real_time = current_epoch0 + time_unit * config%TU
        call da_compute_acceleration(position, velocity, real_time, acceleration)
        
        ! 3. 导数去量纲化
        do i = 1, 3
            derivatives_unit%elements(i) = velocity%elements(i) / config%VU
            derivatives_unit%elements(i+3) = acceleration%elements(i) / config%AccU
        end do
        
    end subroutine da_compute_derivatives
    
    ! =========================================================
    ! 经典 RK4 积分器 (DA版)
    ! =========================================================
    subroutine da_rk4_integrate(state, dt, time, new_state)
        type(AlgebraicVector), intent(in) :: state
        real(DP), intent(in) :: dt, time
        type(AlgebraicVector), intent(inout) :: new_state
        
        type(AlgebraicVector) :: k1, k2, k3, k4, temp_state
        real(DP) :: half_dt
        
        ! 局部 DA 向量初始化
        call k1%init(6); call k2%init(6); call k3%init(6); call k4%init(6)
        call temp_state%init(6)
        if (new_state%size /= 6) call new_state%init(6)
        
        half_dt = 0.5_DP * dt
        
        ! k1 = f(t, y)
        call da_compute_derivatives(state, time, k1)
        
        ! k2 = f(t + dt/2, y + dt*k1/2)
        temp_state = state + half_dt * k1
        call da_compute_derivatives(temp_state, time + half_dt, k2)
        
        ! k3 = f(t + dt/2, y + dt*k2/2)
        temp_state = state + half_dt * k2
        call da_compute_derivatives(temp_state, time + half_dt, k3)
        
        ! k4 = f(t + dt, y + dt*k3)
        temp_state = state + dt * k3
        call da_compute_derivatives(temp_state, time + dt, k4)
        
        ! y(t + dt) = y(t) + dt*(k1 + 2*k2 + 2*k3 + k4)/6
        new_state = state + (dt / 6.0_DP) * (k1 + 2.0_DP * k2 + 2.0_DP * k3 + k4)
        
    end subroutine da_rk4_integrate

      subroutine da_rkf45_integrate(state, dt, time, new_state)
        type(AlgebraicVector), intent(in) :: state
        real(DP), intent(in) :: dt, time
        type(AlgebraicVector), intent(out) :: new_state
        
        type(AlgebraicVector) :: temp_4th, temp_5th
        type(AlgebraicVector) :: err_est
        
        ! 1. 直接调用核心计算引擎
        call da_rkf45_step(state, dt, time, temp_4th, temp_5th, err_est)
        
        ! 2. 丢弃4阶解和误差，只把最高精度的5阶解作为结果返回
        new_state = temp_5th
    end subroutine da_rkf45_integrate

    subroutine da_rkf78_integrate(state, dt, time, new_state)
        type(AlgebraicVector), intent(in) :: state
        real(DP), intent(in) :: dt, time
        type(AlgebraicVector), intent(out) :: new_state
        
        type(AlgebraicVector) :: state_7th, state_8th
        type(AlgebraicVector) :: error_estimate
        real(DP) :: tolerance
        
        tolerance = config%propagation_abs_tol
        
        call da_rkf78_step(state, dt, time, state_7th, state_8th, error_estimate)
        
        ! 返回高阶解
        new_state = state_8th
        
        ! 简单的误差检查逻辑
        if (any(error_estimate%cons() > tolerance)) then
            ! 这里可以添加警告，但在自适应循环中，这个逻辑会被上层接管
        end if
    end subroutine da_rkf78_integrate
    ! =========================================================
    ! 自适应步长控制主循环 (DA版 WRMS)
    ! =========================================================
    subroutine da_adaptive_step_integrate(state, t_start, t_end, integrator_method, &
                                          times, states, n_steps, &
                                          max_steps_in, rel_tol_in, abs_tol_in, dt_min_in, dt_max_in)
        ! 核心物理参数
        type(AlgebraicVector), intent(in) :: state
        real(DP), intent(in) :: t_start, t_end
        integer, intent(in) :: integrator_method
        
        ! 输出参数
        real(DP), allocatable, dimension(:), intent(out) :: times
        type(AlgebraicVector), allocatable, dimension(:), intent(out) :: states
        integer, intent(out) :: n_steps
        
        ! 可选的控制参数
        integer, intent(in), optional :: max_steps_in
        real(DP), intent(in), optional :: rel_tol_in, abs_tol_in, dt_min_in, dt_max_in
        
        ! 内部工作变量
        integer :: max_steps
        real(DP) :: rel_tol, abs_tol, dt_min, dt_max
        real(DP) :: current_time, dt
        type(AlgebraicVector) :: current_state
        type(AlgebraicVector) :: next_state_4th, next_state_5th, next_state_7th, next_state_8th
        type(AlgebraicVector) :: next_state_high, error_estimate_vector
        
        ! 用于 WRMS 误差计算的 DP 数组
        real(DP), dimension(6) :: err_est_cons, current_state_cons, scale_vector
        real(DP) :: wrms_error, safety_factor, exp_power
        integer :: i
        
        ! 1. 智能参数绑定
        max_steps = config%max_propagation_steps
        if (present(max_steps_in)) max_steps = max_steps_in
        
        if (integrator_method == METHOD_RKF45) then
            abs_tol = config%rkf45_abs_tol
            rel_tol = config%rkf45_rel_tol
            dt_min = config%rkf45_min_step / config%TU
            dt_max = config%rkf45_max_step / config%TU
        else if (integrator_method == METHOD_RKF78) then
            abs_tol = config%rkf78_abs_tol
            rel_tol = config%rkf78_rel_tol
            dt_min = config%rkf78_min_step / config%TU
            dt_max = config%rkf78_max_step / config%TU
        else
            print *, "Error: Unsupported integrator method!"
            stop
        end if
        
        if (present(rel_tol_in)) rel_tol = rel_tol_in
        if (present(abs_tol_in)) abs_tol = abs_tol_in
        if (present(dt_min_in))  dt_min  = dt_min_in/config%TU
        if (present(dt_max_in))  dt_max  = dt_max_in/config%TU
        
        safety_factor = 0.9_DP 
        
        ! 2. 内存分配
        if (allocated(times)) deallocate(times)
        if (allocated(states)) deallocate(states)
        allocate(times(max_steps))
        allocate(states(max_steps))
        
        call current_state%init(6)
        call next_state_4th%init(6); call next_state_5th%init(6)
        call next_state_7th%init(6); call next_state_8th%init(6)
        call error_estimate_vector%init(6)
        
        current_state = state
        current_time = t_start
        times(1) = current_time
        call states(1)%init(6)
        states(1) = current_state
        
        dt = min(dt_max, (t_end - t_start) / 100.0_DP)
        n_steps = 1
        
        do i = 2, max_steps
            if (current_time >= t_end) exit
            
            ! 3. 积分器路由分发
            if (integrator_method == METHOD_RKF78) then
                call da_rkf78_step(current_state, dt, current_time, next_state_7th, next_state_8th, error_estimate_vector)
                next_state_high = next_state_8th
                exp_power = 0.125_DP              ! 1/8
            else if (integrator_method == METHOD_RKF45) then
                call da_rkf45_step(current_state, dt, current_time, next_state_4th, next_state_5th, error_estimate_vector)
                next_state_high = next_state_5th
                exp_power = 0.20_DP               ! 1/5 (DOPRI5)
            end if
            
            ! 4. WRMS 误差评估 (提取常数部分进行标称轨迹误差控制)
            current_state_cons = current_state%cons()
            err_est_cons = error_estimate_vector%cons()
            
            scale_vector = abs_tol + rel_tol * abs(current_state_cons)
            wrms_error = sqrt(sum((err_est_cons / scale_vector)**2) / 6.0_DP)
            ! wrms_error = norm2(err_est_cons)/abs_tol   

            if (ieee_is_nan(wrms_error)) then 
                print *, "FATAL: WRMS 误差产生 NaN，强行退出积分。"
                exit
            end if

            if (wrms_error <= 1.0_DP) then
                ! 接受这一步
                current_state = next_state_high
                current_time = current_time + dt
                n_steps = n_steps + 1
                
                times(n_steps) = current_time
                call states(n_steps)%init(6)
                states(n_steps) = current_state
                
                dt = safety_factor * dt * (1.0_DP / max(wrms_error, 1.0e-15_DP))**exp_power
                dt = max(dt_min, min(dt_max, dt))
            else
                ! 拒绝这一步
                if (dt <= dt_min) then
                    ! 防抱死
                    current_state = next_state_high
                    current_time = current_time + dt
                    n_steps = n_steps + 1
                    
                    times(n_steps) = current_time
                    call states(n_steps)%init(6)
                    states(n_steps) = current_state
                else
                    dt = safety_factor * dt * (1.0_DP / max(wrms_error, 1.0e-15_DP))**exp_power
                    dt = max(dt_min, dt)
                end if
            end if
            
            if (current_time + dt > t_end) then
                dt = t_end - current_time
            end if
        end do
        
        ! 调整数组大小
        if (n_steps < max_steps) then
            block
                real(DP), allocatable, dimension(:) :: temp_times
                type(AlgebraicVector), allocatable, dimension(:) :: temp_states
                
                allocate(temp_times(n_steps))
                allocate(temp_states(n_steps))
                
                temp_times = times(1:n_steps)
                ! DA数组需要手动转移/重新赋值以防止底层句柄丢失
                do i = 1, n_steps
                    call temp_states(i)%init(6)
                    temp_states(i) = states(i)
                end do
                
                call move_alloc(temp_times, times)
                call move_alloc(temp_states, states)
            end block
        end if
    end subroutine da_adaptive_step_integrate

    ! =========================================================
    ! RKF45 单步积分器 (DOPRI5 修复版)
    ! =========================================================
    subroutine da_rkf45_step(state, dt, time, state_4th, state_5th, error_estimate_vector)
        type(AlgebraicVector), intent(in) :: state
        real(DP), intent(in) :: dt, time
        type(AlgebraicVector), intent(inout) :: state_4th, state_5th
        type(AlgebraicVector), intent(inout) :: error_estimate_vector
        
        type(AlgebraicVector) :: k1, k2, k3, k4, k5, k6, k7
        type(AlgebraicVector) :: temp_state
        
        ! DOPRI5 系数
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
        
        call k1%init(6); call k2%init(6); call k3%init(6)
        call k4%init(6); call k5%init(6); call k6%init(6); call k7%init(6)
        call temp_state%init(6)
        if (state_4th%size /= 6) call state_4th%init(6)
        if (state_5th%size /= 6) call state_5th%init(6)
        if (error_estimate_vector%size /= 6) call error_estimate_vector%init(6)
        
        call da_compute_derivatives(state, time, k1)
        
        temp_state = state + (dt * b21) * k1
        call da_compute_derivatives(temp_state, time + dt * b21, k2)
        
        temp_state = state + dt * (b31 * k1 + b32 * k2)
        call da_compute_derivatives(temp_state, time + dt * (b31 + b32), k3)
        
        temp_state = state + dt * (b41 * k1 + b42 * k2 + b43 * k3)
        call da_compute_derivatives(temp_state, time + dt * (b41 + b42 + b43), k4)
        
        temp_state = state + dt * (b51 * k1 + b52 * k2 + b53 * k3 + b54 * k4)
        call da_compute_derivatives(temp_state, time + dt * (b51 + b52 + b53 + b54), k5)
        
        temp_state = state + dt * (b61 * k1 + b62 * k2 + b63 * k3 + b64 * k4 + b65 * k5)
        call da_compute_derivatives(temp_state, time + dt * (b61 + b62 + b63 + b64 + b65), k6)

        ! 1. 5阶主解
        state_5th = state + dt * (c1 * k1 + c3 * k3 + c4 * k4 + c5 * k5 + c6 * k6)
        
        ! 2. 补齐 k7 
        call da_compute_derivatives(state_5th, time + dt, k7)

        ! 3. 4阶解
        state_4th = state + dt * (d1 * k1 + d3 * k3 + d4 * k4 + d5 * k5 + d6 * k6 + d7 * k7)

        ! 4. 输出代数差异向量，交给上层提取 cons()
        error_estimate_vector = state_5th - state_4th
        
    end subroutine da_rkf45_step

    ! ======================================================================
    ! RKF 7(8) 单步推进核心
    ! ======================================================================
    subroutine da_rkf78_step(state, dt, time, state_7th, state_8th, error_estimate_vector)
        type(AlgebraicVector), intent(in) :: state
        real(DP), intent(in) :: dt, time
        type(AlgebraicVector), intent(inout) :: state_7th, state_8th
        type(AlgebraicVector), intent(inout) :: error_estimate_vector
        
        type(AlgebraicVector) :: f0, f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12
        
        ! --- RKF78 系数 ---
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

        call f0%init(6); call f1%init(6); call f2%init(6); call f3%init(6); call f4%init(6)
        call f5%init(6); call f6%init(6); call f7%init(6); call f8%init(6); call f9%init(6)
        call f10%init(6); call f11%init(6); call f12%init(6)
        
        if (state_7th%size /= 6) call state_7th%init(6)
        if (state_8th%size /= 6) call state_8th%init(6)
        if (error_estimate_vector%size /= 6) call error_estimate_vector%init(6)

        if (dt <= 1.0e-15_DP) then
            state_7th = state
            state_8th = state
            error_estimate_vector = 0.0_DP
            return
        end if

        call da_compute_derivatives(state, time, f0)
        call da_compute_derivatives(state + dt*(f0*b10), time + dt*a1, f1)
        call da_compute_derivatives(state + dt*(f0*b20 + f1*b21), time + dt*a2, f2)
        call da_compute_derivatives(state + dt*(f0*b30 + f2*b32), time + dt*a3, f3)
        call da_compute_derivatives(state + dt*(f0*b40 + f2*b42 + f3*b43), time + dt*a4, f4)
        call da_compute_derivatives(state + dt*(f0*b50 + f3*b53 + f4*b54), time + dt*a5, f5)
        call da_compute_derivatives(state + dt*(f0*b60 + f3*b63 + f4*b64 + f5*b65), time + dt*a6, f6)
        call da_compute_derivatives(state + dt*(f0*b70 + f4*b74 + f5*b75 + f6*b76), time + dt*a7, f7)
        call da_compute_derivatives(state + dt*(f0*b80 + f3*b83 + f4*b84 + f5*b85 + f6*b86 + f7*b87), time + dt*a8, f8)
        call da_compute_derivatives(state + dt*(f0*b90 + f3*b93 + f4*b94 + f5*b95 + f6*b96 + f7*b97 + f8*b98), time + dt*a9, f9)
        call da_compute_derivatives(state + dt*(f0*b100 + f3*b103 + f4*b104 + f5*b105 + f6*b106 + f7*b107 + f8*b108 + f9*b109),&
                                     time + dt, f10)
        call da_compute_derivatives(state + dt*(f0*b110 + f5*b115 + f6*b116 + f7*b117 + f8*b118 + f9*b119),&
                                     time, f11)
        call da_compute_derivatives(state + dt*(f0*b120 + f3*b123 + f4*b124 + f5*b125 + f6*b126 + f7*b127 + f8*b128 + f9*b129 + &
                                    f11), time + dt, f12)

        ! 8阶主解
        state_8th = state + dt*(f5*c5 + f6*c6 + f7*c7 + f8*c8 + f9*c9 + f11*c11 + f12*c12)
        
        ! 7阶解
        state_7th = state_8th - (err_factor * (f0 + f10 - f11 - f12) * dt)
        
        ! 输出代数差异向量
        error_estimate_vector = state_8th - state_7th
        
    end subroutine da_rkf78_step

end module pod_da_integrator_module