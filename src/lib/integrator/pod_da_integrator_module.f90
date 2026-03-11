module pod_da_integrator_module
    use pod_global, only: DP
    use pod_config, only: config
    ! 引入你的 DA 版力学模型 (请确保你有这个函数)
    use pod_da_force_model_module, only: da_compute_acceleration
    ! 引入我们的 DA 核心类
    use pod_dace_classes
    
    implicit none
    
contains

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
        
        ! 销毁局部 DA 向量
        call k1%destroy(); call k2%destroy(); call k3%destroy(); call k4%destroy()
        call temp_state%destroy()
    end subroutine da_rk4_integrate
    

    ! =========================================================
    ! 计算导数 (DA版)
    ! =========================================================
    subroutine da_compute_derivatives(state, time, derivatives)
        type(AlgebraicVector), intent(in) :: state
        real(DP), intent(in) :: time
        type(AlgebraicVector), intent(inout) :: derivatives
        
        type(AlgebraicVector) :: position, velocity, acceleration
        integer :: i
        
        call position%init(3)
        call velocity%init(3)
        call acceleration%init(3)
        
        ! 提取位置和速度
        do i = 1, 3
            position%elements(i) = state%elements(i)
            velocity%elements(i) = state%elements(i+3)
        end do
        
        ! 计算加速度 (调用专用的 DA 加速度计算函数)
        call da_compute_acceleration(position, velocity, time, acceleration)
        
        ! 构建导数向量 [dx/dt, dy/dt, dz/dt, dvx/dt, dvy/dt, dvz/dt]
        do i = 1, 3
            derivatives%elements(i) = velocity%elements(i)
            derivatives%elements(i+3) = acceleration%elements(i)
        end do
        
        call position%destroy()
        call velocity%destroy()
        call acceleration%destroy()
    end subroutine da_compute_derivatives
    

    ! =========================================================
    ! 自适应步长控制主循环 (DA版)
    ! =========================================================
    subroutine da_adaptive_step_integrate(state, t_start, t_end, max_steps, tolerance, &
                                          times, states, n_steps)
        type(AlgebraicVector), intent(in) :: state
        real(DP), intent(in) :: t_start, t_end, tolerance
        integer, intent(in) :: max_steps
        real(DP), allocatable, dimension(:), intent(out) :: times
        type(AlgebraicVector), allocatable, dimension(:), intent(out) :: states
        integer, intent(out) :: n_steps
        
        real(DP) :: current_time, dt, dt_min, dt_max
        type(AlgebraicVector) :: current_state, next_state_4th, next_state_5th
        real(DP) :: error_estimate, safety_factor
        integer :: i, j
        
        dt_min = 1.0_DP
        dt_max = 3600.0_DP
        safety_factor = 0.9_DP
        
        allocate(times(max_steps))
        allocate(states(max_steps))
        
        do j = 1, max_steps
            call states(j)%init(6)
        end do
        
        call current_state%init(6)
        call next_state_4th%init(6)
        call next_state_5th%init(6)
        
        current_state = state
        current_time = t_start
        times(1) = current_time
        states(1) = current_state
        
        dt = min(dt_max, (t_end - t_start) / 100.0_DP)
        n_steps = 1
        
        do i = 2, max_steps
            if (current_time >= t_end) exit
            
            call da_rkf45_step(current_state, dt, current_time, next_state_4th, next_state_5th, error_estimate)
            
            if (error_estimate <= tolerance) then
                current_state = next_state_5th
                current_time = current_time + dt
                n_steps = n_steps + 1
                
                times(n_steps) = current_time
                states(n_steps) = current_state
                
                if (error_estimate > 0.0_DP) then
                    dt = safety_factor * dt * (tolerance / error_estimate)**0.25_DP
                end if
                dt = max(dt_min, min(dt_max, dt))
            else
                dt = safety_factor * dt * (tolerance / error_estimate)**0.25_DP
                dt = max(dt_min, dt)
            end if
            
            if (current_time + dt > t_end) then
                dt = t_end - current_time
            end if
        end do
        
        call current_state%destroy()
        call next_state_4th%destroy()
        call next_state_5th%destroy()
    end subroutine da_adaptive_step_integrate
    

    ! =========================================================
    ! RKF45 单步积分器 (DA版)
    ! =========================================================
    subroutine da_rkf45_step(state, dt, time, state_4th, state_5th, error_estimate)
        type(AlgebraicVector), intent(in) :: state
        real(DP), intent(in) :: dt, time
        type(AlgebraicVector), intent(inout) :: state_4th, state_5th
        real(DP), intent(out) :: error_estimate
        
        type(AlgebraicVector) :: k1, k2, k3, k4, k5, k6, temp_state
        integer :: i
        
        ! RKF45 系数
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
        call k4%init(6); call k5%init(6); call k6%init(6)
        call temp_state%init(6)
        if (state_4th%size /= 6) call state_4th%init(6)
        if (state_5th%size /= 6) call state_5th%init(6)
        
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
        
        state_4th = state + dt * (c1 * k1 + c3 * k3 + c4 * k4 + c5 * k5 + c6 * k6)
        state_5th = state + dt * (d1 * k1 + d3 * k3 + d4 * k4 + d5 * k5 + d6 * k6 + d7 * k6)
        
        ! 利用常数项(标称轨迹)进行误差估计
        error_estimate = 0.0_DP
        do i = 1, 6
            error_estimate = error_estimate + &
                             (state_5th%elements(i)%cons() - state_4th%elements(i)%cons())**2
        end do
        error_estimate = sqrt(error_estimate)
        
        call k1%destroy(); call k2%destroy(); call k3%destroy()
        call k4%destroy(); call k5%destroy(); call k6%destroy()
        call temp_state%destroy()
    end subroutine da_rkf45_step

end module pod_da_integrator_module