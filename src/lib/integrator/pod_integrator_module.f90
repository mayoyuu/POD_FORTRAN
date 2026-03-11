module pod_integrator_module
    use pod_global, only: DP
    use pod_config, only: config
    use pod_force_model_module, only: compute_acceleration
    
    implicit none
    
contains

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
        
        real(DP), dimension(6) :: k1, k2, k3, k4, k5, k6
        real(DP), dimension(6) :: temp_state
        real(DP), dimension(6) :: state_4th, state_5th
        real(DP) :: error_estimate, tolerance
        
        tolerance = config%propagation_tolerance
        
        ! RKF45 积分方法 (Runge-Kutta-Fehlberg 4(5))
        ! k1 = f(t, y)
        call compute_derivatives(state, time, k1)
        
        ! k2 = f(t + dt/5, y + dt*k1/5)
        temp_state = state + dt * k1 / 5.0_DP
        call compute_derivatives(temp_state, time + dt/5.0_DP, k2)
        
        ! k3 = f(t + 3*dt/10, y + dt*(3*k1 + 9*k2)/40)
        temp_state = state + dt * (3.0_DP * k1 + 9.0_DP * k2) / 40.0_DP
        call compute_derivatives(temp_state, time + 3.0_DP*dt/10.0_DP, k3)
        
        ! k4 = f(t + 4*dt/5, y + dt*(44*k1 - 168*k2 + 160*k3)/45)
        temp_state = state + dt * (44.0_DP * k1 - 168.0_DP * k2 + 160.0_DP * k3) / 45.0_DP
        call compute_derivatives(temp_state, time + 4.0_DP*dt/5.0_DP, k4)
        
        ! k5 = f(t + 8*dt/9, y + dt*(19372*k1 - 76032*k2 + 64448*k3 - 212*k4)/6561)
        temp_state = state + dt * (19372.0_DP * k1 - 76032.0_DP * k2 + &
                                  64448.0_DP * k3 - 212.0_DP * k4) / 6561.0_DP
        call compute_derivatives(temp_state, time + 8.0_DP*dt/9.0_DP, k5)
        
        ! k6 = f(t + dt, y + dt*(9017*k1 - 355*k2 + 46732*k3 + 49*k4 - 5103*k5)/6561)
        temp_state = state + dt * (9017.0_DP * k1 - 355.0_DP * k2 + 46732.0_DP * k3 + &
                                  49.0_DP * k4 - 5103.0_DP * k5) / 6561.0_DP
        call compute_derivatives(temp_state, time + dt, k6)
        
        ! 4阶解
        state_4th = state + dt * (25.0_DP * k1 / 216.0_DP + 1408.0_DP * k3 / 2565.0_DP + &
                                  2197.0_DP * k4 / 4104.0_DP - k5 / 5.0_DP)
        
        ! 5阶解
        state_5th = state + dt * (16.0_DP * k1 / 135.0_DP + 6656.0_DP * k3 / 12825.0_DP + &
                                  28561.0_DP * k4 / 56430.0_DP - 9.0_DP * k5 / 50.0_DP + &
                                  2.0_DP * k6 / 55.0_DP)
        
        ! 误差估计
        error_estimate = sqrt(sum((state_5th - state_4th)**2))
        
        ! 使用5阶解作为最终结果
        new_state = state_5th
        
        ! 如果误差太大，可以在这里实现自适应步长控制
        if (error_estimate > tolerance) then
            ! 可以减小步长或发出警告
            ! 这里简化处理，实际应用中需要更复杂的自适应控制
        end if
    end subroutine rkf45_integrate
    
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
    
    subroutine adaptive_step_integrate(state, t_start, t_end, max_steps, tolerance, &
                                     times, states, n_steps)
        real(DP), dimension(6), intent(in) :: state
        real(DP), intent(in) :: t_start, t_end, tolerance
        integer, intent(in) :: max_steps
        real(DP), allocatable, dimension(:), intent(out) :: times
        real(DP), allocatable, dimension(:,:), intent(out) :: states
        integer, intent(out) :: n_steps
        
        real(DP) :: current_time, dt, dt_min, dt_max
        real(DP), dimension(6) :: current_state, next_state_4th, next_state_5th
        real(DP) :: error_estimate, safety_factor
        integer :: i
        
        ! 初始化
        dt_min = 1.0_DP  ! 最小步长 (秒)
        dt_max = 3600.0_DP  ! 最大步长 (秒)
        safety_factor = 0.9_DP
        
        ! 分配内存
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
            
            ! 使用RKF45方法进行一步积分
            call rkf45_step(current_state, dt, current_time, next_state_4th, next_state_5th, error_estimate)
            
            ! 检查误差
            if (error_estimate <= tolerance) then
                ! 接受这一步
                current_state = next_state_5th
                current_time = current_time + dt
                n_steps = n_steps + 1
                
                times(n_steps) = current_time
                states(n_steps, :) = current_state
                
                ! 调整步长
                if (error_estimate > 0.0_DP) then
                    dt = safety_factor * dt * (tolerance / error_estimate)**0.25_DP
                end if
                dt = max(dt_min, min(dt_max, dt))
                
            else
                ! 检查是否已经触底
                if (dt <= dt_min) then
                    ! 【防抱死补丁】已经缩小到极小步长依然无法满足容差，强行接受步长并前进！
                    ! 否则会陷入不推进时间的死循环
                    current_state = next_state_5th
                    current_time = current_time + dt
                    n_steps = n_steps + 1
                    
                    times(n_steps) = current_time
                    states(n_steps, :) = current_state
                    
                    ! 保持 dt 为 dt_min 继续尝试下一步
                else
                    ! 正常拒绝这一步，减小步长，重新计算当前时刻
                    dt = safety_factor * dt * (tolerance / error_estimate)**0.25_DP
                    dt = max(dt_min, dt)
                end if
            end if
            !     ! 拒绝这一步，减小步长
            !     dt = safety_factor * dt * (tolerance / error_estimate)**0.25_DP
            !     dt = max(dt_min, dt)
            ! end if
            
            ! 确保不超过结束时间
            if (current_time + dt > t_end) then
                dt = t_end - current_time
            end if
        end do
        
        ! 调整数组大小
        if (n_steps < max_steps) then
            times = times(1:n_steps)
            states = states(1:n_steps, :)
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

end module pod_integrator_module
