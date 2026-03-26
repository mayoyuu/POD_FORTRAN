module pod_da_integrator_module
    use pod_global, only: DP
    use pod_config, only: config
    ! 引入你的 DA 版力学模型 (请确保你有这个函数)
    use pod_da_force_model_module, only: da_compute_acceleration
    ! 引入我们的 DA 核心类
    use pod_dace_classes
    
    implicit none
    ! --- 积分器类型常量定义 (暴露给外部用户使用) ---
    integer, parameter, public :: METHOD_RKF45 = 1
    integer, parameter, public :: METHOD_RKF78 = 2
    
contains
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
        
        ! ! 销毁局部 DA 向量
        ! call k1%destroy(); call k2%destroy(); call k3%destroy(); call k4%destroy()
        ! call temp_state%destroy()
    end subroutine da_rk4_integrate

      subroutine da_rkf45_integrate(state, dt, time, new_state)
        type(AlgebraicVector), intent(in) :: state
        real(DP), intent(in) :: dt, time
        type(AlgebraicVector), intent(inout) :: new_state
        
        type(AlgebraicVector) :: temp_4th, temp_5th
        real(DP) :: err_est
        
        ! 1. 直接调用核心计算引擎
        call da_rkf45_step(state, dt, time, temp_4th, temp_5th, err_est)
        
        ! 2. 丢弃4阶解和误差，只把最高精度的5阶解作为结果返回
        new_state = temp_5th
    end subroutine da_rkf45_integrate

    subroutine da_rkf78_integrate(state, dt, time, new_state)
        type(AlgebraicVector), intent(in) :: state
        real(DP), intent(in) :: dt, time
        type(AlgebraicVector), intent(inout) :: new_state
        
        type(AlgebraicVector) :: state_7th, state_8th
        real(DP) :: error_estimate, tolerance
        
        tolerance = config%propagation_tolerance
        
        call da_rkf78_step(state, dt, time, state_7th, state_8th, error_estimate)
        
        ! 返回高阶解
        new_state = state_8th
        
        ! 简单的误差检查逻辑
        if (error_estimate > tolerance) then
            ! 这里可以添加警告，但在自适应循环中，这个逻辑会被上层接管
        end if
    end subroutine da_rkf78_integrate
    


    ! =========================================================
    ! 自适应步长控制主循环 (DA版)
    ! =========================================================
    subroutine da_adaptive_step_integrate(state, t_start, t_end, max_steps, tolerance, &
                                          integrator_method, times, states, n_steps)
        type(AlgebraicVector), intent(in) :: state
        real(DP), intent(in) :: t_start, t_end, tolerance
        integer, intent(in) :: max_steps
        integer, intent(in) :: integrator_method  ! 新增参数：选择使用 RKF45 还是 RKF78
        real(DP), allocatable, dimension(:), intent(out) :: times
        type(AlgebraicVector), allocatable, dimension(:), intent(out) :: states
        integer, intent(out) :: n_steps
        
        real(DP) :: current_time, dt, dt_min, dt_max
        type(AlgebraicVector) :: next_state_4th, next_state_5th, next_state_7th, next_state_8th
        type(AlgebraicVector) :: current_state, next_state_high
        real(DP) :: error_estimate, safety_factor, exp_power
        integer :: i
        
        dt_min = 1.0_DP
        dt_max = 3600.0_DP
        safety_factor = 0.9_DP
        
        allocate(times(max_steps))
        allocate(states(max_steps))
        
        ! do j = 1, max_steps
        !     call states(j)%init(6)
        ! end do
        
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
            
            ! call da_rkf45_step(current_state, dt, current_time, next_state_4th, next_state_5th, error_estimate)
             ! ==========================================
            ! 1. 积分器路由分发 (使用易读的常量)
            ! ==========================================
            if (integrator_method == METHOD_RKF78) then
                call da_rkf78_step(current_state, dt, current_time, next_state_7th, next_state_8th, error_estimate)
                next_state_high = next_state_8th
                exp_power = 0.125_DP              ! 1/8
            else if (integrator_method == METHOD_RKF45) then
                call da_rkf45_step(current_state, dt, current_time, next_state_4th, next_state_5th, error_estimate)
                next_state_high = next_state_5th
                exp_power = 0.25_DP               ! 1/4
            else
                ! 防御性编程：用户传了不支持的常量
                print *, "Error: Unsupported integrator method!"
                stop
            end if
            
            if (error_estimate <= tolerance) then
                current_state = next_state_high
                current_time = current_time + dt
                n_steps = n_steps + 1
                
                times(n_steps) = current_time
                ! 【在这里按需分配！】为当前步初始化 DA 容器
                call states(n_steps)%init(6)
                states(n_steps) = current_state
                
                
                if (error_estimate > 0.0_DP) then
                    dt = safety_factor * dt * (tolerance / error_estimate)**exp_power
                end if
                dt = max(dt_min, min(dt_max, dt))
            else
                ! 检查是否已经触底
                if (dt <= dt_min) then
                    ! 【防抱死补丁】已经缩小到极小步长依然无法满足容差，强行接受步长并前进！
                    ! 否则会陷入不推进时间的死循环
                    current_state = next_state_high
                    current_time = current_time + dt
                    n_steps = n_steps + 1
                    
                    times(n_steps) = current_time
                    
                    call states(n_steps)%init(6)
                    states(n_steps) = current_state
                    
                    ! 保持 dt 为 dt_min 继续尝试下一步
                else
                    ! 正常拒绝这一步，减小步长，重新计算当前时刻
                    dt = safety_factor * dt * (tolerance / error_estimate)**exp_power
                    dt = max(dt_min, dt)
                end if
            end if
            !     dt = safety_factor * dt * (tolerance / error_estimate)**0.25_DP
            !     dt = max(dt_min, dt)
            ! end if
            
            if (current_time + dt > t_end) then
                dt = t_end - current_time
            end if
        end do
        
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
    end subroutine da_rkf45_step


    subroutine da_rkf78_step(state, dt, time, state_7th, state_8th, error_estimate)
        type(AlgebraicVector), intent(in) :: state
        real(DP), intent(in) :: dt, time
        type(AlgebraicVector), intent(inout) :: state_7th, state_8th
        type(AlgebraicVector) :: state_diff
        real(DP), intent(out) :: error_estimate
        
        type(AlgebraicVector) :: f0, f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12
        
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

        call f0%init(6); call f1%init(6); call f2%init(6); call f3%init(6); call f4%init(6)
        call f5%init(6); call f6%init(6); call f7%init(6); call f8%init(6); call f9%init(6)
        call f10%init(6); call f11%init(6); call f12%init(6)
        if (state_7th%size /= 6) call state_7th%init(6)
        if (state_8th%size /= 6) call state_8th%init(6)

        if (dt == 0.0_DP) then
            state_7th = state
            state_8th = state
            error_estimate = 0.0_DP
            return
        end if

        ! 13 阶段导数计算
        call da_compute_derivatives(state, time, f0)
        
        call da_compute_derivatives(state + dt*(f0*b10), time + dt*a1, f1)
        
        call da_compute_derivatives(state + dt*(f0*b20 + f1*b21), time + dt*a2, f2)
        
        call da_compute_derivatives(state + dt*(f0*b30 + f2*b32), time + dt*a3, f3)
        
        call da_compute_derivatives(state + dt*(f0*b40 + f2*b42 + f3*b43), time + dt*a4, f4)
        
        call da_compute_derivatives(state + dt*(f0*b50 + f3*b53 + f4*b54), time + dt*a5, f5)
        
        call da_compute_derivatives(state + dt*(f0*b60 + f3*b63 + f4*b64 + f5*b65), &
                                 time + dt*a6, f6)
        
        call da_compute_derivatives(state + dt*(f0*b70 + f4*b74 + f5*b75 + f6*b76), &
                                 time + dt*a7, f7)
        
        call da_compute_derivatives(state + dt*(f0*b80 + f3*b83 + f4*b84 + f5*b85 + &
                                 f6*b86 + f7*b87), time + dt*a8, f8)
        
        call da_compute_derivatives(state + dt*(f0*b90 + f3*b93 + f4*b94 + f5*b95 + &
                                 f6*b96 + f7*b97 + f8*b98), time + dt*a9, f9)
        
        call da_compute_derivatives(state + dt*(f0*b100 + f3*b103 + f4*b104 + f5*b105 + &
                                 f6*b106 + f7*b107 + f8*b108 + f9*b109), &
                                 time + dt, f10)
        
        call da_compute_derivatives(state + dt*(f0*b110 + f5*b115 + f6*b116 + f7*b117 + &
                                 f8*b118 + f9*b119), time, f11)
        
        call da_compute_derivatives(state + dt*(f0*b120 + f3*b123 + f4*b124 + f5*b125 + &
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
        state_diff = state_8th - state_7th
        error_estimate = norm2(state_diff%cons())
        
    end subroutine da_rkf78_step

end module pod_da_integrator_module