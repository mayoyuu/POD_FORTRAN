!--------------------------------------------------------------------------------------------------------------
!> STT 误差传播模块
!>
!> 实现基于 State Transition Tensor 的非线性不确定性传播。
!> 支持可配置的 STT 阶数 (2~6) 和联合积分。
!>
!> 架构:
!>   - 变维度 RKF78 积分器 (从 pod_integrator_module 泛化)
!>   - uq_stt_propagator 类型 (继承 uq_propagator_base)
!>   - 增广 ODE 系统: 标称轨道 + STT 联合积分
!>
!> 依赖:
!>   - pod_global, pod_uq_base_module, pod_uq_state_module
!>   - pod_stt_tensor_module, pod_crtbp_derivatives_module
!>   - pod_integrator_module (仅常量)
!--------------------------------------------------------------------------------------------------------------
module pod_uq_prop_stt_module
    use pod_global, only: DP, MAX_STRING_LEN
    use pod_uq_base_module, only: uq_propagator_base
    use pod_uq_state_module, only: uq_state_type
    use pod_integrator_module, only: METHOD_RKF45, METHOD_RKF78
    use pod_stt_tensor_module, only: STT_MAX_ORDER, STT_DIM, stt_sizes, &
                                      init_stt_indexing, stt_order_p_size, &
                                      tuple_to_sym_index, sym_index_to_tuple, &
                                      compute_stt_rhs, compute_stt_rhs_order, &
                                      compute_stt_moments, stt_store_type, &
                                      stt_flat_offset, generate_partitions, &
                                      gen_block_assignments, factorial
    use pod_crtbp_derivatives_module, only: crtbp_force_derivatives, crtbp_derivatives_init
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    implicit none
    private

    public :: uq_stt_propagator, stt_propagate_deviates

    ! =========================================================================
    ! STT 传播器类型
    ! =========================================================================
    type, extends(uq_propagator_base), public :: uq_stt_propagator
        integer :: stt_order = 2
        integer :: total_dim = 0
        real(DP) :: stt_abs_tol = 1.0d-12
        real(DP) :: stt_rel_tol = 1.0d-12
        real(DP) :: stt_dt_min = 1.0d-6
        real(DP) :: stt_dt_max = 3600.0_DP
        integer :: stt_max_steps = 100000
        real(DP) :: mu = 0.012153614091892_DP
        real(DP), allocatable :: aug_state(:)
    contains
        procedure :: propagate => stt_propagate
        procedure :: get_method_name => stt_get_method_name
        procedure :: set_stt_order
        procedure :: set_stt_tolerances
        procedure :: set_stt_mu
    end type uq_stt_propagator

    ! =========================================================================
    ! 增广 RHS 回调的抽象接口
    ! =========================================================================
    abstract interface
        subroutine aug_derivs_interface(y, t, dydt, mu, order, total_dim)
            import :: DP
            real(DP), intent(in) :: y(:), t
            real(DP), intent(out) :: dydt(:)
            real(DP), intent(in) :: mu
            integer, intent(in) :: order, total_dim
        end subroutine
    end interface

contains

    ! =====================================================================
    ! 方法名
    ! =====================================================================
    function stt_get_method_name(this) result(name)
        class(uq_stt_propagator), intent(in) :: this
        character(len=MAX_STRING_LEN) :: name
        write(name, '(A,I0)') 'STT-', this%stt_order
    end function stt_get_method_name

    ! =====================================================================
    ! 设置 STT 阶数
    ! =====================================================================
    subroutine set_stt_order(this, order)
        class(uq_stt_propagator), intent(inout) :: this
        integer, intent(in) :: order
        if (order >= 2 .and. order <= STT_MAX_ORDER) then
            this%stt_order = order
        else
            write(*,*) '[STT] ERROR: order must be 2..6, got:', order
        end if
    end subroutine set_stt_order

    ! =====================================================================
    ! 设置积分容差
    ! =====================================================================
    subroutine set_stt_tolerances(this, abs_tol, rel_tol, dt_min, dt_max, max_steps)
        class(uq_stt_propagator), intent(inout) :: this
        real(DP), optional, intent(in) :: abs_tol, rel_tol, dt_min, dt_max
        integer, optional, intent(in) :: max_steps
        if (present(abs_tol)) this%stt_abs_tol = abs_tol
        if (present(rel_tol)) this%stt_rel_tol = rel_tol
        if (present(dt_min)) this%stt_dt_min = dt_min
        if (present(dt_max)) this%stt_dt_max = dt_max
        if (present(max_steps)) this%stt_max_steps = max_steps
    end subroutine set_stt_tolerances

    ! =====================================================================
    ! 设置 CRTBP mu 参数
    ! =====================================================================
    subroutine set_stt_mu(this, mu)
        class(uq_stt_propagator), intent(inout) :: this
        real(DP), intent(in) :: mu
        this%mu = mu
    end subroutine set_stt_mu

    ! =====================================================================
    ! 计算 STT 增广 ODE 右端项
    !
    ! aug_state 布局:
    !   [ x(1:6)                            ] — 标称轨道
    !   [ Φ¹_{1,1..6}, ..., Φ¹_{6,1..6}    ] — STM 扁平化 (36)
    !   [ Φ²_{1,1..N₂}, ..., Φ²_{6,1..N₂}  ] — 2 阶 STT 扁平化 (6*21=126)
    !   [ ... 高阶同理 ...                   ]
    ! =====================================================================
    subroutine compute_stt_aug_derivatives(y, t, dydt, mu, order, total_dim)
        real(DP), intent(in) :: y(:), t
        real(DP), intent(out) :: dydt(:)
        real(DP), intent(in) :: mu
        integer, intent(in) :: order, total_dim

        real(DP) :: x_nominal(6)
        real(DP), allocatable :: f_tensors(:, :, :)
        real(DP) :: nominal_f(6)
        integer :: pos, p, np, i_comp, tidx, ia
        real(DP) :: mu_

        mu_ = mu

        ! ---- 1. 解包标称状态 ----
        x_nominal = y(1:6)

        ! ---- 2. 标称导数 ----
        call crtbp_force_real(x_nominal, mu_, nominal_f)
        dydt(1:6) = nominal_f

        ! ---- 3. 计算 f* 张量 ----
        call crtbp_force_derivatives(x_nominal, order, f_tensors)

        ! ---- 4. STM RHS (order 1): Φ̇_{i,a} = f*_{i,α} Φ_{α,a} ----
        ! 使用扁平数组 y(7:) 直接访问, 无需 stt_store 解包
        pos = 7
        np = stt_sizes(1)  ! = 6
        do i_comp = 1, 6
            do tidx = 1, np
                dydt(pos) = 0.0_DP
                do ia = 1, 6
                    dydt(pos) = dydt(pos) + f_tensors(i_comp, ia, 1) &
                                * y(7 + stt_flat_offset(ia, tidx, 1) - 1)
                end do
                pos = pos + 1
            end do
        end do

        ! ---- 5. 高阶 STT RHS (order 2..max_order) ----
        do p = 2, order
            call compute_stt_rhs_order(p, y(7:), f_tensors, order, dydt, pos)
            pos = pos + 6 * stt_sizes(p)
        end do

        if (allocated(f_tensors)) deallocate(f_tensors)
    end subroutine compute_stt_aug_derivatives

    ! =====================================================================
    ! CRTBP 标称力 (局部副本)
    ! =====================================================================
    subroutine crtbp_force_real(x, mu, f)
        real(DP), intent(in) :: x(6), mu
        real(DP), intent(out) :: f(6)
        real(DP) :: r1, r2, r13, r23, omx, omy, omz

        r1 = sqrt((x(1) + mu)**2 + x(2)**2 + x(3)**2)
        r2 = sqrt((x(1) - 1.0_DP + mu)**2 + x(2)**2 + x(3)**2)
        r13 = r1**3
        r23 = r2**3

        omx = x(1) - (1.0_DP - mu)*(x(1) + mu)/r13 - mu*(x(1) - 1.0_DP + mu)/r23
        omy = x(2) - (1.0_DP - mu)*x(2)/r13 - mu*x(2)/r23
        omz = -(1.0_DP - mu)*x(3)/r13 - mu*x(3)/r23

        f(1) = x(4)
        f(2) = x(5)
        f(3) = x(6)
        f(4) = omx + 2.0_DP*x(5)
        f(5) = omy - 2.0_DP*x(4)
        f(6) = omz
    end subroutine crtbp_force_real

    ! =====================================================================
    ! 主传播接口
    ! =====================================================================
    subroutine stt_propagate(this, t_start, t_end, input_state, output_state)
        class(uq_stt_propagator), intent(inout) :: this
        real(DP), intent(in) :: t_start, t_end
        type(uq_state_type), intent(in) :: input_state
        type(uq_state_type), intent(inout) :: output_state

        real(DP), allocatable :: times(:), states(:,:)
        integer :: n_steps, dim, p, np, pos, i_comp, tidx
        type(stt_store_type) :: stt_store
        real(DP) :: P0(6,6), mean_final(6), cov_final(6,6), x_star(6)

        character(len=64) :: msg

        ! ---- 1. 验证初始化 ----
        if (this%stt_order < 2 .or. this%stt_order > STT_MAX_ORDER) then
            write(*,*) '[STT] ERROR: call set_stt_order() before propagate()'
            return
        end if

        ! ---- 2. 确保 tensor 索引系统已初始化 ----
        call init_stt_indexing(this%stt_order)
        call crtbp_derivatives_init(this%mu)

        ! ---- 3. 计算增广系统维度 ----
        this%total_dim = 6
        do p = 1, this%stt_order
            this%total_dim = this%total_dim + 6 * stt_sizes(p)
        end do

        if (this%verbose) then
            write(msg, '(A,I0,A,I0)') '[STT] order=', this%stt_order, &
                ', augmented dimension=', this%total_dim
            write(*,*) trim(msg)
        end if

        ! ---- 4. 分配并初始化增广状态 ----
        allocate(this%aug_state(this%total_dim))
        this%aug_state = 0.0_DP

        ! 标称初始状态
        if (.not. allocated(input_state%mean)) then
            write(*,*) '[STT] ERROR: input_state%mean not allocated'
            return
        end if
        this%aug_state(1:6) = input_state%mean(1:6)

        ! STM 初始化为 I₆
        do i_comp = 1, 6
            this%aug_state(7 + (i_comp-1)*6 + (i_comp-1)) = 1.0_DP
        end do
        ! 高阶 STT 初始化为 0 (已由整体置零完成)

        ! ---- 5. 联合积分 ----
        call adaptive_integrate_var_dim(this%aug_state, this%total_dim, &
            t_start, t_end, &
            this%stt_abs_tol, this%stt_rel_tol, &
            this%stt_dt_min, this%stt_dt_max, this%stt_max_steps, &
            this%stt_order, this%mu, &
            times, states, n_steps)

        if (this%verbose) then
            write(*,'(A,I0,A)') '[STT] integration complete: ', n_steps, ' steps'
        end if

        ! ---- 6. 提取终点 STT 并计算矩 ----
        call stt_store%init(this%stt_order)

        ! 终点标称状态
        x_star = states(n_steps, 1:6)

        ! 解包 STT
        pos = 7
        do p = 1, this%stt_order
            np = stt_sizes(p)
            do i_comp = 1, 6
                do tidx = 1, np
                    call stt_store%set(i_comp, tidx, p, states(n_steps, pos))
                    pos = pos + 1
                end do
            end do
        end do

        ! 提取初始协方差
        P0 = 0.0_DP
        if (allocated(input_state%cov)) then
            P0(1:6, 1:6) = input_state%cov(1:6, 1:6)
        else
            ! 默认单位对角协方差
            do i_comp = 1, 6
                P0(i_comp, i_comp) = 1.0d-10
            end do
        end if

        call compute_stt_moments(x_star, stt_store, P0, this%stt_order, &
                                  mean_final, cov_final)

        ! ---- 7. 填充输出状态 ----
        call output_state%deallocate_memory()
        call output_state%allocate_memory(6, 0)
        if (allocated(output_state%mean)) deallocate(output_state%mean)
        if (allocated(output_state%cov)) deallocate(output_state%cov)
        allocate(output_state%mean(6), output_state%cov(6,6))
        output_state%mean = mean_final
        output_state%cov = cov_final

        ! ---- 8. 清理 ----
        call stt_store%destroy()
        deallocate(times, states, this%aug_state)
    end subroutine stt_propagate

    ! =====================================================================
    ! 变维度 RKF78 单步
    !
    ! 从 pod_integrator_module%rkf78_step 泛化而来,
    ! 将固定 dimension(6) 替换为 dimension(n)。
    ! RKF78 系数不变。
    ! =====================================================================
    subroutine rkf78_var_dim(state, n, dt, time, mu, order, total_dim, &
                              state_7th, state_8th, error_est)
        real(DP), intent(in) :: state(n), dt, time
        integer, intent(in) :: n
        real(DP), intent(in) :: mu
        integer, intent(in) :: order, total_dim
        real(DP), intent(out) :: state_7th(n), state_8th(n), error_est(n)

        real(DP), allocatable :: f0(:), f1(:), f2(:), f3(:), f4(:), f5(:), f6(:)
        real(DP), allocatable :: f7(:), f8(:), f9(:), f10(:), f11(:), f12(:)
        real(DP), allocatable :: temp(:)

        ! ---- RKF78 系数 (与经典实现相同) ----
        real(DP), parameter :: a1 = 2.0_DP/27.0_DP,   a2 = 1.0_DP/9.0_DP
        real(DP), parameter :: a3 = 1.0_DP/6.0_DP,    a4 = 5.0_DP/12.0_DP
        real(DP), parameter :: a5 = 1.0_DP/2.0_DP,    a6 = 5.0_DP/6.0_DP
        real(DP), parameter :: a7 = 1.0_DP/6.0_DP,    a8 = 2.0_DP/3.0_DP
        real(DP), parameter :: a9 = 1.0_DP/3.0_DP

        real(DP), parameter :: b10 = 2.0_DP/27.0_DP
        real(DP), parameter :: b20 = 1.0_DP/36.0_DP,    b21 = 1.0_DP/12.0_DP
        real(DP), parameter :: b30 = 1.0_DP/24.0_DP,    b32 = 1.0_DP/8.0_DP
        real(DP), parameter :: b40 = 5.0_DP/12.0_DP,    b42 = -25.0_DP/16.0_DP, &
                               b43 = 25.0_DP/16.0_DP
        real(DP), parameter :: b50 = 1.0_DP/20.0_DP,    b53 = 1.0_DP/4.0_DP, &
                               b54 = 1.0_DP/5.0_DP
        real(DP), parameter :: b60 = -25.0_DP/108.0_DP, b63 = 125.0_DP/108.0_DP, &
                               b64 = -65.0_DP/27.0_DP,  b65 = 125.0_DP/54.0_DP
        real(DP), parameter :: b70 = 31.0_DP/300.0_DP,  b74 = 61.0_DP/225.0_DP, &
                               b75 = -2.0_DP/9.0_DP,    b76 = 13.0_DP/900.0_DP
        real(DP), parameter :: b80 = 2.0_DP,            b83 = -53.0_DP/6.0_DP, &
                               b84 = 704.0_DP/45.0_DP,  b85 = -107.0_DP/9.0_DP, &
                               b86 = 67.0_DP/90.0_DP,   b87 = 3.0_DP
        real(DP), parameter :: b90 = -91.0_DP/108.0_DP, b93 = 23.0_DP/108.0_DP, &
                               b94 = -976.0_DP/135.0_DP, b95 = 311.0_DP/54.0_DP, &
                               b96 = -19.0_DP/60.0_DP,  b97 = 17.0_DP/6.0_DP, &
                               b98 = -1.0_DP/12.0_DP
        real(DP), parameter :: b100 = 2383.0_DP/4100.0_DP, b103 = -341.0_DP/164.0_DP, &
                               b104 = 4496.0_DP/1025.0_DP, b105 = -301.0_DP/82.0_DP, &
                               b106 = 2133.0_DP/4100.0_DP, b107 = 45.0_DP/82.0_DP, &
                               b108 = 45.0_DP/164.0_DP,    b109 = 18.0_DP/41.0_DP
        real(DP), parameter :: b110 = 3.0_DP/205.0_DP,   b115 = -6.0_DP/41.0_DP, &
                               b116 = -3.0_DP/205.0_DP,  b117 = -3.0_DP/41.0_DP, &
                               b118 = 3.0_DP/41.0_DP,    b119 = 6.0_DP/41.0_DP
        real(DP), parameter :: b120 = -1777.0_DP/4100.0_DP, b123 = -341.0_DP/164.0_DP, &
                               b124 = 4496.0_DP/1025.0_DP,  b125 = -289.0_DP/82.0_DP, &
                               b126 = 2193.0_DP/4100.0_DP,  b127 = 51.0_DP/82.0_DP, &
                               b128 = 33.0_DP/164.0_DP,     b129 = 12.0_DP/41.0_DP

        real(DP), parameter :: c5 = 34.0_DP/105.0_DP,  c6 = 9.0_DP/35.0_DP, &
                               c7 = 9.0_DP/35.0_DP,    c8 = 9.0_DP/280.0_DP, &
                               c9 = 9.0_DP/280.0_DP,   c11 = 41.0_DP/840.0_DP, &
                               c12 = 41.0_DP/840.0_DP
        real(DP), parameter :: err_factor = 41.0_DP/840.0_DP

        integer :: i

        if (dt <= 1.0e-15_DP) then
            state_7th = state
            state_8th = state
            error_est = 0.0_DP
            return
        end if

        allocate(f0(n), f1(n), f2(n), f3(n), f4(n), f5(n), f6(n))
        allocate(f7(n), f8(n), f9(n), f10(n), f11(n), f12(n))
        allocate(temp(n))

        ! 13 阶段计算
        call compute_stt_aug_derivatives(state, time, f0, mu, order, total_dim)

        temp = state + dt*(f0*b10)
        call compute_stt_aug_derivatives(temp, time + dt*a1, f1, mu, order, total_dim)

        temp = state + dt*(f0*b20 + f1*b21)
        call compute_stt_aug_derivatives(temp, time + dt*a2, f2, mu, order, total_dim)

        temp = state + dt*(f0*b30 + f2*b32)
        call compute_stt_aug_derivatives(temp, time + dt*a3, f3, mu, order, total_dim)

        temp = state + dt*(f0*b40 + f2*b42 + f3*b43)
        call compute_stt_aug_derivatives(temp, time + dt*a4, f4, mu, order, total_dim)

        temp = state + dt*(f0*b50 + f3*b53 + f4*b54)
        call compute_stt_aug_derivatives(temp, time + dt*a5, f5, mu, order, total_dim)

        temp = state + dt*(f0*b60 + f3*b63 + f4*b64 + f5*b65)
        call compute_stt_aug_derivatives(temp, time + dt*a6, f6, mu, order, total_dim)

        temp = state + dt*(f0*b70 + f4*b74 + f5*b75 + f6*b76)
        call compute_stt_aug_derivatives(temp, time + dt*a7, f7, mu, order, total_dim)

        temp = state + dt*(f0*b80 + f3*b83 + f4*b84 + f5*b85 + f6*b86 + f7*b87)
        call compute_stt_aug_derivatives(temp, time + dt*a8, f8, mu, order, total_dim)

        temp = state + dt*(f0*b90 + f3*b93 + f4*b94 + f5*b95 + f6*b96 + f7*b97 + f8*b98)
        call compute_stt_aug_derivatives(temp, time + dt*a9, f9, mu, order, total_dim)

        temp = state + dt*(f0*b100 + f3*b103 + f4*b104 + f5*b105 + &
                           f6*b106 + f7*b107 + f8*b108 + f9*b109)
        call compute_stt_aug_derivatives(temp, time + dt, f10, mu, order, total_dim)

        temp = state + dt*(f0*b110 + f5*b115 + f6*b116 + f7*b117 + f8*b118 + f9*b119)
        call compute_stt_aug_derivatives(temp, time, f11, mu, order, total_dim)

        temp = state + dt*(f0*b120 + f3*b123 + f4*b124 + f5*b125 + &
                           f6*b126 + f7*b127 + f8*b128 + f9*b129 + f11)
        call compute_stt_aug_derivatives(temp, time + dt, f12, mu, order, total_dim)

        ! 8 阶解
        state_8th = state + dt*(f5*c5 + f6*c6 + f7*c7 + f8*c8 + f9*c9 + &
                                f11*c11 + f12*c12)

        ! 7 阶解 (用于误差估计)
        state_7th = state_8th - (err_factor * (f0 + f10 - f11 - f12) * dt)

        error_est = abs(state_8th - state_7th)

        deallocate(f0, f1, f2, f3, f4, f5, f6)
        deallocate(f7, f8, f9, f10, f11, f12)
        deallocate(temp)
    end subroutine rkf78_var_dim

    ! =====================================================================
    ! 变维度自适应积分器
    ! =====================================================================
    subroutine adaptive_integrate_var_dim(state, n, t_start, t_end, &
                                           abs_tol, rel_tol, dt_min, dt_max, max_steps, &
                                           order, mu, &
                                           times, states, n_steps)
        real(DP), intent(inout) :: state(n)
        integer, intent(in) :: n
        real(DP), intent(in) :: t_start, t_end
        real(DP), intent(in) :: abs_tol, rel_tol, dt_min, dt_max
        integer, intent(in) :: max_steps
        integer, intent(in) :: order
        real(DP), intent(in) :: mu
        real(DP), allocatable, intent(out) :: times(:), states(:,:)
        integer, intent(out) :: n_steps

        real(DP) :: current_time, dt
        real(DP), allocatable :: current_state(:)
        real(DP), allocatable :: next_state_7th(:), next_state_8th(:)
        real(DP), allocatable :: error_est(:), scale_vector(:)
        real(DP) :: wrms_error, safety_factor, exp_power
        integer :: i

        safety_factor = 0.9_DP
        exp_power = 0.125_DP  ! 1/8 for RKF78

        allocate(current_state(n))
        allocate(next_state_7th(n), next_state_8th(n))
        allocate(error_est(n), scale_vector(n))

        ! 输出缓冲区
        allocate(times(max_steps))
        allocate(states(max_steps, n))

        current_state = state
        current_time = t_start
        times(1) = current_time
        states(1, :) = current_state
        n_steps = 1

        dt = min(dt_max, (t_end - t_start) / 100.0_DP)

        do i = 2, max_steps
            if (current_time >= t_end) exit

            call rkf78_var_dim(current_state, n, dt, current_time, &
                                mu, order, n, next_state_7th, next_state_8th, error_est)

            ! WRMS 误差评估
            scale_vector = abs_tol + rel_tol * abs(current_state)
            wrms_error = sqrt(sum((error_est / scale_vector)**2) / real(n, DP))

            if (ieee_is_nan(wrms_error)) then
                write(*,*) '[STT] FATAL: WRMS error is NaN, aborting integration'
                exit
            end if

            if (wrms_error <= 1.0_DP) then
                ! 接受此步
                current_state = next_state_8th
                current_time = current_time + dt
                n_steps = n_steps + 1
                times(n_steps) = current_time
                states(n_steps, :) = current_state

                dt = safety_factor * dt * (1.0_DP / max(wrms_error, 1.0e-15_DP))**exp_power
                dt = max(dt_min, min(dt_max, dt))
            else
                if (dt <= dt_min) then
                    current_state = next_state_8th
                    current_time = current_time + dt
                    n_steps = n_steps + 1
                    times(n_steps) = current_time
                    states(n_steps, :) = current_state
                else
                    dt = safety_factor * dt * (1.0_DP / max(wrms_error, 1.0e-15_DP))**exp_power
                    dt = max(dt_min, dt)
                end if
            end if

            if (current_time + dt > t_end) then
                dt = t_end - current_time
            end if
        end do

        ! 截断输出数组
        if (n_steps < max_steps) then
            block
                real(DP), allocatable :: temp_times(:), temp_states(:,:)
                allocate(temp_times(n_steps), temp_states(n_steps, n))
                temp_times = times(1:n_steps)
                temp_states = states(1:n_steps, :)
                call move_alloc(temp_times, times)
                call move_alloc(temp_states, states)
            end block
        end if

        deallocate(current_state, next_state_7th, next_state_8th)
        deallocate(error_est, scale_vector)
    end subroutine adaptive_integrate_var_dim

    ! =====================================================================
    ! 使用给定初始偏移量传播 (deviates 版本)
    !
    ! 不依赖协方差矩阵，直接用 STT Taylor 多项式映射每个偏离量。
    ! =====================================================================
    subroutine stt_propagate_deviates(this, t_start, t_end, &
                                       nominal_state, deviates, &
                                       output_samples, output_mean, output_cov, &
                                       stt_store_out)
        class(uq_stt_propagator), intent(inout) :: this
        real(DP), intent(in) :: t_start, t_end
        real(DP), intent(in) :: nominal_state(6)
        real(DP), intent(in) :: deviates(:,:)
        real(DP), allocatable, intent(out) :: output_samples(:,:)
        real(DP), intent(out) :: output_mean(6)
        real(DP), intent(out) :: output_cov(6,6)
        type(stt_store_type), intent(out), optional :: stt_store_out

        real(DP), allocatable :: times(:), states(:,:)
        integer :: n_steps, p, np, pos, i, j, n_dev
        type(stt_store_type) :: stt_store
        real(DP) :: x_star(6)

        n_dev = size(deviates, 2)

        if (this%stt_order < 2 .or. this%stt_order > STT_MAX_ORDER) then
            write(*,*) '[STT] ERROR: call set_stt_order() before propagate()'
            return
        end if

        call init_stt_indexing(this%stt_order)
        call crtbp_derivatives_init(this%mu)

        ! ---- 1. 计算增广维度 ----
        this%total_dim = 6
        do p = 1, this%stt_order
            this%total_dim = this%total_dim + 6 * stt_sizes(p)
        end do

        ! ---- 2. 初始化增广状态 ----
        allocate(this%aug_state(this%total_dim))
        this%aug_state = 0.0_DP
        this%aug_state(1:6) = nominal_state
        do i = 1, 6
            this%aug_state(7 + (i-1)*6 + (i-1)) = 1.0_DP
        end do

        ! ---- 3. 联合积分 ----
        call adaptive_integrate_var_dim(this%aug_state, this%total_dim, &
            t_start, t_end, &
            this%stt_abs_tol, this%stt_rel_tol, &
            this%stt_dt_min, this%stt_dt_max, this%stt_max_steps, &
            this%stt_order, this%mu, times, states, n_steps)

        ! ---- 4. 提取终点 STT ----
        x_star = states(n_steps, 1:6)
        call stt_store%init(this%stt_order)
        pos = 7
        do p = 1, this%stt_order
            np = stt_sizes(p)
            do i = 1, 6
                do j = 1, np
                    call stt_store%set(i, j, p, states(n_steps, pos))
                    pos = pos + 1
                end do
            end do
        end do

        ! ---- 5. 用 STT 多项式映射所有偏离量 ----
        allocate(output_samples(6, n_dev))
        call stt_eval_deviates(x_star, stt_store, this%stt_order, &
                                deviates, n_dev, output_samples)

        ! ---- 6. 计算经验矩 ----
        output_mean = 0.0_DP
        do i = 1, n_dev
            output_mean = output_mean + output_samples(:, i)
        end do
        output_mean = output_mean / real(n_dev, DP)

        output_cov = 0.0_DP
        do i = 1, n_dev
            do j = 1, 6
                output_cov(:, j) = output_cov(:, j) + &
                    (output_samples(:, i) - output_mean) * (output_samples(j, i) - output_mean(j))
            end do
        end do
        output_cov = output_cov / real(n_dev - 1, DP)

        if (present(stt_store_out)) then
            stt_store_out%order = this%stt_order
            if (allocated(stt_store_out%stt)) deallocate(stt_store_out%stt)
            allocate(stt_store_out%stt(6, size(stt_store%stt, 2), 0:this%stt_order))
            stt_store_out%stt = stt_store%stt
        end if

        call stt_store%destroy()
        deallocate(times, states, this%aug_state)
    end subroutine stt_propagate_deviates

    ! =====================================================================
    ! 用 STT Taylor 多项式映射偏离量样本
    !
    !   dx_f(i) = Σ_{p=1}^{order} Σ_{c} Φ^p_{i,c} * Π_j dx0(tup_c(j)) / Π(cnt_v!)
    ! =====================================================================
    subroutine stt_eval_deviates(x_star, stt_store, order, deviates, n_dev, mapped)
        real(DP), intent(in) :: x_star(6)
        type(stt_store_type), intent(in) :: stt_store
        integer, intent(in) :: order, n_dev
        real(DP), intent(in) :: deviates(6, n_dev)
        real(DP), intent(out) :: mapped(6, n_dev)

        integer :: p, cidx, np, i_comp, k, j, tup(6)
        real(DP) :: prod, phi_val
        integer, allocatable, save :: eq_factors(:,:)

        if (.not. allocated(eq_factors)) then
            allocate(eq_factors(2:STT_MAX_ORDER, maxval(stt_sizes(2:STT_MAX_ORDER))))
            do p = 2, STT_MAX_ORDER
                np = stt_sizes(p)
                do cidx = 1, np
                    call sym_index_to_tuple(cidx, p, tup(1:p))
                    eq_factors(p, cidx) = compute_eq_factor(tup, p)
                end do
            end do
        end if

        do k = 1, n_dev
            mapped(:, k) = x_star

            ! Order 1: STM 直接矩阵乘法 (p=1: factor=1/1=1)
            do i_comp = 1, 6
                do j = 1, 6
                    mapped(i_comp, k) = mapped(i_comp, k) + &
                        stt_store%get(i_comp, j, 1) * deviates(j, k)
                end do
            end do

            ! Orders 2+: dx_i += Σ_c Φ^p_{i,c} · Π dx0(tup) / Π cnt_v!
            ! (p! in Taylor expansion cancels with permutation multiplicity)
            do p = 2, order
                np = stt_sizes(p)
                do i_comp = 1, 6
                    do cidx = 1, np
                        phi_val = stt_store%get(i_comp, cidx, p)
                        if (phi_val == 0.0_DP) cycle
                        call sym_index_to_tuple(cidx, p, tup(1:p))
                        prod = product(deviates(tup(1:p), k))
                        mapped(i_comp, k) = mapped(i_comp, k) + &
                            phi_val * prod / real(eq_factors(p, cidx), DP)
                    end do
                end do
            end do
        end do
    end subroutine stt_eval_deviates

    ! =====================================================================
    ! 计算对称元组的等式因子: Π_v cnt(v)!
    ! 即: 每个不同值出现次数的阶乘之积
    ! =====================================================================
    integer function compute_eq_factor(tup, p) result(fac)
        integer, intent(in) :: tup(p), p
        integer :: cnt(6), i, j

        cnt = 0
        do i = 1, p
            cnt(tup(i)) = cnt(tup(i)) + 1
        end do

        fac = 1
        do i = 1, 6
            do j = 2, cnt(i)
                fac = fac * j
            end do
        end do
    end function compute_eq_factor

end module pod_uq_prop_stt_module
