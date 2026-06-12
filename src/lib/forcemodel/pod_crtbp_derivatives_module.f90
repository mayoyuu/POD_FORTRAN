!--------------------------------------------------------------------------------------------------------------
!> CRTBP 解析力场导数模块
!>
!> 计算 CRTBP 力场沿标称轨道的 1~m 阶偏导数 f*_{i, k₁...kₚ}。
!>
!> 核心递推:
!>   U_n = 1/r^(2n+1)
!>   ∂U_n/∂q = -(2n+1)·(q-q₀)·U_{n+1}
!>
!> 位置偏导数使用项表法计算: 每个项 (a,b,c,n,coeff) 表示
!>   coeff · dx^a dy^b dz^c · U_n
!>
!> 依赖:
!>   - pod_global: DP 精度类型
!>   - pod_stt_tensor_module: 对称索引、stt_sizes
!>   - pod_crtbp_module: crtbp_derivatives_real (标称轨道)
!--------------------------------------------------------------------------------------------------------------
module pod_crtbp_derivatives_module
    use pod_global, only: DP
    use pod_stt_tensor_module, only: STT_MAX_ORDER, STT_DIM, stt_sizes, &
                                      tuple_to_sym_index, sym_index_to_tuple, &
                                      init_stt_indexing
    implicit none
    private

    public :: crtbp_force_derivatives, crtbp_derivatives_init

    ! ---- 模块级存储 ----
    real(DP), save :: mu_stored = 0.012153614091892_DP

    ! ---- 项表类型: coeff * dx^a dy^b dz^c * U_n ----
    integer, parameter :: MAX_TERMS = 30
    type :: term_list_type
        integer :: n_terms = 0
        integer :: a(MAX_TERMS) = 0
        integer :: b(MAX_TERMS) = 0
        integer :: c(MAX_TERMS) = 0
        integer :: n(MAX_TERMS) = 0
        real(DP) :: coeff(MAX_TERMS) = 0.0_DP
    end type

    ! ---- (项表类型定义已在上面) ----

contains

    ! =========================================================================
    ! 初始化 (设置 mu)
    ! =========================================================================
    subroutine crtbp_derivatives_init(mu)
        real(DP), intent(in) :: mu
        mu_stored = mu
    end subroutine crtbp_derivatives_init

    ! =========================================================================
    ! 主入口: 计算 CRTBP 力场解析导数
    !
    ! 输入: x(6) — 标称状态 (位置+速度)
    !       max_order — 最高导数阶数 (1..6)
    ! 输出: f_tensors(6, stt_sizes(max_order), 0:max_order)
    !       f_tensors(:, :, 0) — 标称力 (f₁..f₆)
    !       f_tensors(:, :, 1) — Jacobian (6×6 STM)
    !       f_tensors(:, :, p) — p 阶导数压缩格式
    ! =========================================================================
    subroutine crtbp_force_derivatives(x, max_order, f_tensors)
        real(DP), intent(in) :: x(6)
        integer, intent(in) :: max_order
        real(DP), allocatable, intent(out) :: f_tensors(:,:,:)

        integer :: p, idx, max_sz, i_comp
        integer :: tup(6)
        real(DP) :: r1, r2, dx1, dy1, dz1, dx2, dy2, dz2
        real(DP) :: U1(0:STT_MAX_ORDER+1), U2(0:STT_MAX_ORDER+1)

        call init_stt_indexing(max_order)

        max_sz = stt_sizes(max_order)
        allocate(f_tensors(6, max_sz, 0:max_order))
        f_tensors = 0.0_DP

        ! ---- 标称力 (order 0) ----
        call crtbp_force_real(x, mu_stored, f_tensors(:, 1, 0))

        ! ---- 位置偏移 ----
        dx1 = x(1) + mu_stored
        dy1 = x(2)
        dz1 = x(3)
        dx2 = x(1) - 1.0_DP + mu_stored
        dy2 = x(2)
        dz2 = x(3)

        r1 = sqrt(dx1**2 + dy1**2 + dz1**2)
        r2 = sqrt(dx2**2 + dy2**2 + dz2**2)

        ! ---- 预计算 U_n = 1/r^(2n+1) ----
        call compute_u_series(r1, max_order + 1, U1)
        call compute_u_series(r2, max_order + 1, U2)

        ! ---- 填充各阶力场导数 ----
        do p = 1, max_order
            do idx = 1, stt_sizes(p)
                call sym_index_to_tuple(idx, p, tup(1:p))
                do i_comp = 1, 6
                    f_tensors(i_comp, idx, p) = &
                        eval_force_derivative(i_comp, p, tup(1:p), &
                                              dx1, dy1, dz1, dx2, dy2, dz2, &
                                              U1, U2)
                end do
            end do
        end do
    end subroutine crtbp_force_derivatives

    ! =========================================================================
    ! 计算标称 CRTBP 力 (不依赖 pod_crtbp_module 以避免循环依赖)
    ! =========================================================================
    subroutine crtbp_force_real(x, mu, f)
        real(DP), intent(in) :: x(6), mu
        real(DP), intent(out) :: f(6)
        real(DP) :: r1, r2, r13, r23
        real(DP) :: omx, omy, omz

        r1 = sqrt((x(1) + mu)**2 + x(2)**2 + x(3)**2)
        r2 = sqrt((x(1) - 1.0_DP + mu)**2 + x(2)**2 + x(3)**2)
        r13 = r1**3
        r23 = r2**3

        omx = x(1) - (1.0_DP - mu)*(x(1) + mu)/r13 - mu*(x(1) - 1.0_DP + mu)/r23
        omy = x(2) - (1.0_DP - mu)*x(2)/r13 - mu*x(2)/r23
        omz = 0.0_DP - (1.0_DP - mu)*x(3)/r13 - mu*x(3)/r23

        f(1) = x(4)
        f(2) = x(5)
        f(3) = x(6)
        f(4) = omx + 2.0_DP*x(5)
        f(5) = omy - 2.0_DP*x(4)
        f(6) = omz
    end subroutine crtbp_force_real

    ! =========================================================================
    ! 预计算 U_n = 1/r^(2n+1) 数列
    ! =========================================================================
    subroutine compute_u_series(r, max_n, U)
        real(DP), intent(in) :: r
        integer, intent(in) :: max_n
        real(DP), intent(out) :: U(0:max_n)
        integer :: i
        real(DP) :: r_inv

        r_inv = 1.0_DP / r
        U(0) = r_inv
        do i = 1, max_n
            U(i) = U(i-1) * r_inv * r_inv  ! U_n = U_{n-1} / r²
        end do
    end subroutine compute_u_series

    ! =========================================================================
    ! 对项表应用 ∂/∂x
    ! =========================================================================
    subroutine apply_dx(terms)
        type(term_list_type), intent(inout) :: terms
        integer :: t, n_old
        integer :: old_a(MAX_TERMS), old_b(MAX_TERMS), old_c(MAX_TERMS)
        integer :: old_n(MAX_TERMS)
        real(DP) :: old_coeff(MAX_TERMS)

        n_old = terms%n_terms
        old_a(1:n_old) = terms%a(1:n_old)
        old_b(1:n_old) = terms%b(1:n_old)
        old_c(1:n_old) = terms%c(1:n_old)
        old_n(1:n_old) = terms%n(1:n_old)
        old_coeff(1:n_old) = terms%coeff(1:n_old)

        terms%n_terms = 0
        do t = 1, n_old
            ! 降幂项: a·coeff * dx^{a-1} dy^b dz^c * U_n
            if (old_a(t) > 0) then
                call add_term(terms, old_a(t)-1, old_b(t), old_c(t), &
                              old_n(t), old_a(t) * old_coeff(t))
            end if
            ! 升幂项: -(2n+1)·coeff * dx^{a+1} dy^b dz^c * U_{n+1}
            call add_term(terms, old_a(t)+1, old_b(t), old_c(t), &
                          old_n(t)+1, -(2*old_n(t)+1) * old_coeff(t))
        end do
    end subroutine apply_dx

    subroutine apply_dy(terms)
        type(term_list_type), intent(inout) :: terms
        integer :: t, n_old
        integer :: old_a(MAX_TERMS), old_b(MAX_TERMS), old_c(MAX_TERMS)
        integer :: old_n(MAX_TERMS)
        real(DP) :: old_coeff(MAX_TERMS)

        n_old = terms%n_terms
        old_a(1:n_old) = terms%a(1:n_old)
        old_b(1:n_old) = terms%b(1:n_old)
        old_c(1:n_old) = terms%c(1:n_old)
        old_n(1:n_old) = terms%n(1:n_old)
        old_coeff(1:n_old) = terms%coeff(1:n_old)

        terms%n_terms = 0
        do t = 1, n_old
            if (old_b(t) > 0) then
                call add_term(terms, old_a(t), old_b(t)-1, old_c(t), &
                              old_n(t), old_b(t) * old_coeff(t))
            end if
            call add_term(terms, old_a(t), old_b(t)+1, old_c(t), &
                          old_n(t)+1, -(2*old_n(t)+1) * old_coeff(t))
        end do
    end subroutine apply_dy

    subroutine apply_dz(terms)
        type(term_list_type), intent(inout) :: terms
        integer :: t, n_old
        integer :: old_a(MAX_TERMS), old_b(MAX_TERMS), old_c(MAX_TERMS)
        integer :: old_n(MAX_TERMS)
        real(DP) :: old_coeff(MAX_TERMS)

        n_old = terms%n_terms
        old_a(1:n_old) = terms%a(1:n_old)
        old_b(1:n_old) = terms%b(1:n_old)
        old_c(1:n_old) = terms%c(1:n_old)
        old_n(1:n_old) = terms%n(1:n_old)
        old_coeff(1:n_old) = terms%coeff(1:n_old)

        terms%n_terms = 0
        do t = 1, n_old
            if (old_c(t) > 0) then
                call add_term(terms, old_a(t), old_b(t), old_c(t)-1, &
                              old_n(t), old_c(t) * old_coeff(t))
            end if
            call add_term(terms, old_a(t), old_b(t), old_c(t)+1, &
                          old_n(t)+1, -(2*old_n(t)+1) * old_coeff(t))
        end do
    end subroutine apply_dz

    ! =========================================================================
    ! 向项表添加一项 (合并同类项)
    ! =========================================================================
    subroutine add_term(terms, a, b, c, n, coeff)
        type(term_list_type), intent(inout) :: terms
        integer, intent(in) :: a, b, c, n
        real(DP), intent(in) :: coeff
        integer :: t

        if (abs(coeff) < 1.0d-30) return

        ! 查找是否已有相同的 (a,b,c,n)
        do t = 1, terms%n_terms
            if (terms%a(t) == a .and. terms%b(t) == b .and. &
                terms%c(t) == c .and. terms%n(t) == n) then
                terms%coeff(t) = terms%coeff(t) + coeff
                if (abs(terms%coeff(t)) < 1.0d-30) then
                    ! 系数归零, 移除该项
                    terms%a(t) = terms%a(terms%n_terms)
                    terms%b(t) = terms%b(terms%n_terms)
                    terms%c(t) = terms%c(terms%n_terms)
                    terms%n(t) = terms%n(terms%n_terms)
                    terms%coeff(t) = terms%coeff(terms%n_terms)
                    terms%n_terms = terms%n_terms - 1
                end if
                return
            end if
        end do

        ! 新项
        terms%n_terms = terms%n_terms + 1
        t = terms%n_terms
        terms%a(t) = a; terms%b(t) = b; terms%c(t) = c
        terms%n(t) = n; terms%coeff(t) = coeff
    end subroutine add_term

    ! =========================================================================
    ! 计算单个力场导数分量 f*_{i, tup}
    !
    ! CRTBP 力场: f₁=vx, f₂=vy, f₃=vz
    !              f₄=Ω_x+2vy, f₅=Ω_y-2vx, f₆=Ω_z
    !
    ! Ω = (x²+y²)/2 + (1-μ)/r₁ + μ/r₂
    !
    ! 速度导数的规则:
    !   - f₁ 仅 ∂/∂vx=1 (一阶), 其余为零
    !   - f₂ 仅 ∂/∂vy=1 (一阶), 其余为零
    !   - f₃ 仅 ∂/∂vz=1 (一阶), 其余为零
    !   - f₄ 仅 ∂/∂vy=2 (一阶), 其余速度导数为零
    !   - f₅ 仅 ∂/∂vx=-2 (一阶), 其余速度导数为零
    !   - f₆ 所有速度导数为零
    ! =========================================================================
    real(DP) function eval_force_derivative(i_comp, p, tup, &
                                             dx1, dy1, dz1, dx2, dy2, dz2, &
                                             U1, U2) result(val)
        integer, intent(in) :: i_comp, p
        integer, intent(in) :: tup(p)
        real(DP), intent(in) :: dx1, dy1, dz1, dx2, dy2, dz2
        real(DP), intent(in) :: U1(0:), U2(0:)

        integer :: pos_counts(3)  ! 各位置坐标出现的次数
        integer :: vel_counts(3)  ! 各速度坐标出现的次数
        integer :: i, comp

        val = 0.0_DP

        ! 统计位置 (1-3) 和速度 (4-6) 指标的出现次数
        pos_counts = 0
        vel_counts = 0
        do i = 1, p
            comp = tup(i)
            if (comp >= 1 .and. comp <= 3) then
                pos_counts(comp) = pos_counts(comp) + 1
            else if (comp >= 4 .and. comp <= 6) then
                vel_counts(comp - 3) = vel_counts(comp - 3) + 1
            end if
        end do

        ! ---- 速度导数规则 ----
        ! 如果总速度导数阶数 > 1，全为零
        if (sum(vel_counts) > 1) then
            val = 0.0_DP
            return
        end if

        ! 如果恰好一个速度导数
        if (sum(vel_counts) == 1) then
            ! p=1 且速度分量匹配时才有非零值
            if (p == 1) then
                select case (i_comp)
                case (1)
                    if (vel_counts(1) == 1) val = 1.0_DP  ! ∂f₁/∂vx
                case (2)
                    if (vel_counts(2) == 1) val = 1.0_DP  ! ∂f₂/∂vy
                case (3)
                    if (vel_counts(3) == 1) val = 1.0_DP  ! ∂f₃/∂vz
                case (4)
                    if (vel_counts(2) == 1) val = 2.0_DP  ! ∂f₄/∂vy
                case (5)
                    if (vel_counts(1) == 1) val = -2.0_DP ! ∂f₅/∂vx
                end select
            end if
            return
        end if

        ! ---- 纯位置导数 ----
        ! p = i+j+k, 其中 i=pos_counts(1), j=pos_counts(2), k=pos_counts(3)
        val = eval_position_force_derivative(i_comp, &
            pos_counts(1), pos_counts(2), pos_counts(3), &
            dx1, dy1, dz1, dx2, dy2, dz2, U1, U2)
    end function eval_force_derivative

    ! =========================================================================
    ! 计算纯位置力场导数 ∂^{i+j+k} f / ∂x^i ∂y^j ∂z^k
    ! =========================================================================
    real(DP) function eval_position_force_derivative(i_comp, i, j, k, &
                                                      dx1, dy1, dz1, dx2, dy2, dz2, &
                                                      U1, U2) result(val)
        integer, intent(in) :: i_comp, i, j, k
        real(DP), intent(in) :: dx1, dy1, dz1, dx2, dy2, dz2
        real(DP), intent(in) :: U1(0:), U2(0:)
        real(DP) :: omu, mu

        omu = 1.0_DP - mu_stored
        mu = mu_stored
        val = 0.0_DP

        select case (i_comp)
        case (1)  ! f₁ = vx — 纯位置导数为零
            val = 0.0_DP

        case (2)  ! f₂ = vy — 纯位置导数为零
            val = 0.0_DP

        case (3)  ! f₃ = vz — 纯位置导数为零
            val = 0.0_DP

        case (4)  ! f₄ = Ω_x + 2vy = x - omu*(x+μ)/r₁³ - mu*(x-1+μ)/r₂³ + 2vy
            ! 离心项 x: i==1,j==0,k==0 时 ∂/∂x(x)=1; 二阶以上为零
            if (i == 1 .and. j == 0 .and. k == 0) then
                val = val + 1.0_DP
            end if
            ! (1-μ)/r₁ 项的 x-导数: ∂/∂x((1-μ)/r₁) = -(1-μ)(x+μ)/r₁³
            ! 高阶导数: ∂^{i,j,k}/∂x^i∂y^j∂z^k [-(1-μ)·(x+μ)/r₁³]
            ! = -(1-μ) · ∂^{i,j,k}/∂x^i∂y^j∂z^k [(x+μ)/r₁³]
            ! (x+μ)/r₁³ = dx₁/r₁³ = dx₁·U₁(1)
            ! 对于高阶: dx₁·U₁(1)
            ! ∂/∂x[dx₁·U₁] = U₁ + dx₁·∂U₁/∂x = U₁ - 3·dx₁²·U₂
            ! 使用项表法: 从 dx₁·U_1 出发，应用剩余导数
            val = val + omu * eval_linear_r_term(i, j, k, dx1, dy1, dz1, U1, -1.0_DP, 1, 0, 0)

            ! μ/r₂ 项: 同理，注意符号
            val = val + mu * eval_linear_r_term(i, j, k, dx2, dy2, dz2, U2, -1.0_DP, 1, 0, 0)

        case (5)  ! f₅ = Ω_y - 2vx = y - omu*y/r₁³ - mu*y/r₂³ - 2vx
            ! 离心项 y
            if (i == 0 .and. j == 1 .and. k == 0) then
                val = val + 1.0_DP
            end if
            ! y/r₁³ = dy₁·U₁
            val = val + omu * eval_linear_r_term(i, j, k, dx1, dy1, dz1, U1, -1.0_DP, 0, 1, 0)
            val = val + mu * eval_linear_r_term(i, j, k, dx2, dy2, dz2, U2, -1.0_DP, 0, 1, 0)

        case (6)  ! f₆ = Ω_z = -omu*z/r₁³ - mu*z/r₂³
            ! z/r₁³ = dz₁·U₁
            val = val + omu * eval_linear_r_term(i, j, k, dx1, dy1, dz1, U1, -1.0_DP, 0, 0, 1)
            val = val + mu * eval_linear_r_term(i, j, k, dx2, dy2, dz2, U2, -1.0_DP, 0, 0, 1)
        end select
    end function eval_position_force_derivative

    ! =========================================================================
    ! 计算 ∂^{i,j,k}/∂x^i∂y^j∂z^k [dx * U_1] (用于 f₄ 中的 x 项)
    !
    ! 即: 从项 dx^1 dy^0 dz^0 U_1 出发，应用剩余 i-1 次 ∂/∂x, j 次 ∂/∂y, k 次 ∂/∂z
    ! =========================================================================
    real(DP) function eval_linear_r_term(i, j, k, dx, dy, dz, U, sign_factor, init_a, init_b, init_c) result(val)
        integer, intent(in) :: i, j, k
        real(DP), intent(in) :: dx, dy, dz
        real(DP), intent(in) :: U(0:)
        real(DP), intent(in) :: sign_factor
        integer, intent(in) :: init_a, init_b, init_c
        type(term_list_type) :: terms
        integer :: step, t
        real(DP) :: dx_pow, dy_pow, dz_pow

        ! 零阶: val = sign_factor * q * U_1
        if (i == 0 .and. j == 0 .and. k == 0) then
            select case (init_a + init_b*2 + init_c*4)
            case (1); val = sign_factor * dx * U(1)
            case (2); val = sign_factor * dy * U(1)
            case (4); val = sign_factor * dz * U(1)
            end select
            return
        end if

        ! 初始化项表: sign_factor * dx^init_a dy^init_b dz^init_c * U_1
        terms%n_terms = 1
        terms%a(1) = init_a; terms%b(1) = init_b; terms%c(1) = init_c
        terms%n(1) = 1; terms%coeff(1) = sign_factor

        do step = 1, i
            call apply_dx(terms)
        end do
        do step = 1, j
            call apply_dy(terms)
        end do
        do step = 1, k
            call apply_dz(terms)
        end do

        val = 0.0_DP
        do t = 1, terms%n_terms
            dx_pow = 1.0_DP; dy_pow = 1.0_DP; dz_pow = 1.0_DP
            if (terms%a(t) > 0) dx_pow = dx ** terms%a(t)
            if (terms%b(t) > 0) dy_pow = dy ** terms%b(t)
            if (terms%c(t) > 0) dz_pow = dz ** terms%c(t)
            val = val + terms%coeff(t) * dx_pow * dy_pow * dz_pow * U(terms%n(t))
        end do
    end function eval_linear_r_term

end module pod_crtbp_derivatives_module
