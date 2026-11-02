!> Adaptive RKF78 integrator for the CRTBP + absolute-beta SRP study force
!! model (Q1-Q4). Mirrors pod_crtbp_integrator_module's RKF78 structure but
!! calls crtbp_srp_study_derivatives_{real,da} with the study config and
!! parameters. Parameters [beta,chi,u,v] are constant over one integration arc.
module pod_crtbp_srp_study_integrator_module
    use pod_global, only: DP
    use pod_dace_classes
    use pod_crtbp_srp_study_module, only: CrtbpSrpStudyConfig, SrpStudyTempPool, &
        crtbp_srp_study_derivatives_real, crtbp_srp_study_derivatives_da
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    implicit none
    private

    public :: StudyIntegratorTempPool
    public :: study_adaptive_integrate_real
    public :: study_adaptive_integrate_da

    !> Caller-owned DA temp pool: SRP force pool + integrator AlgebraicVector
    !! scratch. Reuse one pool across many segments; do not re-init mid-arc.
    type, public :: StudyIntegratorTempPool
        type(AlgebraicVector) :: tmp1, tmp2, tmp3, t4, t5, t6, t7, t8, t9, temp_state
        type(SrpStudyTempPool) :: force_pool
    contains
        procedure :: init => study_int_pool_init
        procedure :: destroy => study_int_pool_destroy
    end type

contains

    subroutine study_int_pool_init(this)
        class(StudyIntegratorTempPool), intent(inout) :: this
        call this%tmp1%init(6); call this%tmp2%init(6); call this%tmp3%init(6)
        call this%t4%init(6); call this%t5%init(6); call this%t6%init(6)
        call this%t7%init(6); call this%t8%init(6); call this%t9%init(6)
        call this%temp_state%init(6)
        call this%force_pool%init()
    end subroutine study_int_pool_init

    subroutine study_int_pool_destroy(this)
        class(StudyIntegratorTempPool), intent(inout) :: this
        call this%tmp1%destroy(); call this%tmp2%destroy(); call this%tmp3%destroy()
        call this%t4%destroy(); call this%t5%destroy(); call this%t6%destroy()
        call this%t7%destroy(); call this%t8%destroy(); call this%t9%destroy()
        call this%temp_state%destroy()
        call this%force_pool%destroy()
    end subroutine study_int_pool_destroy

    ! ========================================================================
    !  REAL RKF78
    ! ========================================================================
    subroutine study_rkf78_real_step(state, dt, t, state_8th, state_7th, error_est, cfg, p)
        real(DP), dimension(6), intent(in)  :: state
        real(DP), intent(in)  :: dt, t
        real(DP), dimension(6), intent(out) :: state_8th, state_7th
        real(DP), dimension(6), intent(out) :: error_est
        type(CrtbpSrpStudyConfig), intent(in) :: cfg
        real(DP), intent(in) :: p(4)

        real(DP), dimension(6) :: f0, f1, f2, f3, f4, f5, f6
        real(DP), dimension(6) :: f7, f8, f9, f10, f11, f12
        real(DP), dimension(6) :: temp_state

        real(DP), parameter :: a1 = 2.0_DP/27.0_DP,  a2 = 1.0_DP/9.0_DP
        real(DP), parameter :: a3 = 1.0_DP/6.0_DP,   a4 = 5.0_DP/12.0_DP
        real(DP), parameter :: a5 = 1.0_DP/2.0_DP,   a6 = 5.0_DP/6.0_DP
        real(DP), parameter :: a7 = 1.0_DP/6.0_DP,   a8 = 2.0_DP/3.0_DP
        real(DP), parameter :: a9 = 1.0_DP/3.0_DP

        real(DP), parameter :: b10 = 2.0_DP/27.0_DP
        real(DP), parameter :: b20 = 1.0_DP/36.0_DP,    b21 = 1.0_DP/12.0_DP
        real(DP), parameter :: b30 = 1.0_DP/24.0_DP,    b32 = 1.0_DP/8.0_DP
        real(DP), parameter :: b40 = 5.0_DP/12.0_DP,    b42 = -25.0_DP/16.0_DP, &
                               b43 = 25.0_DP/16.0_DP
        real(DP), parameter :: b50 = 1.0_DP/20.0_DP,    b53 = 1.0_DP/4.0_DP,    &
                               b54 = 1.0_DP/5.0_DP
        real(DP), parameter :: b60 = -25.0_DP/108.0_DP, b63 = 125.0_DP/108.0_DP, &
                               b64 = -65.0_DP/27.0_DP,  b65 = 125.0_DP/54.0_DP
        real(DP), parameter :: b70 = 31.0_DP/300.0_DP,  b74 = 61.0_DP/225.0_DP, &
                               b75 = -2.0_DP/9.0_DP,    b76 = 13.0_DP/900.0_DP
        real(DP), parameter :: b80 = 2.0_DP,            b83 = -53.0_DP/6.0_DP,   &
                               b84 = 704.0_DP/45.0_DP,  b85 = -107.0_DP/9.0_DP, &
                               b86 = 67.0_DP/90.0_DP,   b87 = 3.0_DP
        real(DP), parameter :: b90 = -91.0_DP/108.0_DP, b93 = 23.0_DP/108.0_DP, &
                               b94 = -976.0_DP/135.0_DP,b95 = 311.0_DP/54.0_DP, &
                               b96 = -19.0_DP/60.0_DP,  b97 = 17.0_DP/6.0_DP,   &
                               b98 = -1.0_DP/12.0_DP
        real(DP), parameter :: b100 = 2383.0_DP/4100.0_DP, &
                               b103 = -341.0_DP/164.0_DP,  &
                               b104 = 4496.0_DP/1025.0_DP, &
                               b105 = -301.0_DP/82.0_DP,   &
                               b106 = 2133.0_DP/4100.0_DP, &
                               b107 = 45.0_DP/82.0_DP,     &
                               b108 = 45.0_DP/164.0_DP,    &
                               b109 = 18.0_DP/41.0_DP
        real(DP), parameter :: b110 = 3.0_DP/205.0_DP,   b115 = -6.0_DP/41.0_DP, &
                               b116 = -3.0_DP/205.0_DP,  b117 = -3.0_DP/41.0_DP, &
                               b118 = 3.0_DP/41.0_DP,    b119 = 6.0_DP/41.0_DP
        real(DP), parameter :: b120 = -1777.0_DP/4100.0_DP, &
                               b123 = -341.0_DP/164.0_DP,  &
                               b124 = 4496.0_DP/1025.0_DP,  &
                               b125 = -289.0_DP/82.0_DP,   &
                               b126 = 2193.0_DP/4100.0_DP,  &
                               b127 = 51.0_DP/82.0_DP,     &
                               b128 = 33.0_DP/164.0_DP,    &
                               b129 = 12.0_DP/41.0_DP

        real(DP), parameter :: c5 = 34.0_DP/105.0_DP,  c6 = 9.0_DP/35.0_DP
        real(DP), parameter :: c7 = 9.0_DP/35.0_DP,    c8 = 9.0_DP/280.0_DP
        real(DP), parameter :: c9 = 9.0_DP/280.0_DP,   c11 = 41.0_DP/840.0_DP
        real(DP), parameter :: c12 = 41.0_DP/840.0_DP
        real(DP), parameter :: err_factor = 41.0_DP/840.0_DP

        if (dt <= 1.0e-15_DP) then
            state_8th = state
            state_7th = state
            error_est = 0.0_DP
            return
        end if

        call crtbp_srp_study_derivatives_real(state, f0, t, cfg, p)
        temp_state = state + dt * b10 * f0
        call crtbp_srp_study_derivatives_real(temp_state, f1, t + a1*dt, cfg, p)
        temp_state = state + dt * (b20*f0 + b21*f1)
        call crtbp_srp_study_derivatives_real(temp_state, f2, t + a2*dt, cfg, p)
        temp_state = state + dt * (b30*f0 + b32*f2)
        call crtbp_srp_study_derivatives_real(temp_state, f3, t + a3*dt, cfg, p)
        temp_state = state + dt * (b40*f0 + b42*f2 + b43*f3)
        call crtbp_srp_study_derivatives_real(temp_state, f4, t + a4*dt, cfg, p)
        temp_state = state + dt * (b50*f0 + b53*f3 + b54*f4)
        call crtbp_srp_study_derivatives_real(temp_state, f5, t + a5*dt, cfg, p)
        temp_state = state + dt * (b60*f0 + b63*f3 + b64*f4 + b65*f5)
        call crtbp_srp_study_derivatives_real(temp_state, f6, t + a6*dt, cfg, p)
        temp_state = state + dt * (b70*f0 + b74*f4 + b75*f5 + b76*f6)
        call crtbp_srp_study_derivatives_real(temp_state, f7, t + a7*dt, cfg, p)
        temp_state = state + dt * (b80*f0 + b83*f3 + b84*f4 + b85*f5 + b86*f6 + b87*f7)
        call crtbp_srp_study_derivatives_real(temp_state, f8, t + a8*dt, cfg, p)
        temp_state = state + dt * (b90*f0 + b93*f3 + b94*f4 + b95*f5 + b96*f6 + b97*f7 + b98*f8)
        call crtbp_srp_study_derivatives_real(temp_state, f9, t + a9*dt, cfg, p)
        temp_state = state + dt * (b100*f0 + b103*f3 + b104*f4 + b105*f5 &
            + b106*f6 + b107*f7 + b108*f8 + b109*f9)
        call crtbp_srp_study_derivatives_real(temp_state, f10, t + dt, cfg, p)
        temp_state = state + dt * (b110*f0 + b115*f5 + b116*f6 + b117*f7 + b118*f8 + b119*f9)
        call crtbp_srp_study_derivatives_real(temp_state, f11, t, cfg, p)
        temp_state = state + dt * (b120*f0 + b123*f3 + b124*f4 + b125*f5 &
            + b126*f6 + b127*f7 + b128*f8 + b129*f9 + f11)
        call crtbp_srp_study_derivatives_real(temp_state, f12, t + dt, cfg, p)

        state_8th = state + dt * (c5*f5 + c6*f6 + c7*f7 + c8*f8 + c9*f9 + c11*f11 + c12*f12)
        state_7th = state_8th - dt * err_factor * (f0 + f10 - f11 - f12)
        error_est = abs(state_8th - state_7th)
    end subroutine study_rkf78_real_step

    subroutine study_adaptive_integrate_real(state, t_start, t_end, &
            rel_tol, abs_tol, dt_min, dt_max, max_steps, cfg, p, final_state)
        real(DP), dimension(6), intent(in) :: state
        real(DP), intent(in) :: t_start, t_end
        real(DP), intent(in) :: rel_tol, abs_tol, dt_min, dt_max
        integer,  intent(in) :: max_steps
        type(CrtbpSrpStudyConfig), intent(in) :: cfg
        real(DP), intent(in) :: p(4)
        real(DP), dimension(6), intent(out) :: final_state

        real(DP), dimension(6) :: current_state, next_8th, next_7th
        real(DP), dimension(6) :: error_vector, scale_vector
        real(DP) :: current_time, dt, wrms_error
        real(DP), parameter :: safety_factor = 0.9_DP
        real(DP), parameter :: exp_power = 0.125_DP
        integer :: i

        current_state = state
        current_time = t_start
        dt = min(dt_max, (t_end - t_start) / 100.0_DP)

        do i = 2, max_steps
            if (current_time >= t_end) exit
            if (current_time + dt > t_end) dt = t_end - current_time

            call study_rkf78_real_step(current_state, dt, current_time, &
                next_8th, next_7th, error_vector, cfg, p)

            scale_vector = abs_tol + rel_tol * abs(current_state)
            wrms_error = sqrt(sum((error_vector / scale_vector)**2) / 6.0_DP)

            if (ieee_is_nan(wrms_error)) then
                final_state = current_state
                return
            end if

            if (wrms_error <= 1.0_DP) then
                current_state = next_8th
                current_time = current_time + dt
                dt = safety_factor * dt * (1.0_DP / max(wrms_error, 1.0e-15_DP))**exp_power
                dt = max(dt_min, min(dt_max, dt))
            else
                if (dt <= dt_min) then
                    current_state = next_8th
                    current_time = current_time + dt
                else
                    dt = safety_factor * dt * (1.0_DP / max(wrms_error, 1.0e-15_DP))**exp_power
                    dt = max(dt_min, dt)
                end if
            end if
        end do

        final_state = current_state
    end subroutine study_adaptive_integrate_real

    ! ========================================================================
    !  DA RKF78
    ! ========================================================================
    subroutine study_rkf78_da_step(state, dt, t, state_8th, state_7th, error_est, cfg, p, pool)
        type(AlgebraicVector), intent(in)    :: state
        real(DP), intent(in) :: dt, t
        type(AlgebraicVector), intent(inout) :: state_8th, state_7th
        type(AlgebraicVector), intent(inout) :: error_est
        type(CrtbpSrpStudyConfig), intent(in) :: cfg
        type(DA), intent(in) :: p(4)
        type(StudyIntegratorTempPool), intent(inout) :: pool

        type(AlgebraicVector) :: f0, f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12

        real(DP), parameter :: a1 = 2.0_DP/27.0_DP,  a2 = 1.0_DP/9.0_DP
        real(DP), parameter :: a3 = 1.0_DP/6.0_DP,   a4 = 5.0_DP/12.0_DP
        real(DP), parameter :: a5 = 1.0_DP/2.0_DP,   a6 = 5.0_DP/6.0_DP
        real(DP), parameter :: a7 = 1.0_DP/6.0_DP,   a8 = 2.0_DP/3.0_DP
        real(DP), parameter :: a9 = 1.0_DP/3.0_DP

        real(DP), parameter :: b10 = 2.0_DP/27.0_DP
        real(DP), parameter :: b20 = 1.0_DP/36.0_DP,    b21 = 1.0_DP/12.0_DP
        real(DP), parameter :: b30 = 1.0_DP/24.0_DP,    b32 = 1.0_DP/8.0_DP
        real(DP), parameter :: b40 = 5.0_DP/12.0_DP,    b42 = -25.0_DP/16.0_DP, &
                               b43 = 25.0_DP/16.0_DP
        real(DP), parameter :: b50 = 1.0_DP/20.0_DP,    b53 = 1.0_DP/4.0_DP,    &
                               b54 = 1.0_DP/5.0_DP
        real(DP), parameter :: b60 = -25.0_DP/108.0_DP, b63 = 125.0_DP/108.0_DP, &
                               b64 = -65.0_DP/27.0_DP,  b65 = 125.0_DP/54.0_DP
        real(DP), parameter :: b70 = 31.0_DP/300.0_DP,  b74 = 61.0_DP/225.0_DP, &
                               b75 = -2.0_DP/9.0_DP,    b76 = 13.0_DP/900.0_DP
        real(DP), parameter :: b80 = 2.0_DP,            b83 = -53.0_DP/6.0_DP,   &
                               b84 = 704.0_DP/45.0_DP,  b85 = -107.0_DP/9.0_DP, &
                               b86 = 67.0_DP/90.0_DP,   b87 = 3.0_DP
        real(DP), parameter :: b90 = -91.0_DP/108.0_DP, b93 = 23.0_DP/108.0_DP, &
                               b94 = -976.0_DP/135.0_DP,b95 = 311.0_DP/54.0_DP, &
                               b96 = -19.0_DP/60.0_DP,  b97 = 17.0_DP/6.0_DP,   &
                               b98 = -1.0_DP/12.0_DP
        real(DP), parameter :: b100 = 2383.0_DP/4100.0_DP, &
                               b103 = -341.0_DP/164.0_DP,  &
                               b104 = 4496.0_DP/1025.0_DP, &
                               b105 = -301.0_DP/82.0_DP,   &
                               b106 = 2133.0_DP/4100.0_DP, &
                               b107 = 45.0_DP/82.0_DP,     &
                               b108 = 45.0_DP/164.0_DP,    &
                               b109 = 18.0_DP/41.0_DP
        real(DP), parameter :: b110 = 3.0_DP/205.0_DP,   b115 = -6.0_DP/41.0_DP, &
                               b116 = -3.0_DP/205.0_DP,  b117 = -3.0_DP/41.0_DP, &
                               b118 = 3.0_DP/41.0_DP,    b119 = 6.0_DP/41.0_DP
        real(DP), parameter :: b120 = -1777.0_DP/4100.0_DP, &
                               b123 = -341.0_DP/164.0_DP,  &
                               b124 = 4496.0_DP/1025.0_DP,  &
                               b125 = -289.0_DP/82.0_DP,   &
                               b126 = 2193.0_DP/4100.0_DP,  &
                               b127 = 51.0_DP/82.0_DP,     &
                               b128 = 33.0_DP/164.0_DP,    &
                               b129 = 12.0_DP/41.0_DP

        real(DP), parameter :: c5 = 34.0_DP/105.0_DP,  c6 = 9.0_DP/35.0_DP
        real(DP), parameter :: c7 = 9.0_DP/35.0_DP,    c8 = 9.0_DP/280.0_DP
        real(DP), parameter :: c9 = 9.0_DP/280.0_DP,   c11 = 41.0_DP/840.0_DP
        real(DP), parameter :: c12 = 41.0_DP/840.0_DP
        real(DP), parameter :: err_factor = 41.0_DP/840.0_DP

        call f0%init(6); call f1%init(6); call f2%init(6); call f3%init(6)
        call f4%init(6); call f5%init(6); call f6%init(6); call f7%init(6)
        call f8%init(6); call f9%init(6); call f10%init(6); call f11%init(6)
        call f12%init(6)

        if (state_8th%size /= 6) call state_8th%init(6)
        if (state_7th%size /= 6) call state_7th%init(6)
        if (error_est%size /= 6) call error_est%init(6)

        if (dt <= 1.0e-15_DP) then
            state_8th = state
            state_7th = state
            error_est = 0.0_DP
            go to 999
        end if

        call crtbp_srp_study_derivatives_da(state%elements, f0%elements, t, cfg, p, pool%force_pool)

        call real_mul_vector_sub(dt*b10, f0, pool%tmp1)
        call vec_add(state, pool%tmp1, pool%temp_state)
        call crtbp_srp_study_derivatives_da(pool%temp_state%elements, f1%elements, t + a1*dt, cfg, p, pool%force_pool)

        call real_mul_vector_sub(dt*b20, f0, pool%tmp1)
        call real_mul_vector_sub(dt*b21, f1, pool%tmp2)
        call vec_add(pool%tmp1, pool%tmp2, pool%tmp3)
        call vec_add(state, pool%tmp3, pool%temp_state)
        call crtbp_srp_study_derivatives_da(pool%temp_state%elements, f2%elements, t + a2*dt, cfg, p, pool%force_pool)

        call real_mul_vector_sub(dt*b30, f0, pool%tmp1)
        call real_mul_vector_sub(dt*b32, f2, pool%tmp2)
        call vec_add(pool%tmp1, pool%tmp2, pool%tmp3)
        call vec_add(state, pool%tmp3, pool%temp_state)
        call crtbp_srp_study_derivatives_da(pool%temp_state%elements, f3%elements, t + a3*dt, cfg, p, pool%force_pool)

        call real_mul_vector_sub(dt*b40, f0, pool%tmp1)
        call real_mul_vector_sub(dt*b42, f2, pool%tmp2)
        call real_mul_vector_sub(dt*b43, f3, pool%tmp3)
        call vec_add(pool%tmp1, pool%tmp2, pool%t4)
        call vec_add(pool%t4, pool%tmp3, pool%t5)
        call vec_add(state, pool%t5, pool%temp_state)
        call crtbp_srp_study_derivatives_da(pool%temp_state%elements, f4%elements, t + a4*dt, cfg, p, pool%force_pool)

        call real_mul_vector_sub(dt*b50, f0, pool%tmp1)
        call real_mul_vector_sub(dt*b53, f3, pool%tmp2)
        call real_mul_vector_sub(dt*b54, f4, pool%tmp3)
        call vec_add(pool%tmp1, pool%tmp2, pool%t4)
        call vec_add(pool%t4, pool%tmp3, pool%t5)
        call vec_add(state, pool%t5, pool%temp_state)
        call crtbp_srp_study_derivatives_da(pool%temp_state%elements, f5%elements, t + a5*dt, cfg, p, pool%force_pool)

        call real_mul_vector_sub(dt*b60, f0, pool%tmp1)
        call real_mul_vector_sub(dt*b63, f3, pool%tmp2)
        call real_mul_vector_sub(dt*b64, f4, pool%tmp3)
        call real_mul_vector_sub(dt*b65, f5, pool%t4)
        call vec_add(pool%tmp1, pool%tmp2, pool%t5)
        call vec_add(pool%t5, pool%tmp3, pool%t6)
        call vec_add(pool%t6, pool%t4, pool%t7)
        call vec_add(state, pool%t7, pool%temp_state)
        call crtbp_srp_study_derivatives_da(pool%temp_state%elements, f6%elements, t + a6*dt, cfg, p, pool%force_pool)

        call real_mul_vector_sub(dt*b70, f0, pool%tmp1)
        call real_mul_vector_sub(dt*b74, f4, pool%tmp2)
        call real_mul_vector_sub(dt*b75, f5, pool%tmp3)
        call real_mul_vector_sub(dt*b76, f6, pool%t4)
        call vec_add(pool%tmp1, pool%tmp2, pool%t5)
        call vec_add(pool%t5, pool%tmp3, pool%t6)
        call vec_add(pool%t6, pool%t4, pool%t7)
        call vec_add(state, pool%t7, pool%temp_state)
        call crtbp_srp_study_derivatives_da(pool%temp_state%elements, f7%elements, t + a7*dt, cfg, p, pool%force_pool)

        call real_mul_vector_sub(dt*b80, f0, pool%tmp1)
        call real_mul_vector_sub(dt*b83, f3, pool%tmp2)
        call real_mul_vector_sub(dt*b84, f4, pool%tmp3)
        call real_mul_vector_sub(dt*b85, f5, pool%t4)
        call real_mul_vector_sub(dt*b86, f6, pool%t5)
        call real_mul_vector_sub(dt*b87, f7, pool%t6)
        call vec_add(pool%tmp1, pool%tmp2, pool%t7)
        call vec_add(pool%t7, pool%tmp3, pool%t8)
        call vec_add(pool%t8, pool%t4, pool%t7)
        call vec_add(pool%t7, pool%t5, pool%t8)
        call vec_add(pool%t8, pool%t6, pool%t7)
        call vec_add(state, pool%t7, pool%temp_state)
        call crtbp_srp_study_derivatives_da(pool%temp_state%elements, f8%elements, t + a8*dt, cfg, p, pool%force_pool)

        call real_mul_vector_sub(dt*b90, f0, pool%tmp1)
        call real_mul_vector_sub(dt*b93, f3, pool%tmp2)
        call real_mul_vector_sub(dt*b94, f4, pool%tmp3)
        call real_mul_vector_sub(dt*b95, f5, pool%t4)
        call real_mul_vector_sub(dt*b96, f6, pool%t5)
        call real_mul_vector_sub(dt*b97, f7, pool%t6)
        call real_mul_vector_sub(dt*b98, f8, pool%t7)
        call vec_add(pool%tmp1, pool%tmp2, pool%t8)
        call vec_add(pool%t8, pool%tmp3, pool%t9)
        call vec_add(pool%t9, pool%t4, pool%t8)
        call vec_add(pool%t8, pool%t5, pool%t9)
        call vec_add(pool%t9, pool%t6, pool%t8)
        call vec_add(pool%t8, pool%t7, pool%t9)
        call vec_add(pool%t9, state, pool%temp_state)
        call crtbp_srp_study_derivatives_da(pool%temp_state%elements, f9%elements, t + a9*dt, cfg, p, pool%force_pool)

        call real_mul_vector_sub(dt*b100, f0, pool%tmp1)
        call real_mul_vector_sub(dt*b103, f3, pool%tmp2)
        call real_mul_vector_sub(dt*b104, f4, pool%tmp3)
        call real_mul_vector_sub(dt*b105, f5, pool%t4)
        call real_mul_vector_sub(dt*b106, f6, pool%t5)
        call real_mul_vector_sub(dt*b107, f7, pool%t6)
        call real_mul_vector_sub(dt*b108, f8, pool%t7)
        call real_mul_vector_sub(dt*b109, f9, pool%t8)
        call vec_add(pool%tmp1, pool%tmp2, pool%t9)
        call vec_add(pool%t9, pool%tmp3, pool%tmp1)
        call vec_add(pool%tmp1, pool%t4, pool%t9)
        call vec_add(pool%t9, pool%t5, pool%tmp1)
        call vec_add(pool%tmp1, pool%t6, pool%t9)
        call vec_add(pool%t9, pool%t7, pool%tmp1)
        call vec_add(pool%tmp1, pool%t8, pool%t9)
        call vec_add(state, pool%t9, pool%temp_state)
        call crtbp_srp_study_derivatives_da(pool%temp_state%elements, f10%elements, t + dt, cfg, p, pool%force_pool)

        call real_mul_vector_sub(dt*b110, f0, pool%tmp1)
        call real_mul_vector_sub(dt*b115, f5, pool%tmp2)
        call real_mul_vector_sub(dt*b116, f6, pool%tmp3)
        call real_mul_vector_sub(dt*b117, f7, pool%t4)
        call real_mul_vector_sub(dt*b118, f8, pool%t5)
        call real_mul_vector_sub(dt*b119, f9, pool%t6)
        call vec_add(pool%tmp1, pool%tmp2, pool%t7)
        call vec_add(pool%t7, pool%tmp3, pool%t8)
        call vec_add(pool%t8, pool%t4, pool%t7)
        call vec_add(pool%t7, pool%t5, pool%t8)
        call vec_add(pool%t8, pool%t6, pool%t7)
        call vec_add(state, pool%t7, pool%temp_state)
        call crtbp_srp_study_derivatives_da(pool%temp_state%elements, f11%elements, t, cfg, p, pool%force_pool)

        call real_mul_vector_sub(dt*b120, f0, pool%tmp1)
        call real_mul_vector_sub(dt*b123, f3, pool%tmp2)
        call real_mul_vector_sub(dt*b124, f4, pool%tmp3)
        call real_mul_vector_sub(dt*b125, f5, pool%t4)
        call real_mul_vector_sub(dt*b126, f6, pool%t5)
        call real_mul_vector_sub(dt*b127, f7, pool%t6)
        call real_mul_vector_sub(dt*b128, f8, pool%t7)
        call real_mul_vector_sub(dt*b129, f9, pool%t8)
        call vec_add(pool%tmp1, pool%tmp2, pool%t9)
        call vec_add(pool%t9, pool%tmp3, pool%tmp1)
        call vec_add(pool%tmp1, pool%t4, pool%tmp2)
        call vec_add(pool%tmp2, pool%t5, pool%tmp1)
        call vec_add(pool%tmp1, pool%t6, pool%tmp2)
        call vec_add(pool%tmp2, pool%t7, pool%tmp1)
        call vec_add(pool%tmp1, pool%t8, pool%tmp2)
        call real_mul_vector_sub(dt, f11, pool%tmp1)
        call vec_add(pool%tmp1, pool%tmp2, pool%t8)
        call vec_add(state, pool%t8, pool%temp_state)
        call crtbp_srp_study_derivatives_da(pool%temp_state%elements, f12%elements, t + dt, cfg, p, pool%force_pool)

        call real_mul_vector_sub(dt*c5, f5, pool%tmp1)
        call real_mul_vector_sub(dt*c6, f6, pool%tmp2)
        call real_mul_vector_sub(dt*c7, f7, pool%tmp3)
        call real_mul_vector_sub(dt*c8, f8, pool%t4)
        call real_mul_vector_sub(dt*c9, f9, pool%t5)
        call real_mul_vector_sub(dt*c11, f11, pool%t6)
        call real_mul_vector_sub(dt*c12, f12, pool%t7)
        call vec_add(pool%tmp1, pool%tmp2, pool%t8)
        call vec_add(pool%t8, pool%tmp3, pool%t9)
        call vec_add(pool%t9, pool%t4, pool%t8)
        call vec_add(pool%t8, pool%t5, pool%t9)
        call vec_add(pool%t9, pool%t6, pool%t8)
        call vec_add(pool%t8, pool%t7, pool%t9)
        call vec_add(state, pool%t9, state_8th)

        call vec_add(f0, f10, pool%tmp1)
        call vec_sub(pool%tmp1, f11, pool%tmp2)
        call vec_sub(pool%tmp2, f12, pool%tmp3)
        call real_mul_vector_sub(dt * err_factor, pool%tmp3, pool%t4)
        call vec_sub(state_8th, pool%t4, state_7th)

        call vec_sub(state_8th, state_7th, error_est)

        999 continue
        call f0%destroy(); call f1%destroy(); call f2%destroy(); call f3%destroy()
        call f4%destroy(); call f5%destroy(); call f6%destroy(); call f7%destroy()
        call f8%destroy(); call f9%destroy(); call f10%destroy(); call f11%destroy()
        call f12%destroy()
    end subroutine study_rkf78_da_step

    subroutine study_adaptive_integrate_da(state, t_start, t_end, &
            rel_tol, abs_tol, dt_min, dt_max, max_steps, cfg, p, pool, final_state)
        type(AlgebraicVector), intent(in)    :: state
        real(DP), intent(in) :: t_start, t_end
        real(DP), intent(in) :: rel_tol, abs_tol, dt_min, dt_max
        integer,  intent(in) :: max_steps
        type(CrtbpSrpStudyConfig), intent(in) :: cfg
        type(DA), intent(in) :: p(4)
        type(StudyIntegratorTempPool), intent(inout) :: pool
        type(AlgebraicVector), intent(inout) :: final_state

        type(AlgebraicVector) :: current_state, next_8th, next_7th, error_vec
        real(DP) :: current_time, dt, wrms_error
        real(DP), dimension(6) :: scale_vector, cons_current
        integer :: i, j

        real(DP), parameter :: safety_factor = 0.9_DP
        real(DP), parameter :: exp_power = 0.125_DP

        call current_state%init(6)
        call next_8th%init(6)
        call next_7th%init(6)
        call error_vec%init(6)
        if (final_state%size /= 6) call final_state%init(6)

        do j = 1, 6
            current_state%elements(j) = state%elements(j)
        end do
        current_time = t_start
        dt = min(dt_max, (t_end - t_start) / 100.0_DP)

        do i = 2, max_steps
            if (current_time >= t_end) exit
            if (current_time + dt > t_end) dt = t_end - current_time

            call study_rkf78_da_step(current_state, dt, current_time, &
                next_8th, next_7th, error_vec, cfg, p, pool)

            cons_current = current_state%cons()
            do j = 1, 6
                scale_vector(j) = abs_tol + rel_tol * abs(cons_current(j))
            end do

            wrms_error = 0.0_DP
            do j = 1, 6
                wrms_error = wrms_error + (error_vec%elements(j)%cons() / scale_vector(j))**2
            end do
            wrms_error = sqrt(wrms_error / 6.0_DP)

            if (ieee_is_nan(wrms_error)) then
                do j = 1, 6
                    final_state%elements(j) = current_state%elements(j)
                end do
                call current_state%destroy(); call next_8th%destroy()
                call next_7th%destroy(); call error_vec%destroy()
                return
            end if

            if (wrms_error <= 1.0_DP) then
                do j = 1, 6
                    current_state%elements(j) = next_8th%elements(j)
                end do
                current_time = current_time + dt
                dt = safety_factor * dt * (1.0_DP / max(wrms_error, 1.0e-15_DP))**exp_power
                dt = max(dt_min, min(dt_max, dt))
            else
                if (dt <= dt_min) then
                    do j = 1, 6
                        current_state%elements(j) = next_8th%elements(j)
                    end do
                    current_time = current_time + dt
                else
                    dt = safety_factor * dt * (1.0_DP / max(wrms_error, 1.0e-15_DP))**exp_power
                    dt = max(dt_min, dt)
                end if
            end if
        end do

        do j = 1, 6
            final_state%elements(j) = current_state%elements(j)
        end do

        call current_state%destroy()
        call next_8th%destroy(); call next_7th%destroy()
        call error_vec%destroy()
    end subroutine study_adaptive_integrate_da

end module pod_crtbp_srp_study_integrator_module
