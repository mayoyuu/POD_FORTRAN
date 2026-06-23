module pod_crtbp_integrator_module
    use pod_global, only: DP
    use pod_crtbp_module
    use pod_dace_classes
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    implicit none
    private

    public :: rkf45_crtbp_step
    public :: adaptive_integrate_crtbp
    public :: da_rkf45_crtbp_step
    public :: da_adaptive_integrate_crtbp
    public :: rkf78_crtbp_step
    public :: adaptive_integrate_crtbp_rkf78
    public :: da_rkf78_crtbp_step
    public :: da_adaptive_integrate_crtbp_rkf78
    public :: CrtbpIntegratorTempPool

    type, public :: CrtbpIntegratorTempPool
        type(AlgebraicVector) :: k1, k2, k3, k4, k5, k6, k7
        type(AlgebraicVector) :: tmp1, tmp2, tmp3
        type(AlgebraicVector) :: t4, t5, t6, t7, t8, t9
        type(AlgebraicVector) :: temp_state
        type(CrtbpForceTempPool) :: force_pool
    contains
        procedure :: init    => crtbp_int_pool_init
        procedure :: destroy => crtbp_int_pool_destroy
    end type CrtbpIntegratorTempPool

contains

    ! ---- Temp Pool Init / Destroy ----

    subroutine crtbp_int_pool_init(this)
        class(CrtbpIntegratorTempPool), intent(inout) :: this
        call this%k1%init(6); call this%k2%init(6); call this%k3%init(6)
        call this%k4%init(6); call this%k5%init(6); call this%k6%init(6)
        call this%k7%init(6)
        call this%tmp1%init(6); call this%tmp2%init(6); call this%tmp3%init(6)
        call this%t4%init(6); call this%t5%init(6); call this%t6%init(6)
        call this%t7%init(6); call this%t8%init(6); call this%t9%init(6)
        call this%temp_state%init(6)
        call this%force_pool%init()
    end subroutine crtbp_int_pool_init

    subroutine crtbp_int_pool_destroy(this)
        class(CrtbpIntegratorTempPool), intent(inout) :: this
        call this%k1%destroy(); call this%k2%destroy(); call this%k3%destroy()
        call this%k4%destroy(); call this%k5%destroy(); call this%k6%destroy()
        call this%k7%destroy()
        call this%tmp1%destroy(); call this%tmp2%destroy(); call this%tmp3%destroy()
        call this%t4%destroy(); call this%t5%destroy(); call this%t6%destroy()
        call this%t7%destroy(); call this%t8%destroy(); call this%t9%destroy()
        call this%temp_state%destroy()
        call this%force_pool%destroy()
    end subroutine crtbp_int_pool_destroy

    ! ================================================================
    !  REAL RKF45
    ! ================================================================

    subroutine rkf45_crtbp_step(state, dt, t, state_5th, state_4th, error_est)
        real(DP), dimension(6), intent(in)  :: state
        real(DP), intent(in)  :: dt, t
        real(DP), dimension(6), intent(out) :: state_5th, state_4th
        real(DP), dimension(6), intent(out) :: error_est

        real(DP), dimension(6) :: k1, k2, k3, k4, k5, k6, k7
        real(DP), dimension(6) :: temp_state

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

        call crtbp_derivatives_real(state, k1, t)

        temp_state = state + dt * b21 * k1
        call crtbp_derivatives_real(temp_state, k2, t + a2 * dt)

        temp_state = state + dt * (b31 * k1 + b32 * k2)
        call crtbp_derivatives_real(temp_state, k3, t + a3 * dt)

        temp_state = state + dt * (b41 * k1 + b42 * k2 + b43 * k3)
        call crtbp_derivatives_real(temp_state, k4, t + a4 * dt)

        temp_state = state + dt * (b51 * k1 + b52 * k2 + b53 * k3 + b54 * k4)
        call crtbp_derivatives_real(temp_state, k5, t + a5 * dt)

        temp_state = state + dt * (b61 * k1 + b62 * k2 + b63 * k3 + b64 * k4 + b65 * k5)
        call crtbp_derivatives_real(temp_state, k6, t + dt)

        state_5th = state + dt * (c1 * k1 + c3 * k3 + c4 * k4 + c5 * k5 + c6 * k6)
        call crtbp_derivatives_real(state_5th, k7, t + dt)

        state_4th = state + dt * (d1 * k1 + d3 * k3 + d4 * k4 + d5 * k5 + d6 * k6 + d7 * k7)
        error_est = abs(state_5th - state_4th)
    end subroutine rkf45_crtbp_step

    subroutine adaptive_integrate_crtbp(state, t_start, t_end, &
            rel_tol, abs_tol, dt_min, dt_max, max_steps, &
            times, states, n_steps)
        real(DP), dimension(6), intent(in) :: state
        real(DP), intent(in) :: t_start, t_end
        real(DP), intent(in) :: rel_tol, abs_tol, dt_min, dt_max
        integer,  intent(in) :: max_steps
        real(DP), allocatable, dimension(:), intent(out) :: times
        real(DP), allocatable, dimension(:,:), intent(out) :: states
        integer, intent(out) :: n_steps

        real(DP), dimension(6) :: current_state, next_state_5th, next_state_4th
        real(DP), dimension(6) :: error_vector, scale_vector
        real(DP) :: current_time, dt, wrms_error
        real(DP), parameter :: safety_factor = 0.9_DP
        real(DP), parameter :: exp_power = 0.20_DP
        integer :: i

        if (allocated(times)) deallocate(times)
        if (allocated(states)) deallocate(states)
        allocate(times(max_steps))
        allocate(states(max_steps, 6))

        current_state = state
        current_time = t_start
        times(1) = current_time
        states(1, :) = current_state
        n_steps = 1
        dt = min(dt_max, (t_end - t_start) / 100.0_DP)

        do i = 2, max_steps
            if (current_time >= t_end) exit
            if (current_time + dt > t_end) dt = t_end - current_time

            call rkf45_crtbp_step(current_state, dt, current_time, &
                next_state_5th, next_state_4th, error_vector)

            scale_vector = abs_tol + rel_tol * abs(current_state)
            wrms_error = sqrt(sum((error_vector / scale_vector)**2) / 6.0_DP)

            if (ieee_is_nan(wrms_error)) then
                print *, 'FATAL: WRMS error is NaN, aborting CRTBP integration.'
                exit
            end if

            if (wrms_error <= 1.0_DP) then
                current_state = next_state_5th
                current_time = current_time + dt
                n_steps = n_steps + 1
                times(n_steps) = current_time
                states(n_steps, :) = current_state
                dt = safety_factor * dt * (1.0_DP / max(wrms_error, 1.0e-15_DP))**exp_power
                dt = max(dt_min, min(dt_max, dt))
            else
                if (dt <= dt_min) then
                    current_state = next_state_5th
                    current_time = current_time + dt
                    n_steps = n_steps + 1
                    times(n_steps) = current_time
                    states(n_steps, :) = current_state
                else
                    dt = safety_factor * dt * (1.0_DP / max(wrms_error, 1.0e-15_DP))**exp_power
                    dt = max(dt_min, dt)
                end if
            end if
        end do

        if (n_steps < max_steps) then
            block
                real(DP), allocatable, dimension(:) :: tmp_times
                real(DP), allocatable, dimension(:,:) :: tmp_states
                allocate(tmp_times(n_steps))
                allocate(tmp_states(n_steps, 6))
                tmp_times = times(1:n_steps)
                tmp_states = states(1:n_steps, :)
                call move_alloc(tmp_times, times)
                call move_alloc(tmp_states, states)
            end block
        end if
    end subroutine adaptive_integrate_crtbp

    ! ================================================================
    !  DA RKF45
    ! ================================================================

    subroutine da_rkf45_crtbp_step(state, dt, t, state_5th, state_4th, error_est, pool)
        type(AlgebraicVector), intent(in)    :: state
        real(DP), intent(in) :: dt, t
        type(AlgebraicVector), intent(inout) :: state_5th, state_4th
        type(AlgebraicVector), intent(inout) :: error_est
        type(CrtbpIntegratorTempPool), intent(inout) :: pool

        type(DA), dimension(6) :: d_state, d_temp, d_k1, d_k2, d_k3, d_k4, d_k5, d_k6, d_k7
        type(DA), dimension(6) :: d_5th, d_4th, d_err
        integer :: j

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

        ! Extract state into DA array
        do j = 1, 6
            call d_state(j)%init()
            d_state(j) = state%elements(j)
            call d_k1(j)%init(); call d_k2(j)%init(); call d_k3(j)%init()
            call d_k4(j)%init(); call d_k5(j)%init(); call d_k6(j)%init(); call d_k7(j)%init()
            call d_temp(j)%init()
            call d_5th(j)%init(); call d_4th(j)%init(); call d_err(j)%init()
        end do

        ! k1
        call crtbp_derivatives_da(d_state, d_k1, t, pool%force_pool)

        ! k2: temp = state + dt * b21 * k1
        do j = 1, 6
            d_temp(j) = d_state(j) + (dt * b21) * d_k1(j)
        end do
        call crtbp_derivatives_da(d_temp, d_k2, t + a2 * dt, pool%force_pool)

        ! k3
        do j = 1, 6
            d_temp(j) = d_state(j) + dt * (b31 * d_k1(j) + b32 * d_k2(j))
        end do
        call crtbp_derivatives_da(d_temp, d_k3, t + a3 * dt, pool%force_pool)

        ! k4
        do j = 1, 6
            d_temp(j) = d_state(j) + dt * (b41 * d_k1(j) + b42 * d_k2(j) + b43 * d_k3(j))
        end do
        call crtbp_derivatives_da(d_temp, d_k4, t + a4 * dt, pool%force_pool)

        ! k5
        do j = 1, 6
            d_temp(j) = d_state(j) + dt * (b51 * d_k1(j) + b52 * d_k2(j) + b53 * d_k3(j) + b54 * d_k4(j))
        end do
        call crtbp_derivatives_da(d_temp, d_k5, t + a5 * dt, pool%force_pool)

        ! k6
        do j = 1, 6
            d_temp(j) = d_state(j) + dt * (b61 * d_k1(j) + b62 * d_k2(j) + b63 * d_k3(j) + b64 * d_k4(j) + b65 * d_k5(j))
        end do
        call crtbp_derivatives_da(d_temp, d_k6, t + dt, pool%force_pool)

        ! state_5th
        do j = 1, 6
            d_5th(j) = d_state(j) + dt * (c1 * d_k1(j) + c3 * d_k3(j) + c4 * d_k4(j) + c5 * d_k5(j) + c6 * d_k6(j))
        end do
        call crtbp_derivatives_da(d_5th, d_k7, t + dt, pool%force_pool)

        ! state_4th
        do j = 1, 6
            d_4th(j) = d_state(j) + dt * (d1 * d_k1(j) + d3 * d_k3(j) + d4 * d_k4(j) + d5 * d_k5(j) + d6 * d_k6(j) + d7 * d_k7(j))
            d_err(j) = d_5th(j) - d_4th(j)
        end do

        ! Copy results back to AlgebraicVector
        do j = 1, 6
            state_5th%elements(j) = d_5th(j)
            state_4th%elements(j) = d_4th(j)
            error_est%elements(j)  = d_err(j)
        end do

        ! Cleanup local DA arrays
        do j = 1, 6
            call d_state(j)%destroy()
            call d_k1(j)%destroy(); call d_k2(j)%destroy(); call d_k3(j)%destroy()
            call d_k4(j)%destroy(); call d_k5(j)%destroy(); call d_k6(j)%destroy()
            call d_k7(j)%destroy()
            call d_temp(j)%destroy()
            call d_5th(j)%destroy(); call d_4th(j)%destroy(); call d_err(j)%destroy()
        end do
    end subroutine da_rkf45_crtbp_step

    subroutine da_adaptive_integrate_crtbp(state, t_start, t_end, &
            rel_tol, abs_tol, dt_min, dt_max, max_steps, &
            final_state)
        type(AlgebraicVector), intent(in)    :: state
        real(DP), intent(in) :: t_start, t_end
        real(DP), intent(in) :: rel_tol, abs_tol, dt_min, dt_max
        integer,  intent(in) :: max_steps
        type(AlgebraicVector), intent(inout) :: final_state

        type(CrtbpIntegratorTempPool) :: pool
        type(AlgebraicVector) :: current_state, next_5th, next_4th, error_vec
        real(DP) :: current_time, dt, wrms_error
        real(DP), dimension(6) :: scale_vector, cons_current
        integer :: i, j

        real(DP), parameter :: safety_factor = 0.9_DP
        real(DP), parameter :: exp_power = 0.20_DP

        call pool%init()
        call current_state%init(6)
        call next_5th%init(6)
        call next_4th%init(6)
        call error_vec%init(6)

        do j = 1, 6
            current_state%elements(j) = state%elements(j)
        end do
        current_time = t_start
        dt = min(dt_max, (t_end - t_start) / 100.0_DP)

        do i = 2, max_steps
            if (current_time >= t_end) exit
            if (current_time + dt > t_end) dt = t_end - current_time

            call da_rkf45_crtbp_step(current_state, dt, current_time, &
                next_5th, next_4th, error_vec, pool)

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
                print *, 'FATAL: DA WRMS error is NaN, aborting CRTBP integration.'
                exit
            end if

            if (wrms_error <= 1.0_DP) then
                ! Accept step: assign next_5th into current_state
                do j = 1, 6
                    current_state%elements(j) = next_5th%elements(j)
                end do
                current_time = current_time + dt
                dt = safety_factor * dt * (1.0_DP / max(wrms_error, 1.0e-15_DP))**exp_power
                dt = max(dt_min, min(dt_max, dt))
            else
                if (dt <= dt_min) then
                    do j = 1, 6
                        current_state%elements(j) = next_5th%elements(j)
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
        call next_5th%destroy(); call next_4th%destroy()
        call error_vec%destroy()
        call pool%destroy()
    end subroutine da_adaptive_integrate_crtbp

    ! ================================================================
    !  REAL RKF78
    ! ================================================================

    subroutine rkf78_crtbp_step(state, dt, t, state_8th, state_7th, error_est)
        real(DP), dimension(6), intent(in)  :: state
        real(DP), intent(in)  :: dt, t
        real(DP), dimension(6), intent(out) :: state_8th, state_7th
        real(DP), dimension(6), intent(out) :: error_est

        real(DP), dimension(6) :: f0, f1, f2, f3, f4, f5, f6
        real(DP), dimension(6) :: f7, f8, f9, f10, f11, f12
        real(DP), dimension(6) :: temp_state
        integer :: j

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

        ! f0
        call crtbp_derivatives_real(state, f0, t)

        ! f1: state + dt*b10*f0
        temp_state = state + dt * b10 * f0
        call crtbp_derivatives_real(temp_state, f1, t + a1*dt)

        ! f2: state + dt*(b20*f0 + b21*f1)
        temp_state = state + dt * (b20*f0 + b21*f1)
        call crtbp_derivatives_real(temp_state, f2, t + a2*dt)

        ! f3: state + dt*(b30*f0 + b32*f2)
        temp_state = state + dt * (b30*f0 + b32*f2)
        call crtbp_derivatives_real(temp_state, f3, t + a3*dt)

        ! f4: state + dt*(b40*f0 + b42*f2 + b43*f3)
        temp_state = state + dt * (b40*f0 + b42*f2 + b43*f3)
        call crtbp_derivatives_real(temp_state, f4, t + a4*dt)

        ! f5: state + dt*(b50*f0 + b53*f3 + b54*f4)
        temp_state = state + dt * (b50*f0 + b53*f3 + b54*f4)
        call crtbp_derivatives_real(temp_state, f5, t + a5*dt)

        ! f6: state + dt*(b60*f0 + b63*f3 + b64*f4 + b65*f5)
        temp_state = state + dt * (b60*f0 + b63*f3 + b64*f4 + b65*f5)
        call crtbp_derivatives_real(temp_state, f6, t + a6*dt)

        ! f7: state + dt*(b70*f0 + b74*f4 + b75*f5 + b76*f6)
        temp_state = state + dt * (b70*f0 + b74*f4 + b75*f5 + b76*f6)
        call crtbp_derivatives_real(temp_state, f7, t + a7*dt)

        ! f8: state + dt*(b80*f0 + b83*f3 + b84*f4 + b85*f5 + b86*f6 + b87*f7)
        temp_state = state + dt * (b80*f0 + b83*f3 + b84*f4 + b85*f5 + b86*f6 + b87*f7)
        call crtbp_derivatives_real(temp_state, f8, t + a8*dt)

        ! f9: state + dt*(b90*f0 + b93*f3 + b94*f4 + b95*f5 + b96*f6 + b97*f7 + b98*f8)
        temp_state = state + dt * (b90*f0 + b93*f3 + b94*f4 + b95*f5 + b96*f6 + b97*f7 + b98*f8)
        call crtbp_derivatives_real(temp_state, f9, t + a9*dt)

        ! f10: state + dt*(b100*f0 + b103*f3 + b104*f4 + b105*f5 + b106*f6 + b107*f7 + b108*f8 + b109*f9)
        temp_state = state + dt * (b100*f0 + b103*f3 + b104*f4 + b105*f5 &
            + b106*f6 + b107*f7 + b108*f8 + b109*f9)
        call crtbp_derivatives_real(temp_state, f10, t + dt)

        ! f11: state + dt*(b110*f0 + b115*f5 + b116*f6 + b117*f7 + b118*f8 + b119*f9)
        temp_state = state + dt * (b110*f0 + b115*f5 + b116*f6 + b117*f7 + b118*f8 + b119*f9)
        call crtbp_derivatives_real(temp_state, f11, t)

        ! f12: state + dt*(b120*f0 + b123*f3 + b124*f4 + b125*f5 + b126*f6 + b127*f7 + b128*f8 + b129*f9 + f11)
        temp_state = state + dt * (b120*f0 + b123*f3 + b124*f4 + b125*f5 &
            + b126*f6 + b127*f7 + b128*f8 + b129*f9 + f11)
        call crtbp_derivatives_real(temp_state, f12, t + dt)

        ! 8th-order solution
        state_8th = state + dt * (c5*f5 + c6*f6 + c7*f7 + c8*f8 + c9*f9 + c11*f11 + c12*f12)

        ! 7th-order solution: state_8th - dt*err_factor*(f0 + f10 - f11 - f12)
        state_7th = state_8th - dt * err_factor * (f0 + f10 - f11 - f12)

        error_est = abs(state_8th - state_7th)
    end subroutine rkf78_crtbp_step

    subroutine adaptive_integrate_crtbp_rkf78(state, t_start, t_end, &
            rel_tol, abs_tol, dt_min, dt_max, max_steps, &
            times, states, n_steps)
        real(DP), dimension(6), intent(in) :: state
        real(DP), intent(in) :: t_start, t_end
        real(DP), intent(in) :: rel_tol, abs_tol, dt_min, dt_max
        integer,  intent(in) :: max_steps
        real(DP), allocatable, dimension(:), intent(out) :: times
        real(DP), allocatable, dimension(:,:), intent(out) :: states
        integer, intent(out) :: n_steps

        real(DP), dimension(6) :: current_state, next_state_8th, next_state_7th
        real(DP), dimension(6) :: error_vector, scale_vector
        real(DP) :: current_time, dt, wrms_error
        real(DP), parameter :: safety_factor = 0.9_DP
        real(DP), parameter :: exp_power = 0.125_DP  ! 1/8 for RKF78
        integer :: i

        if (allocated(times)) deallocate(times)
        if (allocated(states)) deallocate(states)
        allocate(times(max_steps))
        allocate(states(max_steps, 6))

        current_state = state
        current_time = t_start
        times(1) = current_time
        states(1, :) = current_state
        n_steps = 1
        dt = min(dt_max, (t_end - t_start) / 100.0_DP)

        do i = 2, max_steps
            if (current_time >= t_end) exit
            if (current_time + dt > t_end) dt = t_end - current_time

            call rkf78_crtbp_step(current_state, dt, current_time, &
                next_state_8th, next_state_7th, error_vector)

            scale_vector = abs_tol + rel_tol * abs(current_state)
            wrms_error = sqrt(sum((error_vector / scale_vector)**2) / 6.0_DP)

            if (ieee_is_nan(wrms_error)) then
                print *, 'FATAL: WRMS error is NaN, aborting CRTBP RKF78 integration.'
                exit
            end if

            if (wrms_error <= 1.0_DP) then
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
        end do

        if (n_steps < max_steps) then
            block
                real(DP), allocatable, dimension(:) :: tmp_times
                real(DP), allocatable, dimension(:,:) :: tmp_states
                allocate(tmp_times(n_steps))
                allocate(tmp_states(n_steps, 6))
                tmp_times = times(1:n_steps)
                tmp_states = states(1:n_steps, :)
                call move_alloc(tmp_times, times)
                call move_alloc(tmp_states, states)
            end block
        end if
    end subroutine adaptive_integrate_crtbp_rkf78

    ! ================================================================
    !  DA RKF78
    ! ================================================================

    subroutine da_rkf78_crtbp_step(state, dt, t, state_8th, state_7th, error_est, pool)
        type(AlgebraicVector), intent(in)    :: state
        real(DP), intent(in) :: dt, t
        type(AlgebraicVector), intent(inout) :: state_8th, state_7th
        type(AlgebraicVector), intent(inout) :: error_est
        type(CrtbpIntegratorTempPool), intent(inout) :: pool

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

        ! f0
        call crtbp_derivatives_da(state%elements, f0%elements, t, pool%force_pool)

        ! f1: state + dt*b10*f0, t + a1*dt
        call real_mul_vector_sub(dt*b10, f0, pool%tmp1)
        call vec_add(state, pool%tmp1, pool%temp_state)
        call crtbp_derivatives_da(pool%temp_state%elements, f1%elements, t + a1*dt, pool%force_pool)

        ! f2: state + dt*(b20*f0 + b21*f1), t + a2*dt
        call real_mul_vector_sub(dt*b20, f0, pool%tmp1)
        call real_mul_vector_sub(dt*b21, f1, pool%tmp2)
        call vec_add(pool%tmp1, pool%tmp2, pool%tmp3)
        call vec_add(state, pool%tmp3, pool%temp_state)
        call crtbp_derivatives_da(pool%temp_state%elements, f2%elements, t + a2*dt, pool%force_pool)

        ! f3: state + dt*(b30*f0 + b32*f2), t + a3*dt
        call real_mul_vector_sub(dt*b30, f0, pool%tmp1)
        call real_mul_vector_sub(dt*b32, f2, pool%tmp2)
        call vec_add(pool%tmp1, pool%tmp2, pool%tmp3)
        call vec_add(state, pool%tmp3, pool%temp_state)
        call crtbp_derivatives_da(pool%temp_state%elements, f3%elements, t + a3*dt, pool%force_pool)

        ! f4: state + dt*(b40*f0 + b42*f2 + b43*f3), t + a4*dt
        call real_mul_vector_sub(dt*b40, f0, pool%tmp1)
        call real_mul_vector_sub(dt*b42, f2, pool%tmp2)
        call real_mul_vector_sub(dt*b43, f3, pool%tmp3)
        call vec_add(pool%tmp1, pool%tmp2, pool%t4)
        call vec_add(pool%t4, pool%tmp3, pool%t5)
        call vec_add(state, pool%t5, pool%temp_state)
        call crtbp_derivatives_da(pool%temp_state%elements, f4%elements, t + a4*dt, pool%force_pool)

        ! f5: state + dt*(b50*f0 + b53*f3 + b54*f4), t + a5*dt
        call real_mul_vector_sub(dt*b50, f0, pool%tmp1)
        call real_mul_vector_sub(dt*b53, f3, pool%tmp2)
        call real_mul_vector_sub(dt*b54, f4, pool%tmp3)
        call vec_add(pool%tmp1, pool%tmp2, pool%t4)
        call vec_add(pool%t4, pool%tmp3, pool%t5)
        call vec_add(state, pool%t5, pool%temp_state)
        call crtbp_derivatives_da(pool%temp_state%elements, f5%elements, t + a5*dt, pool%force_pool)

        ! f6: state + dt*(b60*f0 + b63*f3 + b64*f4 + b65*f5), t + a6*dt
        call real_mul_vector_sub(dt*b60, f0, pool%tmp1)
        call real_mul_vector_sub(dt*b63, f3, pool%tmp2)
        call real_mul_vector_sub(dt*b64, f4, pool%tmp3)
        call real_mul_vector_sub(dt*b65, f5, pool%t4)
        call vec_add(pool%tmp1, pool%tmp2, pool%t5)
        call vec_add(pool%t5, pool%tmp3, pool%t6)
        call vec_add(pool%t6, pool%t4, pool%t7)
        call vec_add(state, pool%t7, pool%temp_state)
        call crtbp_derivatives_da(pool%temp_state%elements, f6%elements, t + a6*dt, pool%force_pool)

        ! f7: state + dt*(b70*f0 + b74*f4 + b75*f5 + b76*f6), t + a7*dt
        call real_mul_vector_sub(dt*b70, f0, pool%tmp1)
        call real_mul_vector_sub(dt*b74, f4, pool%tmp2)
        call real_mul_vector_sub(dt*b75, f5, pool%tmp3)
        call real_mul_vector_sub(dt*b76, f6, pool%t4)
        call vec_add(pool%tmp1, pool%tmp2, pool%t5)
        call vec_add(pool%t5, pool%tmp3, pool%t6)
        call vec_add(pool%t6, pool%t4, pool%t7)
        call vec_add(state, pool%t7, pool%temp_state)
        call crtbp_derivatives_da(pool%temp_state%elements, f7%elements, t + a7*dt, pool%force_pool)

        ! f8: state + dt*(b80*f0 + b83*f3 + b84*f4 + b85*f5 + b86*f6 + b87*f7), t + a8*dt
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
        call crtbp_derivatives_da(pool%temp_state%elements, f8%elements, t + a8*dt, pool%force_pool)

        ! f9: state + dt*(b90*f0 + b93*f3 + b94*f4 + b95*f5 + b96*f6 + b97*f7 + b98*f8), t + a9*dt
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
        call crtbp_derivatives_da(pool%temp_state%elements, f9%elements, t + a9*dt, pool%force_pool)

        ! f10: state + dt*(b100*f0 + b103*f3 + b104*f4 + b105*f5 + b106*f6 + b107*f7 + b108*f8 + b109*f9), t + dt
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
        call crtbp_derivatives_da(pool%temp_state%elements, f10%elements, t + dt, pool%force_pool)

        ! f11: state + dt*(b110*f0 + b115*f5 + b116*f6 + b117*f7 + b118*f8 + b119*f9), t
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
        call crtbp_derivatives_da(pool%temp_state%elements, f11%elements, t, pool%force_pool)

        ! f12: state + dt*(b120*f0 + b123*f3 + b124*f4 + b125*f5 + b126*f6 + b127*f7 + b128*f8 + b129*f9 + f11), t + dt
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
        call crtbp_derivatives_da(pool%temp_state%elements, f12%elements, t + dt, pool%force_pool)

        ! 8th-order: state + dt*(f5*c5 + f6*c6 + f7*c7 + f8*c8 + f9*c9 + f11*c11 + f12*c12)
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

        ! 7th-order: state_8th - dt*err_factor*(f0 + f10 - f11 - f12)
        call vec_add(f0, f10, pool%tmp1)
        call vec_sub(pool%tmp1, f11, pool%tmp2)
        call vec_sub(pool%tmp2, f12, pool%tmp3)
        call real_mul_vector_sub(dt * err_factor, pool%tmp3, pool%t4)
        call vec_sub(state_8th, pool%t4, state_7th)

        ! error_est = state_8th - state_7th
        call vec_sub(state_8th, state_7th, error_est)

        999 continue
        call f0%destroy(); call f1%destroy(); call f2%destroy(); call f3%destroy()
        call f4%destroy(); call f5%destroy(); call f6%destroy(); call f7%destroy()
        call f8%destroy(); call f9%destroy(); call f10%destroy(); call f11%destroy()
        call f12%destroy()
    end subroutine da_rkf78_crtbp_step

    subroutine da_adaptive_integrate_crtbp_rkf78(state, t_start, t_end, &
            rel_tol, abs_tol, dt_min, dt_max, max_steps, &
            final_state)
        type(AlgebraicVector), intent(in)    :: state
        real(DP), intent(in) :: t_start, t_end
        real(DP), intent(in) :: rel_tol, abs_tol, dt_min, dt_max
        integer,  intent(in) :: max_steps
        type(AlgebraicVector), intent(inout) :: final_state

        type(CrtbpIntegratorTempPool) :: pool
        type(AlgebraicVector) :: current_state, next_8th, next_7th, error_vec
        real(DP) :: current_time, dt, wrms_error
        real(DP), dimension(6) :: scale_vector, cons_current
        integer :: i, j

        real(DP), parameter :: safety_factor = 0.9_DP
        real(DP), parameter :: exp_power = 0.125_DP  ! 1/8 for RKF78

        call pool%init()
        call current_state%init(6)
        call next_8th%init(6)
        call next_7th%init(6)
        call error_vec%init(6)

        do j = 1, 6
            current_state%elements(j) = state%elements(j)
        end do
        current_time = t_start
        dt = min(dt_max, (t_end - t_start) / 100.0_DP)

        do i = 2, max_steps
            if (current_time >= t_end) exit
            if (current_time + dt > t_end) dt = t_end - current_time

            call da_rkf78_crtbp_step(current_state, dt, current_time, &
                next_8th, next_7th, error_vec, pool)

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
                print *, 'FATAL: DA WRMS error is NaN, aborting CRTBP RKF78 integration.'
                exit
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
        call pool%destroy()
    end subroutine da_adaptive_integrate_crtbp_rkf78

end module pod_crtbp_integrator_module
