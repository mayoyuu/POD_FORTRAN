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
    public :: CrtbpIntegratorTempPool

    type, public :: CrtbpIntegratorTempPool
        type(AlgebraicVector) :: k1, k2, k3, k4, k5, k6, k7
        type(AlgebraicVector) :: tmp1, tmp2, tmp3
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
        call this%temp_state%init(6)
        call this%force_pool%init()
    end subroutine crtbp_int_pool_init

    subroutine crtbp_int_pool_destroy(this)
        class(CrtbpIntegratorTempPool), intent(inout) :: this
        call this%k1%destroy(); call this%k2%destroy(); call this%k3%destroy()
        call this%k4%destroy(); call this%k5%destroy(); call this%k6%destroy()
        call this%k7%destroy()
        call this%tmp1%destroy(); call this%tmp2%destroy(); call this%tmp3%destroy()
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

end module pod_crtbp_integrator_module
