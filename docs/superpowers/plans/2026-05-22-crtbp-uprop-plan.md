# CRTBP Uncertainty Propagation — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax.

**Goal:** Add CRTBP force model + RKF45 integrator + MC/DA uncertainty propagation + standalone executable with JSON config, all in independent files.

**Architecture:** Five new files, zero modifications to existing code. CRTBP module provides real and DA EOM. Integrator module provides standalone RKF45 (real + DA) with adaptive step control. MC propagator parallel-integrates particle samples. DA propagator does one DA integration → polynomial compile → parallel sample evaluation. Executable reads JSON config and routes to MC or DA.

**Tech Stack:** Fortran 2008+, `pod_global` (DP), `pod_dace_classes` (DA/AlgebraicVector/CompiledDA), `pod_random_module` (multivariate normal), OpenMP.

---

### Task 1: Create `pod_crtbp_module.f90` — CRTBP Equations of Motion

**Files:**
- Create: `POD_Fortran/src/lib/forcemodel/pod_crtbp_module.f90`

- [ ] **Step 1: Write the complete module**

```fortran
module pod_crtbp_module
    use pod_global, only: DP
    use pod_dace_classes
    implicit none
    private

    real(DP), public :: crtbp_mu = 0.01_DP
    public :: set_crtbp_mu
    public :: crtbp_derivatives_real
    public :: crtbp_derivatives_da
    public :: CrtbpForceTempPool

    type, public :: CrtbpForceTempPool
        type(DA) :: r1, r2, v1, v2
        type(DA) :: Omega_x, Omega_y, Omega_z
        type(DA) :: tmp1, tmp2, tmp3
    contains
        procedure :: init   => crtbp_force_pool_init
        procedure :: destroy => crtbp_force_pool_destroy
    end type CrtbpForceTempPool

contains

    subroutine set_crtbp_mu(mu)
        real(DP), intent(in) :: mu
        crtbp_mu = mu
    end subroutine set_crtbp_mu

    subroutine crtbp_force_pool_init(this)
        class(CrtbpForceTempPool), intent(inout) :: this
        call this%r1%init();      call this%r2%init()
        call this%v1%init();      call this%v2%init()
        call this%Omega_x%init(); call this%Omega_y%init(); call this%Omega_z%init()
        call this%tmp1%init();    call this%tmp2%init();    call this%tmp3%init()
    end subroutine crtbp_force_pool_init

    subroutine crtbp_force_pool_destroy(this)
        class(CrtbpForceTempPool), intent(inout) :: this
        call this%r1%destroy();      call this%r2%destroy()
        call this%v1%destroy();      call this%v2%destroy()
        call this%Omega_x%destroy(); call this%Omega_y%destroy(); call this%Omega_z%destroy()
        call this%tmp1%destroy();    call this%tmp2%destroy();    call this%tmp3%destroy()
    end subroutine crtbp_force_pool_destroy

    subroutine crtbp_derivatives_real(x, dxdt, t)
        real(DP), dimension(6), intent(in)  :: x
        real(DP), dimension(6), intent(out) :: dxdt
        real(DP), intent(in) :: t

        real(DP) :: r1, r2, v1, v2
        real(DP) :: Omega_x, Omega_y, Omega_z
        real(DP) :: mu

        mu = crtbp_mu

        r1 = sqrt((x(1) + mu)**2 + x(2)**2 + x(3)**2)
        r2 = sqrt((x(1) - 1.0_DP + mu)**2 + x(2)**2 + x(3)**2)
        v1 = (1.0_DP - mu) / (r1**3)
        v2 = mu / (r2**3)

        Omega_x = x(1) - (x(1) + mu) * v1 - v2 * (x(1) - 1.0_DP + mu)
        Omega_y = x(2) * (1.0_DP - v1 - v2)
        Omega_z = -x(3) * (v1 + v2)

        dxdt(1) = x(4)
        dxdt(2) = x(5)
        dxdt(3) = x(6)
        dxdt(4) = Omega_x + 2.0_DP * x(5)
        dxdt(5) = Omega_y - 2.0_DP * x(4)
        dxdt(6) = Omega_z
    end subroutine crtbp_derivatives_real

    subroutine crtbp_derivatives_da(x, dxdt, t, pool)
        type(DA), dimension(6), intent(in)    :: x
        type(DA), dimension(6), intent(inout) :: dxdt
        real(DP), intent(in) :: t
        type(CrtbpForceTempPool), intent(inout) :: pool

        real(DP) :: mu

        mu = crtbp_mu

        ! r1 = sqrt((x1+mu)^2 + x2^2 + x3^2)
        pool%tmp1 = x(1) + mu
        pool%r1 = pool%tmp1 * pool%tmp1 + x(2) * x(2) + x(3) * x(3)
        call da_sqrt_sub(pool%r1, pool%r1)

        ! r2 = sqrt((x1-1+mu)^2 + x2^2 + x3^2)
        pool%tmp2 = x(1) - 1.0_DP + mu
        pool%r2 = pool%tmp2 * pool%tmp2 + x(2) * x(2) + x(3) * x(3)
        call da_sqrt_sub(pool%r2, pool%r2)

        ! v1 = (1-mu) / r1^3
        pool%tmp1 = pool%r1 * pool%r1 * pool%r1
        pool%v1 = (1.0_DP - mu) / pool%tmp1

        ! v2 = mu / r2^3
        pool%tmp1 = pool%r2 * pool%r2 * pool%r2
        pool%v2 = mu / pool%tmp1

        ! Omega_x = x1 - (x1+mu)*v1 - v2*(x1-1+mu)
        pool%Omega_x = x(1) - (x(1) + mu) * pool%v1 - pool%v2 * (x(1) - 1.0_DP + mu)

        ! Omega_y = x2 * (1 - v1 - v2)
        pool%tmp1 = 1.0_DP - pool%v1 - pool%v2
        pool%Omega_y = x(2) * pool%tmp1

        ! Omega_z = -x3 * (v1 + v2)
        pool%tmp1 = pool%v1 + pool%v2
        pool%Omega_z = -x(3) * pool%tmp1

        dxdt(1) = x(4)
        dxdt(2) = x(5)
        dxdt(3) = x(6)
        dxdt(4) = pool%Omega_x + 2.0_DP * x(5)
        dxdt(5) = pool%Omega_y - 2.0_DP * x(4)
        dxdt(6) = pool%Omega_z
    end subroutine crtbp_derivatives_da

end module pod_crtbp_module
```

- [ ] **Step 2: Verify file compiles (syntax check)**

Run: `cd POD_Fortran && fpm build 2>&1 | head -20`
Expected: No errors from this file. Other files may still reference missing symbols — that's OK at this stage.

- [ ] **Step 3: Commit**

```bash
git add POD_Fortran/src/lib/forcemodel/pod_crtbp_module.f90
git commit -m "feat: add CRTBP equations of motion module (real + DA)"
```

---

### Task 2: Create `pod_crtbp_integrator_module.f90` — RKF45 Integrator

**Files:**
- Create: `POD_Fortran/src/lib/integrator/pod_crtbp_integrator_module.f90`

- [ ] **Step 1: Write the complete module (real RKF45 step + adaptive, DA RKF45 step + adaptive + temp pool)**

```fortran
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
        call crtbp_derivatives_real(temp_state, t + a2 * dt, k2)

        temp_state = state + dt * (b31 * k1 + b32 * k2)
        call crtbp_derivatives_real(temp_state, t + a3 * dt, k3)

        temp_state = state + dt * (b41 * k1 + b42 * k2 + b43 * k3)
        call crtbp_derivatives_real(temp_state, t + a4 * dt, k4)

        temp_state = state + dt * (b51 * k1 + b52 * k2 + b53 * k3 + b54 * k4)
        call crtbp_derivatives_real(temp_state, t + a5 * dt, k5)

        temp_state = state + dt * (b61 * k1 + b62 * k2 + b63 * k3 + b64 * k4 + b65 * k5)
        call crtbp_derivatives_real(temp_state, t + dt, k6)

        state_5th = state + dt * (c1 * k1 + c3 * k3 + c4 * k4 + c5 * k5 + c6 * k6)
        call crtbp_derivatives_real(state_5th, t + dt, k7)

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
            d_state(j) = state%elements(j)
            d_k1(j)%init(); d_k2(j)%init(); d_k3(j)%init()
            d_k4(j)%init(); d_k5(j)%init(); d_k6(j)%init(); d_k7(j)%init()
            d_temp(j)%init()
            d_5th(j)%init(); d_4th(j)%init(); d_err(j)%init()
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

        current_state = state
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
```

- [ ] **Step 2: Verify compilation**

Run: `cd POD_Fortran && fpm build 2>&1 | head -30`
Expected: No compile errors.

- [ ] **Step 3: Commit**

```bash
git add POD_Fortran/src/lib/integrator/pod_crtbp_integrator_module.f90
git commit -m "feat: add CRTBP RKF45 integrator (real + DA) with adaptive step control"
```

---

### Task 3: Create `pod_uq_crtbp_mc_module.f90` — MC Propagator

**Files:**
- Create: `POD_Fortran/src/lib/uncertainty/propagation/pod_uq_crtbp_mc_module.f90`

- [ ] **Step 1: Write the complete module**

```fortran
module pod_uq_crtbp_mc_module
    use pod_global, only: DP, MAX_STRING_LEN
    use pod_crtbp_module, only: set_crtbp_mu
    use pod_crtbp_integrator_module, only: adaptive_integrate_crtbp
    use pod_random_module, only: generate_multivariate_normal, init_random_seed
    implicit none
    private

    public :: crtbp_mc_propagate

contains

    subroutine crtbp_mc_propagate(nominal_state, cov, mu, t_end, &
            n_particles, rel_tol, abs_tol, dt_min, dt_max, max_steps, &
            final_samples, final_mean, final_cov, verbose)
        real(DP), intent(in) :: nominal_state(6)
        real(DP), intent(in) :: cov(6,6)
        real(DP), intent(in) :: mu, t_end
        integer,  intent(in) :: n_particles
        real(DP), intent(in) :: rel_tol, abs_tol, dt_min, dt_max
        integer,  intent(in) :: max_steps
        real(DP), allocatable, intent(out) :: final_samples(:,:)
        real(DP), intent(out) :: final_mean(6)
        real(DP), intent(out) :: final_cov(6,6)
        logical, intent(in) :: verbose

        integer :: dim, i, n_steps, j
        real(DP), allocatable :: temp_times(:), temp_states(:,:)
        real(DP), dimension(6) :: current_state

        dim = 6

        allocate(final_samples(dim, n_particles))

        call set_crtbp_mu(mu)
        call init_random_seed(.true.)
        call generate_multivariate_normal(nominal_state, cov, final_samples)

        if (verbose) write(*, '(A,I0)') '[MC CRTBP] Propagating ', n_particles

        !$omp parallel do default(none) &
        !$omp private(i, current_state, temp_times, temp_states, n_steps) &
        !$omp shared(n_particles, final_samples, t_end, rel_tol, abs_tol, dt_min, dt_max, max_steps)
        do i = 1, n_particles
            current_state = final_samples(:, i)

            call adaptive_integrate_crtbp(current_state, 0.0_DP, t_end, &
                rel_tol, abs_tol, dt_min, dt_max, max_steps, &
                temp_times, temp_states, n_steps)

            final_samples(:, i) = temp_states(n_steps, :)

            if (allocated(temp_times)) deallocate(temp_times)
            if (allocated(temp_states)) deallocate(temp_states)
        end do
        !$omp end parallel do

        ! Compute mean
        final_mean = 0.0_DP
        do i = 1, n_particles
            final_mean = final_mean + final_samples(:, i)
        end do
        final_mean = final_mean / real(n_particles, DP)

        ! Compute covariance
        final_cov = 0.0_DP
        do i = 1, n_particles
            do j = 1, dim
                final_cov(:, j) = final_cov(:, j) + &
                    (final_samples(:, i) - final_mean) * (final_samples(j, i) - final_mean(j))
            end do
        end do
        final_cov = final_cov / real(n_particles - 1, DP)

        if (verbose) write(*, '(A)') '[MC CRTBP] Propagation complete.'
    end subroutine crtbp_mc_propagate

end module pod_uq_crtbp_mc_module
```

- [ ] **Step 2: Verify compilation**

Run: `cd POD_Fortran && fpm build 2>&1 | head -20`
Expected: No errors.

- [ ] **Step 3: Commit**

```bash
git add POD_Fortran/src/lib/uncertainty/propagation/pod_uq_crtbp_mc_module.f90
git commit -m "feat: add CRTBP Monte Carlo uncertainty propagator"
```

---

### Task 4: Create `pod_uq_crtbp_da_module.f90` — DA Propagator

**Files:**
- Create: `POD_Fortran/src/lib/uncertainty/propagation/pod_uq_crtbp_da_module.f90`

- [ ] **Step 1: Write the complete module**

```fortran
module pod_uq_crtbp_da_module
    use pod_global, only: DP
    use pod_crtbp_module, only: set_crtbp_mu
    use pod_crtbp_integrator_module, only: da_adaptive_integrate_crtbp
    use pod_dace_classes
    use pod_random_module, only: generate_multivariate_normal, init_random_seed
    implicit none
    private

    public :: crtbp_da_propagate

contains

    subroutine crtbp_da_propagate(nominal_state, cov, mu, t_end, &
            da_order, n_particles, rel_tol, abs_tol, dt_min, dt_max, max_steps, &
            final_samples, final_mean, final_cov, propagated_ref, verbose)
        real(DP), intent(in) :: nominal_state(6)
        real(DP), intent(in) :: cov(6,6)
        real(DP), intent(in) :: mu, t_end
        integer,  intent(in) :: da_order, n_particles
        real(DP), intent(in) :: rel_tol, abs_tol, dt_min, dt_max
        integer,  intent(in) :: max_steps
        real(DP), allocatable, intent(out) :: final_samples(:,:)
        real(DP), intent(out) :: final_mean(6)
        real(DP), intent(out) :: final_cov(6,6)
        real(DP), intent(out) :: propagated_ref(6)
        logical, intent(in) :: verbose

        type(AlgebraicVector) :: state_da_0, state_da_f
        type(CompiledDA) :: compiled
        integer :: i, j, dim
        real(DP) :: eval_inputs(6), eval_results(6)
        real(DP), allocatable :: particles(:,:)

        dim = 6

        ! 1. Generate particles (also used as deviation samples for DA eval)
        allocate(particles(dim, n_particles))
        call init_random_seed(.true.)
        call generate_multivariate_normal(nominal_state, cov, particles)

        ! 2. Set mu and push DA order
        call set_crtbp_mu(mu)
        call dace_push_to(da_order)

        ! 3. Build initial DA state: nominal + da_var(i)
        call state_da_0%init(dim)
        do i = 1, dim
            state_da_0%elements(i) = nominal_state(i) + da_var(i)
        end do

        ! 4. Single DA adaptive integration (returns final state only)
        call state_da_f%init(dim)
        if (verbose) write(*, '(A)') '[DA CRTBP] Starting DA integration...'

        call da_adaptive_integrate_crtbp(state_da_0, 0.0_DP, t_end, &
            rel_tol, abs_tol, dt_min, dt_max, max_steps, &
            state_da_f)

        ! 5. Extract constant part (propagated reference)
        propagated_ref = state_da_f%cons()

        ! 6. Compile polynomial
        compiled = state_da_f%compile()

        ! 7. Evaluate all particles
        allocate(final_samples(dim, n_particles))
        eval_inputs = 0.0_DP

        !$omp parallel do default(none) &
        !$omp private(i, eval_inputs, eval_results) &
        !$omp shared(n_particles, particles, nominal_state, compiled, final_samples, dim)
        do i = 1, n_particles
            eval_inputs(:) = particles(:, i) - nominal_state(:)
            eval_results = compiled%eval(eval_inputs)
            final_samples(:, i) = eval_results
        end do
        !$omp end parallel do

        ! 8. Compute moments
        final_mean = 0.0_DP
        do i = 1, n_particles
            final_mean = final_mean + final_samples(:, i)
        end do
        final_mean = final_mean / real(n_particles, DP)

        final_cov = 0.0_DP
        do i = 1, n_particles
            do j = 1, dim
                final_cov(:, j) = final_cov(:, j) + &
                    (final_samples(:, i) - final_mean) * (final_samples(j, i) - final_mean(j))
            end do
        end do
        final_cov = final_cov / real(n_particles - 1, DP)

        ! 9. Cleanup
        call compiled%destroy()
        call state_da_0%destroy()
        call state_da_f%destroy()
        deallocate(particles)
        call dace_pop_to()

        if (verbose) write(*, '(A)') '[DA CRTBP] Propagation complete.'
    end subroutine crtbp_da_propagate

end module pod_uq_crtbp_da_module
```

- [ ] **Step 2: Verify compilation**

Run: `cd POD_Fortran && fpm build 2>&1 | head -20`
Expected: No errors.

- [ ] **Step 3: Commit**

```bash
git add POD_Fortran/src/lib/uncertainty/propagation/pod_uq_crtbp_da_module.f90
git commit -m "feat: add CRTBP Differential Algebra uncertainty propagator"
```

---

### Task 5: Create example JSON config

**Files:**
- Create: `POD_Fortran/config/crtbp_uprop_example.json`

- [ ] **Step 1: Write the example JSON**

```json
{
    "mu": 0.01,
    "initial_state": [0.5, 0.0, 0.0, 0.0, 1.0, 0.0],
    "initial_covariance": [
        [1.0e-8, 0.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 1.0e-8, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 1.0e-8, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 1.0e-8, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, 1.0e-8, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0, 1.0e-8]
    ],
    "propagation_time": 10.0,
    "method": "MC",
    "num_samples": 10000,
    "da_order": 2,
    "integrator": {
        "rel_tol": 1.0e-12,
        "abs_tol": 1.0e-12,
        "min_step": 1.0e-10,
        "max_step": 0.1,
        "max_steps": 100000
    },
    "output": {
        "prefix": "./output/crtbp_uprop"
    }
}
```

- [ ] **Step 2: Commit**

```bash
git add POD_Fortran/config/crtbp_uprop_example.json
git commit -m "feat: add example JSON config for CRTBP uncertainty propagation"
```

---

### Task 6: Create `run_CRTBP_uprop.f90` — Executable Program

**Files:**
- Create: `POD_Fortran/app/run_CRTBP_uprop.f90`

- [ ] **Step 1: Write the complete executable with embedded JSON parser**

```fortran
program run_CRTBP_uprop
    use pod_global, only: DP, MAX_STRING_LEN
    use pod_uq_crtbp_mc_module, only: crtbp_mc_propagate
    use pod_uq_crtbp_da_module, only: crtbp_da_propagate
    implicit none

    character(len=MAX_STRING_LEN) :: json_path
    character(len=:), allocatable :: content
    character(len=:), allocatable :: obj_content
    logical :: ok

    real(DP) :: mu, t_end
    real(DP) :: nominal_state(6)
    real(DP) :: cov(6,6)
    character(len=10) :: method
    integer  :: n_samples, da_order
    real(DP) :: rel_tol, abs_tol, dt_min, dt_max
    integer  :: max_steps_int
    character(len=MAX_STRING_LEN) :: output_prefix

    real(DP), allocatable :: final_samples(:,:)
    real(DP) :: final_mean(6), final_cov(6,6), propagated_ref(6)

    integer :: i, j, file_unit

    ! ---- CLI parsing ----
    if (command_argument_count() < 2) then
        print *, 'Usage: run_CRTBP_uprop --in <config.json>'
        stop
    end if

    call get_command_argument(1, json_path)
    if (trim(json_path) /= '--in') then
        print *, 'Error: expected --in <config.json>'
        stop
    end if
    call get_command_argument(2, json_path)

    ! ---- Read JSON file ----
    call read_file_content(json_path, content, ok)
    if (.not. ok) then
        print *, 'Error: cannot read file ', trim(json_path)
        stop
    end if

    ! ---- Parse required fields ----
    call json_get_real(content, 'mu', mu, ok)
    if (.not. ok) stop 'Error: missing "mu"'

    call json_get_array(content, 'initial_state', nominal_state, 6, ok)
    if (.not. ok) stop 'Error: missing "initial_state"'

    call json_get_matrix(content, 'initial_covariance', cov, 6, 6, ok)
    if (.not. ok) stop 'Error: missing "initial_covariance"'

    call json_get_real(content, 'propagation_time', t_end, ok)
    if (.not. ok) stop 'Error: missing "propagation_time"'

    call json_get_string(content, 'method', method, ok)
    if (.not. ok) stop 'Error: missing "method"'

    call json_get_int(content, 'num_samples', n_samples, ok)
    if (.not. ok) n_samples = 10000

    call json_get_int(content, 'da_order', da_order, ok)
    if (.not. ok) da_order = 2

    ! ---- Parse nested integrator object ----
    call json_get_object(content, 'integrator', obj_content, ok)
    if (ok) then
        call json_get_real(obj_content, 'rel_tol', rel_tol, ok)
        if (.not. ok) rel_tol = 1.0e-12_DP
        call json_get_real(obj_content, 'abs_tol', abs_tol, ok)
        if (.not. ok) abs_tol = 1.0e-12_DP
        call json_get_real(obj_content, 'min_step', dt_min, ok)
        if (.not. ok) dt_min = 1.0e-10_DP
        call json_get_real(obj_content, 'max_step', dt_max, ok)
        if (.not. ok) dt_max = 0.1_DP
        call json_get_int(obj_content, 'max_steps', max_steps_int, ok)
        if (.not. ok) max_steps_int = 100000
        deallocate(obj_content)
    else
        rel_tol = 1.0e-12_DP
        abs_tol = 1.0e-12_DP
        dt_min  = 1.0e-10_DP
        dt_max  = 0.1_DP
        max_steps_int = 100000
    end if

    ! ---- Parse output prefix ----
    call json_get_object(content, 'output', obj_content, ok)
    if (ok) then
        call json_get_string(obj_content, 'prefix', output_prefix, ok)
        if (.not. ok) output_prefix = './output/crtbp_uprop'
        deallocate(obj_content)
    else
        output_prefix = './output/crtbp_uprop'
    end if

    deallocate(content)

    ! ---- Validate non-dimensional input ----
    if (maxval(abs(nominal_state(1:3))) > 100.0_DP .or. &
        maxval(abs(nominal_state(4:6))) > 10.0_DP) then
        print *, 'ERROR: Input state appears to be dimensional.'
        print *, 'CRTBP expects non-dimensional units (positions ~O(1), velocities ~O(1)).'
        print *, 'Position magnitudes found: ', maxval(abs(nominal_state(1:3)))
        print *, 'Velocity magnitudes found: ', maxval(abs(nominal_state(4:6)))
        print *, 'If input is truly non-dimensional, adjust the threshold in source code.'
        stop
    end if

    ! ---- Print config summary ----
    print *, '========================================'
    print *, '  CRTBP Uncertainty Propagation'
    print *, '========================================'
    print '(A,F8.6)',  '  mu              = ', mu
    print '(A,6F10.6)', '  initial state   = ', nominal_state
    print '(A,F12.4)',  '  propagation T   = ', t_end
    print '(A,A)',      '  method          = ', trim(method)
    if (trim(method) == 'MC') then
        print '(A,I0)', '  num_samples     = ', n_samples
    else
        print '(A,I0)', '  da_order        = ', da_order
        print '(A,I0)', '  num_samples     = ', n_samples
    end if
    print '(A,ES10.3)', '  rel_tol         = ', rel_tol
    print '(A,ES10.3)', '  abs_tol         = ', abs_tol
    print '(A,A)',      '  output prefix   = ', trim(output_prefix)
    print *, '========================================'

    ! ---- Route to MC or DA ----
    if (trim(method) == 'MC') then
        call crtbp_mc_propagate(nominal_state, cov, mu, t_end, &
            n_samples, rel_tol, abs_tol, dt_min, dt_max, max_steps_int, &
            final_samples, final_mean, final_cov, .true.)
    else if (trim(method) == 'DA') then
        call crtbp_da_propagate(nominal_state, cov, mu, t_end, &
            da_order, n_samples, rel_tol, abs_tol, dt_min, dt_max, max_steps_int, &
            final_samples, final_mean, final_cov, propagated_ref, .true.)
    else
        print *, 'Error: method must be "MC" or "DA", got: ', trim(method)
        stop
    end if

    ! ---- Write output CSV files ----
    call write_csv_particles(trim(output_prefix) // '_particles.csv', final_samples)
    call write_csv_stats(trim(output_prefix) // '_stats.csv', final_mean, final_cov)

    print *, 'Results saved to:'
    print *, '  ', trim(output_prefix) // '_particles.csv'
    print *, '  ', trim(output_prefix) // '_stats.csv'

    deallocate(final_samples)

contains

    ! ================================================================
    !  JSON PARSER (embedded, minimal)
    ! ================================================================

    subroutine read_file_content(filename, content, ok)
        character(len=*), intent(in) :: filename
        character(len=:), allocatable, intent(out) :: content
        logical, intent(out) :: ok
        character(len=4096) :: line
        integer :: unit, ios, total_len
        character(len=:), allocatable :: temp

        ok = .false.
        open(newunit=unit, file=trim(filename), status='old', action='read', iostat=ios)
        if (ios /= 0) return

        allocate(character(len=0) :: content)
        total_len = 0
        do
            read(unit, '(A)', iostat=ios) line
            if (ios /= 0) exit
            temp = content // trim(adjustl(line))
            deallocate(content)
            content = temp
        end do
        close(unit)
        total_len = len(content)
        if (total_len > 0) ok = .true.
    end subroutine read_file_content

    function find_key_pos(content, key) result(pos)
        character(len=*), intent(in) :: content, key
        integer :: pos
        character(len=:), allocatable :: search
        search = '"' // trim(key) // '"'
        pos = index(content, search)
    end function find_key_pos

    subroutine json_get_real(content, key, value, found)
        character(len=*), intent(in) :: content, key
        real(DP), intent(out) :: value
        logical, intent(out) :: found
        integer :: p, p_end, ios
        character(len=256) :: token

        found = .false.
        value = 0.0_DP
        p = find_key_pos(content, key)
        if (p == 0) return
        p = p + len('"' // trim(key) // '"')
        p = scan_forward(content, p, ':')
        if (p == 0) return
        p = p + 1
        p = skip_whitespace(content, p)
        p_end = scan_forward(content, p, ',}]')
        if (p_end == 0) p_end = len(content)
        token = content(p:p_end-1)
        token = trim(adjustl(token))
        read(token, *, iostat=ios) value
        if (ios == 0) found = .true.
    end subroutine json_get_real

    subroutine json_get_int(content, key, value, found)
        character(len=*), intent(in) :: content, key
        integer, intent(out) :: value
        logical, intent(out) :: found
        integer :: p, p_end, ios
        character(len=256) :: token

        found = .false.
        value = 0
        p = find_key_pos(content, key)
        if (p == 0) return
        p = p + len('"' // trim(key) // '"')
        p = scan_forward(content, p, ':')
        if (p == 0) return
        p = p + 1
        p = skip_whitespace(content, p)
        p_end = scan_forward(content, p, ',}]')
        if (p_end == 0) p_end = len(content)
        token = content(p:p_end-1)
        token = trim(adjustl(token))
        read(token, *, iostat=ios) value
        if (ios == 0) found = .true.
    end subroutine json_get_int

    subroutine json_get_string(content, key, value, found)
        character(len=*), intent(in) :: content, key
        character(len=*), intent(out) :: value
        logical, intent(out) :: found
        integer :: p, p_start, p_end

        found = .false.
        value = ''
        p = find_key_pos(content, key)
        if (p == 0) return
        p = p + len('"' // trim(key) // '"')
        p = scan_forward(content, p, ':')
        if (p == 0) return
        p = p + 1
        p = skip_whitespace(content, p)
        if (p > len(content)) return
        if (content(p:p) /= '"') return
        p_start = p + 1
        p_end = scan_forward(content, p_start, '"')
        if (p_end == 0) return
        value = content(p_start:p_end-1)
        found = .true.
    end subroutine json_get_string

    subroutine json_get_array(content, key, arr, n, found)
        character(len=*), intent(in) :: content, key
        real(DP), intent(out) :: arr(n)
        integer, intent(in) :: n
        logical, intent(out) :: found
        integer :: p, p_end, i, ios
        character(len=4096) :: segment

        found = .false.
        arr = 0.0_DP
        p = find_key_pos(content, key)
        if (p == 0) return
        p = p + len('"' // trim(key) // '"')
        p = scan_forward(content, p, '[')
        if (p == 0) return
        p = p + 1
        p_end = find_matching_bracket(content, p, '[', ']')
        if (p_end == 0) return
        segment = content(p:p_end-1)
        do i = 1, n
            p = skip_whitespace(segment, 1)
            p_end = scan_forward(segment, p, ',]')
            if (p_end == 0) p_end = len_trim(segment) + 1
            read(segment(p:p_end-1), *, iostat=ios) arr(i)
            if (ios /= 0) return
            if (i < n .and. p_end <= len_trim(segment)) &
                segment = segment(p_end+1:)
        end do
        found = .true.
    end subroutine json_get_array

    subroutine json_get_matrix(content, key, mat, rows, cols, found)
        character(len=*), intent(in) :: content, key
        real(DP), intent(out) :: mat(rows, cols)
        integer, intent(in) :: rows, cols
        logical, intent(out) :: found
        integer :: p, p_end, i, j, ios, p_seg
        character(len=4096) :: outer_seg, row_seg

        found = .false.
        mat = 0.0_DP
        p = find_key_pos(content, key)
        if (p == 0) return
        p = p + len('"' // trim(key) // '"')
        p = scan_forward(content, p, '[')
        if (p == 0) return
        p = p + 1
        p_end = find_matching_bracket(content, p, '[', ']')
        if (p_end == 0) return
        outer_seg = content(p:p_end-1)

        do i = 1, rows
            p = scan_forward(outer_seg, 1, '[')
            if (p == 0) return
            p = p + 1
            p_end = find_matching_bracket(outer_seg, p, '[', ']')
            if (p_end == 0) return
            row_seg = outer_seg(p:p_end-1)
            do j = 1, cols
                p_seg = skip_whitespace(row_seg, 1)
                p_end = scan_forward(row_seg, p_seg, ',]')
                if (p_end == 0) p_end = len_trim(row_seg) + 1
                read(row_seg(p_seg:p_end-1), *, iostat=ios) mat(i, j)
                if (ios /= 0) return
                if (j < cols .and. p_end <= len_trim(row_seg)) &
                    row_seg = row_seg(p_end+1:)
            end do
            p_end = scan_forward(outer_seg, p, ']')
            if (p_end == 0) return
            outer_seg = outer_seg(p_end+1:)
        end do
        found = .true.
    end subroutine json_get_matrix

    subroutine json_get_object(content, key, obj_content, found)
        character(len=*), intent(in) :: content, key
        character(len=:), allocatable, intent(out) :: obj_content
        logical, intent(out) :: found
        integer :: p, p_end

        found = .false.
        allocate(character(len=0) :: obj_content)
        p = find_key_pos(content, key)
        if (p == 0) return
        p = p + len('"' // trim(key) // '"')
        p = scan_forward(content, p, '{')
        if (p == 0) return
        p_end = find_matching_bracket(content, p, '{', '}')
        if (p_end == 0) return
        obj_content = content(p:p_end)
        found = .true.
    end subroutine json_get_object

    ! ---- Utility functions ----

    function scan_forward(str, start, ch) result(pos)
        character(len=*), intent(in) :: str
        integer, intent(in) :: start
        character, intent(in) :: ch
        integer :: pos
        pos = index(str(start:), ch)
        if (pos > 0) pos = pos + start - 1
    end function scan_forward

    function skip_whitespace(str, start) result(pos)
        character(len=*), intent(in) :: str
        integer, intent(in) :: start
        integer :: pos
        pos = start
        do while (pos <= len(str))
            if (str(pos:pos) /= ' ' .and. str(pos:pos) /= achar(9) .and. &
                str(pos:pos) /= achar(10) .and. str(pos:pos) /= achar(13)) exit
            pos = pos + 1
        end do
    end function skip_whitespace

    function find_matching_bracket(str, open_pos, open_ch, close_ch) result(close_pos)
        character(len=*), intent(in) :: str
        integer, intent(in) :: open_pos
        character, intent(in) :: open_ch, close_ch
        integer :: close_pos
        integer :: depth, i

        depth = 1
        close_pos = 0
        do i = open_pos + 1, len(str)
            if (str(i:i) == open_ch) depth = depth + 1
            if (str(i:i) == close_ch) depth = depth - 1
            if (depth == 0) then
                close_pos = i
                return
            end if
        end do
    end function find_matching_bracket

    ! ---- Output writers ----

    subroutine write_csv_particles(filename, samples)
        character(len=*), intent(in) :: filename
        real(DP), intent(in) :: samples(:,:)
        integer :: unit, i

        open(newunit=unit, file=trim(filename), status='replace', action='write')
        write(unit, '(A)') 'x,y,z,vx,vy,vz'
        do i = 1, size(samples, 2)
            write(unit, '(*(ES22.14, :, ","))') samples(:, i)
        end do
        close(unit)
    end subroutine write_csv_particles

    subroutine write_csv_stats(filename, mean_vec, cov_mat)
        character(len=*), intent(in) :: filename
        real(DP), intent(in) :: mean_vec(6), cov_mat(6,6)
        integer :: unit, i

        open(newunit=unit, file=trim(filename), status='replace', action='write')
        write(unit, '(A)') '# Mean'
        write(unit, '(*(ES22.14, :, ","))') mean_vec(:)
        write(unit, '(A)') '# Covariance Matrix'
        do i = 1, 6
            write(unit, '(*(ES22.14, :, ","))') cov_mat(i, :)
        end do
        close(unit)
    end subroutine write_csv_stats

end program run_CRTBP_uprop
```

- [ ] **Step 2: Verify full project compilation**

Run: `cd POD_Fortran && fpm build 2>&1`
Expected: Full build succeeds. If fpm complains about the new executable not being declared, note that user needs to add to fpm.toml.

- [ ] **Step 3: Commit**

```bash
git add POD_Fortran/app/run_CRTBP_uprop.f90
git commit -m "feat: add CRTBP uncertainty propagation executable with JSON config"
```

---

### Task 7: Update `fpm.toml` — Register new executable

**Files:**
- Modify: `POD_Fortran/fpm.toml`

- [ ] **Step 1: Add executable entry**

Append to `fpm.toml` (after the last `[[executable]]` block):

```toml
[[executable]]
name = "run_CRTBP_uprop"
source-dir = "app"
main = "run_CRTBP_uprop.f90"
```

- [ ] **Step 2: Verify full build with new executable**

Run: `cd POD_Fortran && fpm build 2>&1`
Expected: Full build succeeds including the new executable.

- [ ] **Step 3: Final commit**

```bash
git add POD_Fortran/fpm.toml
git commit -m "build: register run_CRTBP_uprop executable in fpm.toml"
```

---
