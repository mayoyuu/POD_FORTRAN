# DA Order Convergence Test Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a Fortran program that determines the minimum DA order (3-6) for 15-day L1 Halo orbit uncertainty propagation with <1% per-particle relative error in all components vs MC.

**Architecture:** Single test program (`app/test_da_convergence.f90`) that reads L1Halo OPM, generates 1000 shared particles, runs MC once as golden reference, then scans DA orders 3-6 comparing each particle against MC using `max(|DA-MC|/|MC|)` per component. Uses existing `run_particle_propagation` API for both methods.

**Tech Stack:** Fortran + fpm, DACE (Differential Algebra Core Engine), SPICE kernels, RKF45 integrator

---

## File Structure

| File | Action | Responsibility |
|------|--------|---------------|
| `app/test_da_convergence.f90` | Create | Main test program — all 5 phases |
| `fpm.toml` | Modify | Add executable entry |

No library changes needed. All APIs (`load_initial_opm`, `run_particle_propagation`, `uq_state_type`, DACE routines, `generate_multivariate_normal`) already exist.

---

### Task 1: Register executable in fpm.toml

**Files:**
- Modify: `fpm.toml`

- [ ] **Step 1: Add executable entry**

Append to `fpm.toml`:

```toml
[[executable]]
name = "test_da_convergence"
source-dir = "app"
main = "test_da_convergence.f90"
```

- [ ] **Step 2: Commit**

```bash
git add fpm.toml
git commit -m "build: add test_da_convergence executable"
```

---

### Task 2: Write the test program skeleton (Phase 0 — System init)

**Files:**
- Create: `app/test_da_convergence.f90`

- [ ] **Step 1: Write program skeleton with Phase 0 init**

```fortran
program test_da_convergence
    use pod_global, only: DP
    use pod_engine_module, only: pod_engine_init
    use pod_config, only: config
    use pod_spice, only: str2et
    use pod_da_force_model_module, only: init_gravity_network
    use pod_dace_classes, only: dace_initialize
    use pod_data_format_module, only: load_initial_opm

    implicit none

    character(len=*), parameter :: CONFIG_FILE = 'dummy_test_config.txt'
    character(len=*), parameter :: OPM_FILE    = 'OPM/L1Halo/L1Halo_init.opm.json'
    integer, parameter :: DA_MAX_ORDER = 6
    integer, parameter :: DA_NVARS     = 6
    integer, parameter :: N_PARTICLES  = 1000
    real(DP), parameter :: DURATION    = 15.0_DP * 86400.0_DP  ! 15 days in seconds

    write(*,*) '============================================================'
    write(*,*) '  DA Order Convergence Test — L1Halo / 15-day / RKF45'
    write(*,*) '============================================================'
    write(*,*) ''

    ! Phase 0: System init
    write(*,*) '>>> Phase 0: System initialization...'
    call pod_engine_init(CONFIG_FILE)
    config%use_earth_nspheric = .true.
    config%earth_degree       = 10
    config%use_moon_nspheric  = .true.
    config%moon_degree        = 10
    call init_gravity_network()
    call dace_initialize(DA_MAX_ORDER, DA_NVARS)
    write(*,*) '  DACE initialized with max_order=', DA_MAX_ORDER, ', nvars=', DA_NVARS
    write(*,*) ''

    stop
end program test_da_convergence
```

- [ ] **Step 2: Build to verify compilation**

```bash
cd POD_Fortran && fpm build test_da_convergence
```
Expected: build succeeds.

- [ ] **Step 3: Commit**

```bash
git add app/test_da_convergence.f90
git commit -m "feat: add test_da_convergence skeleton with Phase 0 init"
```

---

### Task 3: Phase 1 & 2 — Read OPM and generate particles

**Files:**
- Modify: `app/test_da_convergence.f90`

- [ ] **Step 1: Add Phase 1 (read OPM) and Phase 2 (generate particles)**

Replace the `stop` statement and add after Phase 0:

```fortran
    ! Local variables
    real(DP) :: nominal_state(6), initial_cov(6,6)
    real(DP) :: epoch0
    integer  :: i

    type(uq_state_type) :: init_state

    real(DP), allocatable :: mc_golden(:,:)

    ! Phase 1: Read OPM
    write(*,*) '>>> Phase 1: Reading OPM from ', OPM_FILE
    call load_initial_opm(OPM_FILE, epoch0, nominal_state, initial_cov)
    write(*,*) '  Epoch (TDB sec):  ', epoch0
    write(*,*) '  Nominal state:    ', nominal_state
    write(*,'(A,6ES12.4)') '  Position sigma (km):  ', &
        sqrt(initial_cov(1,1)), sqrt(initial_cov(2,2)), sqrt(initial_cov(3,3))
    write(*,'(A,6ES12.4)') '  Velocity sigma (km/s):', &
        sqrt(initial_cov(4,4)), sqrt(initial_cov(5,5)), sqrt(initial_cov(6,6))
    write(*,*) ''

    ! Phase 2: Generate 1000 particles
    write(*,*) '>>> Phase 2: Generating ', N_PARTICLES, ' initial particles...'
    call init_state%allocate_memory(6, N_PARTICLES)
    call init_random_seed(.true.)
    call generate_multivariate_normal(nominal_state, initial_cov, init_state%samples)
    init_state%mean = nominal_state
    write(*,*) '  Particles generated. Sample (first):', init_state%samples(:, 1)
    write(*,*) ''

    stop
```

Also add these `use` statements after `use pod_data_format_module`:

```fortran
    use pod_uq_state_module, only: uq_state_type
    use pod_random_module, only: init_random_seed, generate_multivariate_normal
```

- [ ] **Step 2: Build to verify**

```bash
cd POD_Fortran && fpm build test_da_convergence
```
Expected: build succeeds.

- [ ] **Step 3: Commit**

```bash
git add app/test_da_convergence.f90
git commit -m "feat: add OPM reading and particle generation phases"
```

---

### Task 4: Phase 3 — MC golden reference propagation

**Files:**
- Modify: `app/test_da_convergence.f90`

- [ ] **Step 1: Add MC propagation code**

Replace the `stop` after Phase 2 and add:

```fortran
    ! Phase 3: MC golden reference
    write(*,*) '>>> Phase 3: MC propagation (golden reference)...'
    write(*,*) '  Propagating ', N_PARTICLES, ' particles for 15 days with RKF45...'

    block
        type(uq_state_type) :: init_mc, final_mc
        integer :: integrator_rkf45

        integrator_rkf45 = 1  ! METHOD_RKF45

        ! Copy initial particles for MC
        call init_mc%allocate_memory(6, N_PARTICLES)
        init_mc%samples = init_state%samples
        init_mc%mean = nominal_state

        call run_particle_propagation( &
            initial_state   = init_mc, &
            reference_orbit = nominal_state, &
            epoch0          = epoch0, &
            t_start         = 0.0_DP, &
            t_end           = DURATION, &
            method_switch   = METHOD_MC, &
            final_state     = final_mc, &
            integrator_switch = integrator_rkf45)

        ! Save MC golden reference
        allocate(mc_golden(6, N_PARTICLES))
        mc_golden = final_mc%samples

        write(*,*) '  MC propagation complete.'
        write(*,*) '  MC final mean: ', final_mc%mean
        write(*,*) ''
    end block
```

Also add `use pod_uq_propagation, only: run_particle_propagation, METHOD_MC, METHOD_DA` to the use statements.

- [ ] **Step 2: Build to verify**

```bash
cd POD_Fortran && fpm build test_da_convergence
```
Expected: build succeeds.

- [ ] **Step 3: Commit**

```bash
git add app/test_da_convergence.f90
git commit -m "feat: add MC golden reference propagation phase"
```

---

### Task 5: Phase 4 & 5 — DA scan and summary output

**Files:**
- Modify: `app/test_da_convergence.f90`

- [ ] **Step 1: Add DA scan and summary code**

Append after Phase 3 (before `end program`):

```fortran
    ! Phase 4: DA scan over orders 3,4,5,6
    write(*,*) '>>> Phase 4: DA order convergence scan...'
    write(*,*) ''

    block
        type(uq_state_type) :: init_da, final_da
        integer :: order, integrator_rkf45
        real(DP) :: max_rel_err(6), rel_err, abs_err
        integer :: i_p, comp, n_exceed
        logical :: all_pass
        character(len=20) :: status_str

        integrator_rkf45 = 1  ! METHOD_RKF45

        write(*,'(A)') ' Order |   X(%)  |   Y(%)  |   Z(%)  |  VX(%)  |  VY(%)  |  VZ(%)  | Exceed | Status'
        write(*,'(A)') '-------+---------+---------+---------+---------+---------+---------+--------+--------'

        do order = 3, 6
            call init_da%allocate_memory(6, N_PARTICLES)
            init_da%samples = init_state%samples
            init_da%mean = nominal_state

            call run_particle_propagation( &
                initial_state   = init_da, &
                reference_orbit = nominal_state, &
                epoch0          = epoch0, &
                t_start         = 0.0_DP, &
                t_end           = DURATION, &
                method_switch   = METHOD_DA, &
                final_state     = final_da, &
                integrator_switch = integrator_rkf45, &
                da_order        = order)

            ! Compute max relative error per component
            max_rel_err = 0.0_DP
            n_exceed = 0
            do i_p = 1, N_PARTICLES
                do comp = 1, 6
                    abs_err = abs(final_da%samples(comp, i_p) - mc_golden(comp, i_p))
                    if (abs(mc_golden(comp, i_p)) > 1.0e-30_DP) then
                        rel_err = abs_err / abs(mc_golden(comp, i_p))
                        max_rel_err(comp) = max(max_rel_err(comp), rel_err)
                        if (rel_err > 0.01_DP) n_exceed = n_exceed + 1
                    end if
                end do
            end do

            all_pass = all(max_rel_err < 0.01_DP)
            if (all_pass) then
                status_str = 'PASS'
            else
                status_str = 'FAIL'
            end if

            write(*,'(I4,A,6(F7.3,A),I6,A,A)') &
                order, '   |', &
                max_rel_err(1)*100, ' |', &
                max_rel_err(2)*100, ' |', &
                max_rel_err(3)*100, ' |', &
                max_rel_err(4)*100, ' |', &
                max_rel_err(5)*100, ' |', &
                max_rel_err(6)*100, ' |', &
                n_exceed, ' | ', trim(status_str)
        end do
    end block

    ! Phase 5: Summary
    write(*,*) ''
    write(*,*) '============================================================'
    write(*,*) '  Test complete.'
    write(*,*) '============================================================'
```

- [ ] **Step 2: Build to verify**

```bash
cd POD_Fortran && fpm build test_da_convergence
```
Expected: build succeeds.

- [ ] **Step 3: Commit**

```bash
git add app/test_da_convergence.f90
git commit -m "feat: add DA order scan and summary output phases"
```

---

### Task 6: Build and dry-run verification

**Files:**
- (none modified)

- [ ] **Step 1: Full build**

```bash
cd POD_Fortran && fpm build test_da_convergence
```
Expected: build succeeds with no warnings.

- [ ] **Step 2: Verify executable exists**

```bash
ls -la POD_Fortran/build/*/test_da_convergence*
```
Expected: executable binary found.

- [ ] **Step 3: Run (system check only — full 15-day run takes very long)**

```bash
cd POD_Fortran && timeout 30 fpm run test_da_convergence 2>&1 || true
```
Expected: Phase 0 and Phase 1 output visible before timeout. If it gets to Phase 2 particle generation, the initialization is working.

- [ ] **Step 4: Commit any final fixes**

```bash
git add -A  # only if corrections were needed
git commit -m "chore: final tweaks for da convergence test"
```
