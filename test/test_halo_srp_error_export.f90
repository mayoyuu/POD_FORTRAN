!> @file test_halo_srp_error_export.f90
!! @brief Export SRP-induced state error over one halo period (15 days).
!!
!! Propagates the L1Halo-1 initial condition for 15 days (about one halo
!! period) with SRP disabled (no-SRP baseline) and under four SRP models:
!! cannonball, and box-wing with Sun/Earth/Moon pointing. The error of each
!! model relative to the no-SRP baseline is sampled at every accepted
!! integrator step and written to CSV files under SRP/260902_state_err/.
program test_halo_srp_error_export
    use pod_global, only: DP
    use pod_engine_module, only: pod_engine_init
    use pod_config, only: config, validate_config
    use pod_data_format_module, only: load_initial_opm
    use pod_force_model_module, only: set_propagation_epoch
    use pod_integrator_module, only: adaptive_step_integrate, METHOD_RKF78
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none

    integer, parameter :: N_CASES = 4
    integer, parameter :: MAX_STEPS = 10000
    real(DP), parameter :: DAY_S = 86400.0_DP
    real(DP), parameter :: DURATION_S = 15.0_DP * DAY_S
    real(DP), parameter :: REL_TOL = 1.0e-12_DP
    real(DP), parameter :: ABS_TOL = 1.0e-12_DP
    real(DP), parameter :: DT_MIN_S = 1.0e-6_DP
    real(DP), parameter :: DT_MAX_S = 3600.0_DP
    real(DP), parameter :: FINAL_TIME_TOL_S = 1.0e-6_DP

    real(DP), parameter :: SPACECRAFT_MASS_KG = 1200.0_DP
    real(DP), parameter :: SOLAR_PRESSURE_1AU = 1367.0_DP / 299792458.0_DP

    character(len=*), parameter :: CONFIG_FILE = 'config/config.txt'
    character(len=*), parameter :: HALO_OPM = &
        'OPM/L1Halo-1/L1Halo-1_init.opm.json'
    character(len=*), parameter :: OUT_DIR = 'SRP/260902_state_err'

    character(len=32), parameter :: CASE_NAMES(N_CASES) = [ &
        character(len=32) :: 'cannonball', 'box-wing / Sun pointing', &
        'box-wing / Earth pointing', 'box-wing / Moon pointing' ]
    character(len=32), parameter :: CASE_TAGS(N_CASES) = [ &
        character(len=32) :: 'cannonball', 'boxwing_sun', &
        'boxwing_earth', 'boxwing_moon' ]

    real(DP) :: epoch0
    real(DP) :: state0(6), covariance0(6,6)
    real(DP) :: initial_nd(6)
    real(DP), allocatable :: times_base(:), states_base(:,:)
    real(DP), allocatable :: times_model(:), states_model(:,:)
    real(DP) :: base_interp(6), dr(3), dv(3), dr_norm, dv_norm, time_days
    integer :: n_base, n_model, case_id, i, u
    character(len=256) :: fname

    call pod_engine_init(CONFIG_FILE)
    call load_initial_opm(HALO_OPM, epoch0, state0, covariance0)
    call set_propagation_epoch(epoch0)

    initial_nd(1:3) = state0(1:3) / config%LU
    initial_nd(4:6) = state0(4:6) / config%VU

    call execute_command_line('mkdir -p ' // OUT_DIR)

    ! 1. No-SRP baseline trajectory.
    config%use_srp = .false.
    call adaptive_step_integrate( &
        state=initial_nd, t_start=0.0_DP, t_end=DURATION_S/config%TU, &
        integrator_method=METHOD_RKF78, times=times_base, states=states_base, &
        n_steps=n_base, max_steps_in=MAX_STEPS, rel_tol_in=REL_TOL, &
        abs_tol_in=ABS_TOL, dt_min_in=DT_MIN_S, dt_max_in=DT_MAX_S)
    call assert_true(n_base > 1 .and. n_base <= MAX_STEPS, &
                     'baseline invalid step count')
    call write_baseline_csv(times_base, states_base, n_base)

    write(*,'(a)') '15-day L1 Halo SRP state-error export'
    write(*,'(a)') 'Case                             |   steps | dR_final (km)'

    ! 2. Four SRP models, each sampled at its own accepted steps.
    do case_id = 1, N_CASES
        call configure_srp_case(case_id)

        call adaptive_step_integrate( &
            state=initial_nd, t_start=0.0_DP, t_end=DURATION_S/config%TU, &
            integrator_method=METHOD_RKF78, times=times_model, &
            states=states_model, n_steps=n_model, max_steps_in=MAX_STEPS, &
            rel_tol_in=REL_TOL, abs_tol_in=ABS_TOL, dt_min_in=DT_MIN_S, &
            dt_max_in=DT_MAX_S)
        call assert_true(n_model > 1 .and. n_model <= MAX_STEPS, &
                         trim(CASE_NAMES(case_id))//' invalid step count')
        call assert_close_scalar(times_model(n_model)*config%TU, DURATION_S, &
                                 FINAL_TIME_TOL_S, &
                                 trim(CASE_NAMES(case_id))//' final time')

        fname = OUT_DIR // '/halo_srp_error_' // trim(CASE_TAGS(case_id)) // '.csv'
        open(newunit=u, file=trim(fname), status='replace', action='write')
        write(u,'(a)') 'time_days,dR_x_km,dR_y_km,dR_z_km,dR_norm_km,'// &
                       'dV_x_kms,dV_y_kms,dV_z_kms,dV_norm_kms'

        do i = 1, n_model
            time_days = times_model(i) * config%TU / DAY_S
            call interp_state(times_model(i), times_base, states_base, &
                              n_base, base_interp)
            dr = (states_model(i,1:3) - base_interp(1:3)) * config%LU
            dv = (states_model(i,4:6) - base_interp(4:6)) * config%VU
            dr_norm = sqrt(sum(dr*dr))
            dv_norm = sqrt(sum(dv*dv))

            call assert_true(ieee_is_finite(dr_norm) .and. &
                             ieee_is_finite(dv_norm), &
                             trim(CASE_NAMES(case_id))//' non-finite error')

            write(u,'(es16.8,",",es16.8,",",es16.8,",",es16.8,",",es16.8,",",es16.8,",",es16.8,",",es16.8,",",es16.8)') &
                time_days, dr(1), dr(2), dr(3), dr_norm, &
                dv(1), dv(2), dv(3), dv_norm
        end do
        close(u)

        write(*,'(a32,2x,i6,2x,es14.6)') CASE_NAMES(case_id), n_model, dr_norm

        if (allocated(times_model)) deallocate(times_model)
        if (allocated(states_model)) deallocate(states_model)
    end do

    write(*,'(a)') 'PASS: halo SRP state-error export complete.'
    write(*,'(a)') 'Output written to ' // OUT_DIR // '/'

contains

    !> Configure one deterministic SRP case (cannonball or box-wing).
    subroutine configure_srp_case(id)
        integer, intent(in) :: id

        config%use_srp = .true.
        config%srp_mass_kg = SPACECRAFT_MASS_KG
        config%srp_box_dimensions_m = [2.0_DP, 3.0_DP, 4.0_DP]
        config%srp_box_optical = [0.30_DP, 0.40_DP, 0.30_DP]
        config%srp_array_total_area_m2 = 24.0_DP
        config%srp_array_tracking_mode = 'single_axis'
        config%srp_array_hinge_axis_body = [0.0_DP, 1.0_DP, 0.0_DP]
        config%srp_array_reference_normal_body = [1.0_DP, 0.0_DP, 0.0_DP]
        config%srp_array_front_optical = [0.10_DP, 0.80_DP, 0.10_DP]
        config%srp_array_back_optical = [0.60_DP, 0.20_DP, 0.20_DP]
        config%srp_primary_axis_body = [0.0_DP, 0.0_DP, 1.0_DP]
        config%srp_secondary_axis_body = [0.0_DP, 1.0_DP, 0.0_DP]
        config%srp_roll_reference = 'orbit_normal'
        config%srp_pressure_1au_n_m2 = SOLAR_PRESSURE_1AU
        config%srp_geometry_tolerance = 1.0e-12_DP

        select case (id)
        case (1)
            config%srp_model = 'cannonball'
            config%srp_attitude_mode = 'sun'
        case (2)
            config%srp_model = 'box_wing'
            config%srp_attitude_mode = 'sun'
        case (3)
            config%srp_model = 'box_wing'
            config%srp_attitude_mode = 'earth'
        case (4)
            config%srp_model = 'box_wing'
            config%srp_attitude_mode = 'moon'
        case default
            error stop 'invalid SRP case identifier'
        end select

        call assert_true(validate_config(), &
                         trim(CASE_NAMES(id))//' has invalid configuration')
    end subroutine configure_srp_case

    !> Write the no-SRP baseline trajectory to CSV.
    subroutine write_baseline_csv(times, states, n)
        real(DP), intent(in) :: times(:), states(:,:)
        integer, intent(in) :: n
        integer :: i, u

        open(newunit=u, file=OUT_DIR // '/halo_srp_baseline_nosrp.csv', &
             status='replace', action='write')
        write(u,'(a)') 'time_days,x_km,y_km,z_km,vx_kms,vy_kms,vz_kms'
        do i = 1, n
            write(u,'(es16.8,",",es16.8,",",es16.8,",",es16.8,",",es16.8,",",es16.8,",",es16.8)') &
                times(i)*config%TU/DAY_S, &
                states(i,1)*config%LU, states(i,2)*config%LU, states(i,3)*config%LU, &
                states(i,4)*config%VU, states(i,5)*config%VU, states(i,6)*config%VU
        end do
        close(u)
    end subroutine write_baseline_csv

    !> Linear interpolation of a 6-state trajectory onto query time tq.
    subroutine interp_state(tq, times, states, n, out)
        real(DP), intent(in) :: tq
        real(DP), intent(in) :: times(:), states(:,:)
        integer, intent(in) :: n
        real(DP), intent(out) :: out(6)
        integer :: lo, hi, mid
        real(DP) :: w

        if (tq <= times(1)) then
            out = states(1,:)
            return
        end if
        if (tq >= times(n)) then
            out = states(n,:)
            return
        end if

        lo = 1
        hi = n
        do while (hi - lo > 1)
            mid = (lo + hi) / 2
            if (times(mid) <= tq) then
                lo = mid
            else
                hi = mid
            end if
        end do

        w = (tq - times(lo)) / (times(hi) - times(lo))
        out(:) = states(lo,:) + w * (states(hi,:) - states(lo,:))
    end subroutine interp_state

    !> Assert a logical condition and retain a useful failure message.
    subroutine assert_true(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (.not. condition) then
            write(*,'(a)') 'FAIL: '//trim(message)
            error stop 1
        end if
    end subroutine assert_true

    !> Assert scalar equality within an absolute tolerance.
    subroutine assert_close_scalar(actual, expected, tolerance, message)
        real(DP), intent(in) :: actual, expected, tolerance
        character(len=*), intent(in) :: message

        if (abs(actual - expected) > tolerance) then
            write(*,'(a,2(1x,es24.16),a,es12.4)') &
                'FAIL: '//trim(message)//' actual/expected:', &
                actual, expected, ' tolerance=', tolerance
            error stop 1
        end if
    end subroutine assert_close_scalar

end program test_halo_srp_error_export
