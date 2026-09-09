!> @file test_halo_srp_error_da_consistency.f90
!! @brief Verify Real vs DA consistency of the SRP-induced state error.
!!
!! Propagates the L1Halo-1 initial condition for 15 days (about one halo
!! period) under the same four SRP models as test_halo_srp_error_export:
!! cannonball, and box-wing with Sun/Earth/Moon pointing. Each model is run in
!! both the Real and the DA force model, alongside a no-SRP baseline in each
!! framework. The SRP error (model minus no-SRP baseline) at the final epoch
!! is computed in both frameworks and asserted to agree. DA runs with
!! constants only (no independent DA variable), so its nominal part must
!! reproduce Real arithmetic.
program test_halo_srp_error_da_consistency
    use pod_global, only: DP
    use pod_engine_module, only: pod_engine_init
    use pod_config, only: config, validate_config
    use pod_data_format_module, only: load_initial_opm
    use pod_force_model_module, only: &
        set_propagation_epoch_real => set_propagation_epoch, &
        current_epoch0_real => current_epoch0
    use pod_da_force_model_module, only: &
        set_propagation_epoch_da => set_propagation_epoch, &
        current_epoch0_da => current_epoch0, &
        init_gravity_network_da => init_gravity_network, &
        clear_srp_scale_uncertainty, set_srp_ballistic_parameters
    use pod_integrator_module, only: adaptive_step_integrate, METHOD_RKF78
    use pod_da_integrator_module, only: da_adaptive_step_integrate
    use pod_dace_classes, only: AlgebraicVector, dace_initialize, assignment(=)
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
    real(DP), parameter :: POSITION_TOL_KM = 1.0e-2_DP
    real(DP), parameter :: VELOCITY_TOL_KM_S = 1.0e-8_DP

    real(DP), parameter :: SPACECRAFT_MASS_KG = 1200.0_DP
    real(DP), parameter :: CANNONBALL_AREA_M2 = 9.0_DP
    real(DP), parameter :: CANNONBALL_CR = 1.25_DP
    real(DP), parameter :: CANNONBALL_SMR = &
        CANNONBALL_AREA_M2 / SPACECRAFT_MASS_KG
    real(DP), parameter :: SOLAR_PRESSURE_1AU = 1367.0_DP / 299792458.0_DP

    character(len=*), parameter :: CONFIG_FILE = 'config/config.txt'
    character(len=*), parameter :: HALO_OPM = &
        'OPM/L1Halo-1/L1Halo-1_init.opm.json'
    character(len=32), parameter :: CASE_NAMES(N_CASES) = [ &
        character(len=32) :: 'cannonball', 'box-wing / Sun pointing', &
        'box-wing / Earth pointing', 'box-wing / Moon pointing' ]

    real(DP) :: epoch0
    real(DP) :: state0(6), covariance0(6,6)
    real(DP) :: initial_real_nd(6)
    real(DP), allocatable :: times_real_base(:), states_real_base(:,:)
    real(DP), allocatable :: times_da_base(:), nominal_da_base(:,:)
    real(DP), allocatable :: times_real(:), states_real(:,:)
    real(DP), allocatable :: times_da(:), nominal_da(:,:)
    type(AlgebraicVector) :: initial_da, final_da
    real(DP) :: dR_real(3), dV_real(3), dR_da(3), dV_da(3)
    real(DP) :: pos_diff, vel_diff, srp_pos_norm
    integer :: n_real_base, n_da_base, n_real, n_da, case_id, i

    call pod_engine_init(CONFIG_FILE)
    call init_gravity_network_da()
    call load_initial_opm(HALO_OPM, epoch0, state0, covariance0)

    call set_propagation_epoch_real(epoch0)
    call set_propagation_epoch_da(epoch0)
    call assert_true(current_epoch0_real == epoch0, &
                     'Real propagation epoch was not set')
    call assert_true(current_epoch0_da == epoch0, &
                     'DA propagation epoch was not set')

    ! Constants-only DA: no independent variable is created anywhere here.
    call dace_initialize(1, 6)
    call clear_srp_scale_uncertainty()
    call set_srp_ballistic_parameters(Cr=CANNONBALL_CR, &
                                      SMR=CANNONBALL_SMR, &
                                      RP=SOLAR_PRESSURE_1AU)

    initial_real_nd(1:3) = state0(1:3) / config%LU
    initial_real_nd(4:6) = state0(4:6) / config%VU

    call initial_da%init(6)
    do i = 1, 6
        initial_da%elements(i) = initial_real_nd(i)
    end do

    ! No-SRP baselines, one per framework.
    config%use_srp = .false.
    call adaptive_step_integrate( &
        state=initial_real_nd, t_start=0.0_DP, t_end=DURATION_S/config%TU, &
        integrator_method=METHOD_RKF78, times=times_real_base, &
        states=states_real_base, n_steps=n_real_base, max_steps_in=MAX_STEPS, &
        rel_tol_in=REL_TOL, abs_tol_in=ABS_TOL, dt_min_in=DT_MIN_S, &
        dt_max_in=DT_MAX_S)
    call da_adaptive_step_integrate( &
        state=initial_da, t_start=0.0_DP, t_end=DURATION_S/config%TU, &
        integrator_method=METHOD_RKF78, times=times_da_base, &
        nominal_states=nominal_da_base, final_state=final_da, n_steps=n_da_base, &
        max_steps_in=MAX_STEPS, rel_tol_in=REL_TOL, abs_tol_in=ABS_TOL, &
        dt_min_in=DT_MIN_S, dt_max_in=DT_MAX_S)
    call final_da%destroy()
    call assert_true(n_real_base > 1 .and. n_da_base > 1, &
                     'baseline invalid step count')

    write(*,'(a)') 'Real/DA SRP-error consistency (15-day)'
    write(*,'(a)') 'Case                             | |dR_srp| (km) | dR diff (km) | dV diff (km/s)'

    do case_id = 1, N_CASES
        call configure_srp_case(case_id)

        call adaptive_step_integrate( &
            state=initial_real_nd, t_start=0.0_DP, t_end=DURATION_S/config%TU, &
            integrator_method=METHOD_RKF78, times=times_real, states=states_real, &
            n_steps=n_real, max_steps_in=MAX_STEPS, rel_tol_in=REL_TOL, &
            abs_tol_in=ABS_TOL, dt_min_in=DT_MIN_S, dt_max_in=DT_MAX_S)
        call da_adaptive_step_integrate( &
            state=initial_da, t_start=0.0_DP, t_end=DURATION_S/config%TU, &
            integrator_method=METHOD_RKF78, times=times_da, &
            nominal_states=nominal_da, final_state=final_da, n_steps=n_da, &
            max_steps_in=MAX_STEPS, rel_tol_in=REL_TOL, abs_tol_in=ABS_TOL, &
            dt_min_in=DT_MIN_S, dt_max_in=DT_MAX_S)
        call final_da%destroy()

        call assert_true(n_real > 1 .and. n_real <= MAX_STEPS, &
                         trim(CASE_NAMES(case_id))//' invalid Real step count')
        call assert_true(n_da > 1 .and. n_da <= MAX_STEPS, &
                         trim(CASE_NAMES(case_id))//' invalid DA step count')
        call assert_close_scalar(times_real(n_real)*config%TU, DURATION_S, &
                                 FINAL_TIME_TOL_S, &
                                 trim(CASE_NAMES(case_id))//' Real final time')
        call assert_close_scalar(times_da(n_da)*config%TU, DURATION_S, &
                                 FINAL_TIME_TOL_S, &
                                 trim(CASE_NAMES(case_id))//' DA final time')

        ! Final-epoch SRP error, Real framework (state is (n,6)).
        dR_real = (states_real(n_real,1:3) - states_real_base(n_real_base,1:3)) * config%LU
        dV_real = (states_real(n_real,4:6) - states_real_base(n_real_base,4:6)) * config%VU

        ! Final-epoch SRP error, DA framework (nominal_states is (6,n)).
        dR_da = (nominal_da(1:3,n_da) - nominal_da_base(1:3,n_da_base)) * config%LU
        dV_da = (nominal_da(4:6,n_da) - nominal_da_base(4:6,n_da_base)) * config%VU

        pos_diff = maxval(abs(dR_real - dR_da))
        vel_diff = maxval(abs(dV_real - dV_da))
        srp_pos_norm = sqrt(sum(dR_real*dR_real))

        call assert_true(ieee_is_finite(pos_diff) .and. &
                         ieee_is_finite(vel_diff), &
                         trim(CASE_NAMES(case_id))//' non-finite diff')
        call assert_true(pos_diff <= POSITION_TOL_KM, &
                         trim(CASE_NAMES(case_id))// &
                         ' SRP position-error mismatch Real vs DA')
        call assert_true(vel_diff <= VELOCITY_TOL_KM_S, &
                         trim(CASE_NAMES(case_id))// &
                         ' SRP velocity-error mismatch Real vs DA')

        write(*,'(a32,3(2x,es14.6))') CASE_NAMES(case_id), srp_pos_norm, pos_diff, vel_diff

        if (allocated(times_real)) deallocate(times_real)
        if (allocated(states_real)) deallocate(states_real)
        if (allocated(times_da)) deallocate(times_da)
        if (allocated(nominal_da)) deallocate(nominal_da)
    end do

    call initial_da%destroy()
    write(*,'(a)') 'PASS: Real and DA SRP errors agree.'

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
        config%srp_array_front_optical = [0.85_DP, 0.08_DP, 0.07_DP]
        config%srp_array_back_optical = [0.60_DP, 0.20_DP, 0.20_DP]
        config%srp_primary_axis_body = [0.0_DP, 0.0_DP, 1.0_DP]
        config%srp_secondary_axis_body = [0.0_DP, 1.0_DP, 0.0_DP]
        config%srp_roll_reference = 'orbit_normal'
        config%srp_pressure_1au_n_m2 = SOLAR_PRESSURE_1AU
        config%srp_geometry_tolerance = 1.0e-12_DP

        ! Keep the DA force model strictly nominal in all four cases.
        config%srp_scale_da_span = 0.0_DP
        config%srp_attitude_bias_span_arcsec = 0.0_DP
        config%srp_array_angle_span_deg = 0.0_DP
        call clear_srp_scale_uncertainty()

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

end program test_halo_srp_error_da_consistency
