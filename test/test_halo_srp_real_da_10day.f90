!> @file test_halo_srp_real_da_10day.f90
!! @brief Ten-day Real/DA consistency test for cannonball and box-wing SRP.
!!
!! This deterministic integration test propagates the repository's L1Halo-1
!! initial condition for ten days under four SRP configurations: cannonball,
!! and box-wing with Sun, Earth, or Moon pointing.
!!
!! The Real and DA propagators use the same state, epoch, force model, RKF78
!! method, tolerances, and step limits. The DA state contains constants only:
!! DACE is initialized, but this test never creates an independent DA variable.
!! SRP scale, attitude, and array-angle uncertainty spans are explicitly zero.
!! This checks nominal DA arithmetic against Real arithmetic; it is not an
!! uncertainty propagation test and does not compute a state transition matrix.
!!
!! Spacecraft definition shared by all cases:
!!   mass                    = 1200 kg
!!   box dimensions          = 2 m x 3 m x 4 m
!!   total solar-array area  = 24 m^2
!!
!! A cannonball needs one fixed reference area instead of attitude-dependent
!! projected areas. Its selected area is 9 m^2, so A/m = 9/1200 = 7.5e-3
!! m^2/kg, matching the current Real cannonball default. The DA ballistic
!! parameters are set identically. Real and DA are compared within each model;
!! cannonball and box-wing trajectories are not required to agree.
program test_halo_srp_real_da_10day
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
    use pod_dace_classes, only: AlgebraicVector, dace_initialize, &
        active_da_count, assignment(=)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none

    integer, parameter :: N_CASES = 4
    integer, parameter :: MAX_STEPS = 10000
    real(DP), parameter :: DAY_S = 86400.0_DP
    real(DP), parameter :: DURATION_S = 10.0_DP * DAY_S
    real(DP), parameter :: REL_TOL = 1.0e-12_DP
    real(DP), parameter :: ABS_TOL = 1.0e-12_DP
    real(DP), parameter :: DT_MIN_S = 1.0e-6_DP
    real(DP), parameter :: DT_MAX_S = 3600.0_DP

    real(DP), parameter :: SPACECRAFT_MASS_KG = 1200.0_DP
    real(DP), parameter :: CANNONBALL_AREA_M2 = 9.0_DP
    real(DP), parameter :: CANNONBALL_CR = 1.25_DP
    real(DP), parameter :: CANNONBALL_SMR = &
        CANNONBALL_AREA_M2 / SPACECRAFT_MASS_KG
    real(DP), parameter :: SOLAR_PRESSURE_1AU = &
        1367.0_DP / 299792458.0_DP

    ! Position and velocity need different physical-unit tolerances. These
    ! limits correspond to 1 m and 1 micrometre/s.
    real(DP), parameter :: POSITION_TOL_KM = 1.0e-3_DP
    real(DP), parameter :: VELOCITY_TOL_KM_S = 1.0e-9_DP
    real(DP), parameter :: FINAL_TIME_TOL_S = 1.0e-6_DP

    character(len=*), parameter :: CONFIG_FILE = 'config/config.txt'
    character(len=*), parameter :: HALO_OPM = &
        'OPM/L1Halo-1/L1Halo-1_init.opm.json'
    character(len=32), parameter :: CASE_NAMES(N_CASES) = [ &
        character(len=32) :: 'cannonball', 'box-wing / Sun pointing', &
        'box-wing / Earth pointing', 'box-wing / Moon pointing' ]

    real(DP) :: epoch0
    real(DP) :: state0(6), covariance0(6,6)
    real(DP) :: initial_real_nd(6)
    real(DP) :: final_real(6), final_da_physical(6), final_da_nd(6)
    real(DP) :: position_error, velocity_error
    real(DP), allocatable :: times_real(:), states_real(:,:)
    real(DP), allocatable :: times_da(:), states_da(:,:)
    type(AlgebraicVector) :: initial_da, final_da
    integer :: n_real, n_da, case_id, i
    integer :: handles_before

    ! Initialize production configuration, SPICE, and both gravity networks.
    ! The OPM reader converts UTC epoch to TDB and returns km and km/s.
    call pod_engine_init(CONFIG_FILE)
    call init_gravity_network_da()
    call load_initial_opm(HALO_OPM, epoch0, state0, covariance0)

    call set_propagation_epoch_real(epoch0)
    call set_propagation_epoch_da(epoch0)
    call assert_true(current_epoch0_real == epoch0, &
                     'Real propagation epoch was not set')
    call assert_true(current_epoch0_da == epoch0, &
                     'DA propagation epoch was not set')

    ! A maximum variable count does not create DA variables. Only da_var or
    ! init_var would do that, and neither is called anywhere in this test.
    call dace_initialize(1, 6)
    call clear_srp_scale_uncertainty()
    call set_srp_ballistic_parameters(Cr=CANNONBALL_CR, &
                                      SMR=CANNONBALL_SMR, &
                                      RP=SOLAR_PRESSURE_1AU)
    handles_before = active_da_count()
    write(*,'(a,i0)') 'DA handles at baseline: ', handles_before

    ! Both low-level integrators work in repository nondimensional units.
    initial_real_nd(1:3) = state0(1:3) / config%LU
    initial_real_nd(4:6) = state0(4:6) / config%VU

    call initial_da%init(6)
    do i = 1, 6
        initial_da%elements(i) = initial_real_nd(i)
    end do
    write(*,'(a,i0)') 'DA handles after constant initial state: ', &
        active_da_count()

    write(*,'(a)') '10-day L1 Halo SRP Real/DA consistency test'
    write(*,'(a)') 'Case                             |dR|_inf [km]   |dV|_inf [km/s]'

    do case_id = 1, N_CASES
        call configure_srp_case(case_id)

        ! Adaptive step histories may differ, so compare at the common final
        ! epoch rather than requiring equal accepted-step counts.
        call adaptive_step_integrate( &
            state=initial_real_nd, t_start=0.0_DP, &
            t_end=DURATION_S/config%TU, integrator_method=METHOD_RKF78, &
            times=times_real, states=states_real, n_steps=n_real, &
            max_steps_in=MAX_STEPS, rel_tol_in=REL_TOL, abs_tol_in=ABS_TOL, &
            dt_min_in=DT_MIN_S, dt_max_in=DT_MAX_S)

        call da_adaptive_step_integrate( &
            state=initial_da, t_start=0.0_DP, &
            t_end=DURATION_S/config%TU, integrator_method=METHOD_RKF78, &
            times=times_da, nominal_states=states_da, final_state=final_da, &
            n_steps=n_da, max_steps_in=MAX_STEPS, rel_tol_in=REL_TOL, &
            abs_tol_in=ABS_TOL, dt_min_in=DT_MIN_S, dt_max_in=DT_MAX_S)

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

        ! Restore physical units before comparing nominal states.
        final_real(1:3) = states_real(n_real,1:3) * config%LU
        final_real(4:6) = states_real(n_real,4:6) * config%VU
        final_da_nd = final_da%cons()
        final_da_physical(1:3) = final_da_nd(1:3) * config%LU
        final_da_physical(4:6) = final_da_nd(4:6) * config%VU

        call assert_true(all(ieee_is_finite(final_real)), &
                         trim(CASE_NAMES(case_id))//' Real state is not finite')
        call assert_true(all(ieee_is_finite(final_da_physical)), &
                         trim(CASE_NAMES(case_id))//' DA state is not finite')

        position_error = maxval(abs(final_real(1:3) - &
                                    final_da_physical(1:3)))
        velocity_error = maxval(abs(final_real(4:6) - &
                                    final_da_physical(4:6)))
        write(*,'(a32,2(2x,es14.6))') CASE_NAMES(case_id), &
            position_error, velocity_error
        call print_final_state('  Real final state', final_real)
        call print_final_state('  DA nominal final state', final_da_physical)
        call print_final_state('  Real - DA difference', &
                               final_real - final_da_physical)

        call assert_true(position_error <= POSITION_TOL_KM, &
                         trim(CASE_NAMES(case_id))// &
                         ' Real/DA position mismatch after ten days')
        call assert_true(velocity_error <= VELOCITY_TOL_KM_S, &
                         trim(CASE_NAMES(case_id))// &
                         ' Real/DA velocity mismatch after ten days')

        call cleanup_case_outputs()
        write(*,'(a,a,a,i0)') 'DA handles after ', &
            trim(CASE_NAMES(case_id)), ': ', active_da_count()
    end do

    call initial_da%destroy()
    write(*,'(a,i0)') 'DA handles after final cleanup: ', active_da_count()
    call assert_true(active_da_count() == handles_before, &
                     'DA handle leak in ten-day SRP integration test')
    write(*,'(a)') 'PASS: all four ten-day nominal Real/DA SRP cases agree.'

contains

    !> Configure one deterministic SRP case.
    !! Geometry and optical values are explicit. Only the pointing target
    !! changes among the three box-wing cases.
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

    !> Release all per-case allocations and DA handles.
    subroutine cleanup_case_outputs()
        if (allocated(times_real)) deallocate(times_real)
        if (allocated(states_real)) deallocate(states_real)
        if (allocated(times_da)) deallocate(times_da)
        if (allocated(states_da)) deallocate(states_da)
        call final_da%destroy()
    end subroutine cleanup_case_outputs

    !> Print one Cartesian state with position and velocity units separated.
    !! Keeping the six components visible makes the long integration test useful
    !! both as an automated assertion and as a reproducible numerical reference.
    subroutine print_final_state(label, state)
        character(len=*), intent(in) :: label
        real(DP), intent(in) :: state(6)

        write(*,'(a)') trim(label)
        write(*,'(a,3(1x,es24.16))') '    r [km]    =', state(1:3)
        write(*,'(a,3(1x,es24.16))') '    v [km/s]  =', state(4:6)
    end subroutine print_final_state

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

end program test_halo_srp_real_da_10day
