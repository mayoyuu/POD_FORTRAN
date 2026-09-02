!> @file test_box_wing_top_level_propagation.f90
!! @brief Standard config-driven Real/DA box-wing propagation workflow.
!!
!! This test demonstrates the intended public API sequence:
!!
!!   1. initialize POD and load a physical Halo initial state;
!!   2. set config%srp_model='box_wing' and the spacecraft geometry;
!!   3. call propagate_orbit for the Real trajectory;
!!   4. call propagate_da_orbit for the DA state map;
!!   5. compare the DA nominal trajectory with the Real trajectory;
!!   6. clean both result objects.
!!
!! The orbitprop wrappers do not contain a second SRP implementation. During
!! every RHS evaluation they call the corresponding force-model module, which
!! reads the shared global config and dispatches to cannonball or box-wing.
!!
!! The DA wrapper intentionally creates six initial-state DA variables because
!! that is its standard state-map/STM interface. No SRP model uncertainty is
!! enabled here: scale, attitude, and array-angle spans are all zero. Only the
!! constant term of the DA result is compared with the Real solution.
program test_box_wing_top_level_propagation
    use pod_global, only: DP
    use pod_engine_module, only: pod_engine_init
    use pod_config, only: config, validate_config
    use pod_data_format_module, only: load_initial_opm
    use pod_orbit_propagation, only: orbit_state, propagation_result, &
        propagate_orbit, cleanup_propagation_result
    use pod_da_orbit_propagation, only: da_orbit_state, &
        da_propagation_result, propagate_da_orbit, &
        cleanup_da_propagation_result
    use pod_da_force_model_module, only: init_gravity_network_da => &
        init_gravity_network, cleanup_gravity_network_da => &
        cleanup_gravity_network, clear_srp_scale_uncertainty
    use pod_dace_classes, only: active_da_count
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none

    real(DP), parameter :: DURATION_S = 86400.0_DP
    integer, parameter :: RKF78_CHOICE = 2
    real(DP), parameter :: POSITION_TOL_KM = 1.0e-3_DP
    real(DP), parameter :: VELOCITY_TOL_KM_S = 1.0e-9_DP
    real(DP), parameter :: MIN_MODEL_SEPARATION_KM = 1.0e-6_DP
    character(len=*), parameter :: CONFIG_FILE = 'config/config.txt'
    character(len=*), parameter :: HALO_OPM = &
        'OPM/L1Halo-1/L1Halo-1_init.opm.json'

    real(DP) :: epoch0, state0(6), covariance0(6,6)
    real(DP) :: real_final(6), da_final(6), cannonball_final(6)
    real(DP) :: position_error, velocity_error, model_separation
    type(orbit_state) :: real_initial
    type(da_orbit_state) :: da_initial
    type(propagation_result) :: real_result, cannonball_result
    type(da_propagation_result) :: da_result

    call pod_engine_init(CONFIG_FILE)
    call init_gravity_network_da()
    call load_initial_opm(HALO_OPM, epoch0, state0, covariance0)

    call configure_box_wing()

    ! Real and DA use the same physical initial state and absolute TDB epoch.
    real_initial%state = state0
    real_initial%epoch = epoch0
    da_initial%nominal_state = state0
    da_initial%epoch = epoch0
    da_initial%da_order = 1

    ! This is the standard public propagation path. No low-level integrator or
    ! force-model routine is called directly by the test.
    call propagate_orbit(real_initial, DURATION_S, RKF78_CHOICE, real_result)
    call propagate_da_orbit(da_initial, DURATION_S, RKF78_CHOICE, da_result)

    real_final = real_result%states(real_result%n_steps, :)
    da_final = da_result%nominal_states(:, da_result%n_steps)
    position_error = maxval(abs(real_final(1:3) - da_final(1:3)))
    velocity_error = maxval(abs(real_final(4:6) - da_final(4:6)))

    ! A cannonball reference propagated through the same public Real wrapper
    ! proves that changing config%srp_model changes the selected force model.
    config%srp_model = 'cannonball'
    call propagate_orbit(real_initial, DURATION_S, RKF78_CHOICE, &
                         cannonball_result)
    cannonball_final = cannonball_result%states(cannonball_result%n_steps, :)
    model_separation = maxval(abs(real_final(1:3) - cannonball_final(1:3)))
    config%srp_model = 'box_wing'

    write(*,'(a)') 'Top-level config-driven box-wing propagation'
    write(*,'(a,a)') '  selected SRP model   : ', trim(config%srp_model)
    write(*,'(a,a)') '  box-wing attitude    : ', &
        trim(config%srp_attitude_mode)
    call print_state('  Real box-wing final state', real_final)
    call print_state('  DA nominal box-wing final state', da_final)
    call print_state('  Real cannonball reference state', cannonball_final)
    call print_state('  Real - DA difference', real_final - da_final)
    write(*,'(a,es14.6)') '  max position difference [km]   : ', &
        position_error
    write(*,'(a,es14.6)') '  max velocity difference [km/s] : ', &
        velocity_error
    write(*,'(a,es14.6)') '  box-wing/cannonball separation [km]: ', &
        model_separation

    call assert_true(trim(config%srp_model) == 'box_wing', &
        'canonical box-wing model was not retained')
    call assert_true(real_result%n_steps > 1 .and. da_result%n_steps > 1, &
        'top-level propagation produced too few steps')
    call assert_true(all(ieee_is_finite(real_final)), &
        'Real box-wing final state is not finite')
    call assert_true(all(ieee_is_finite(da_final)), &
        'DA box-wing nominal final state is not finite')
    call assert_true(abs(real_result%times(real_result%n_steps) - &
                          DURATION_S) <= 1.0e-6_DP, &
        'Real wrapper did not reach the requested final time')
    call assert_true(abs(da_result%times(da_result%n_steps) - &
                          DURATION_S) <= 1.0e-6_DP, &
        'DA wrapper did not reach the requested final time')
    call assert_true(position_error <= POSITION_TOL_KM, &
        'top-level Real/DA box-wing position mismatch')
    call assert_true(velocity_error <= VELOCITY_TOL_KM_S, &
        'top-level Real/DA box-wing velocity mismatch')
    call assert_true(model_separation > MIN_MODEL_SEPARATION_KM, &
        'box-wing selection appears to have fallen back to cannonball')

    call cleanup_propagation_result(real_result)
    call cleanup_propagation_result(cannonball_result)
    call cleanup_da_propagation_result(da_result)
    call cleanup_gravity_network_da()
    call assert_true(active_da_count() == 0, &
        'top-level DA propagation left active DA handles')

    write(*,'(a)') 'PASS: top-level Real/DA box-wing workflow is valid.'

contains

    !> Configure a complete Sun-pointing box-wing spacecraft.
    !! The canonical internal selector is box_wing. Neither box-wing nor
    !! boxwing is accepted as an internal model identifier.
    subroutine configure_box_wing()
        config%use_srp = .true.
        config%srp_model = 'box_wing'
        config%srp_attitude_mode = 'sun'
        config%srp_roll_reference = 'orbit_normal'
        config%srp_primary_axis_body = [0.0_DP, 0.0_DP, 1.0_DP]
        config%srp_secondary_axis_body = [0.0_DP, 1.0_DP, 0.0_DP]

        config%srp_mass_kg = 1200.0_DP
        config%srp_box_dimensions_m = [2.0_DP, 3.0_DP, 4.0_DP]
        config%srp_box_optical = [0.30_DP, 0.40_DP, 0.30_DP]
        config%srp_array_total_area_m2 = 24.0_DP
        config%srp_array_tracking_mode = 'single_axis'
        config%srp_array_hinge_axis_body = [0.0_DP, 1.0_DP, 0.0_DP]
        config%srp_array_reference_normal_body = [1.0_DP, 0.0_DP, 0.0_DP]
        config%srp_array_front_optical = [0.10_DP, 0.80_DP, 0.10_DP]
        config%srp_array_back_optical = [0.60_DP, 0.20_DP, 0.20_DP]
        config%srp_pressure_1au_n_m2 = 1367.0_DP / 299792458.0_DP
        config%srp_geometry_tolerance = 1.0e-12_DP

        ! Disable every SRP-specific DA uncertainty source. The six variables
        ! created by propagate_da_orbit belong only to its state-map interface.
        config%srp_scale_da_span = 0.0_DP
        config%srp_attitude_bias_span_arcsec = 0.0_DP
        config%srp_array_angle_span_deg = 0.0_DP
        call clear_srp_scale_uncertainty()

        call assert_true(validate_config(), &
            'invalid box-wing configuration for top-level propagation')
    end subroutine configure_box_wing

    !> Print a Cartesian state using its physical position and velocity units.
    subroutine print_state(label, state)
        character(len=*), intent(in) :: label
        real(DP), intent(in) :: state(6)

        write(*,'(a)') trim(label)
        write(*,'(a,3(1x,es24.16))') '    r [km]   =', state(1:3)
        write(*,'(a,3(1x,es24.16))') '    v [km/s] =', state(4:6)
    end subroutine print_state

    subroutine assert_true(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (.not. condition) then
            write(*,'(a)') 'FAIL: '//trim(message)
            error stop 1
        end if
    end subroutine assert_true

end program test_box_wing_top_level_propagation
