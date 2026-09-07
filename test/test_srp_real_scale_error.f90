!> @file test_srp_real_scale_error.f90
!! @brief Verify the deterministic Real SRP scale override for both models.
!!
!! The override represents a force-model error a=(1+delta_s)*a_nominal.
!! It must not mutate mass, area, optical coefficients, or cannonball inputs.
program test_srp_real_scale_error
    use pod_global, only: DP
    use pod_engine_module, only: pod_engine_init
    use pod_config, only: config, validate_config
    use pod_spice, only: str2et
    use pod_force_model_module, only: compute_solar_radiation_pressure, &
                                      set_srp_scale_error, clear_srp_scale_error
    implicit none

    character(len=*), parameter :: CONFIG_FILE = 'config/config.txt'
    character(len=*), parameter :: TEST_EPOCH = '2027-01-01T00:01:00'
    real(DP), parameter :: CR = 1.25_DP
    real(DP), parameter :: SMR = 7.5e-3_DP
    real(DP), parameter :: RP = 1367.0_DP/299792458.0_DP
    real(DP), parameter :: RTOL = 2.0e-14_DP
    real(DP), parameter :: ATOL = 1.0e-25_DP

    real(DP) :: position(3), velocity(3), epoch
    real(DP) :: nominal(3), plus(3), minus(3), restored(3)

    call pod_engine_init(CONFIG_FILE)
    call str2et(TEST_EPOCH, epoch)
    position = [-315127.0132466179_DP, -99847.53324194982_DP, -78712.98312720712_DP]
    velocity = [0.2469389952601773_DP, -0.670149117406111_DP, -0.2845532515590944_DP]

    call configure_common_geometry()

    ! Cannonball: the override multiplies the already configured physical SRP.
    config%srp_model = 'cannonball'
    call exercise_model('cannonball')

    ! Box-wing: use the same override through the geometry model's scale input.
    config%srp_model = 'box_wing'
    config%srp_attitude_mode = 'sun'
    call assert_true(validate_config(), 'invalid box-wing test configuration')
    call exercise_model('box-wing')

    call clear_srp_scale_error()
    write(*,'(a)') 'PASS: Real SRP scale error is common to cannonball and box-wing.'

contains

    subroutine exercise_model(label)
        character(len=*), intent(in) :: label

        call clear_srp_scale_error()
        call compute_solar_radiation_pressure(position, epoch, nominal, &
                                              Cr=CR, SMR=SMR, RP=RP, velocity=velocity)

        call set_srp_scale_error(0.1_DP)
        call compute_solar_radiation_pressure(position, epoch, plus, &
                                              Cr=CR, SMR=SMR, RP=RP, velocity=velocity)
        call assert_vector_close(plus, 1.1_DP*nominal, label//' +10 percent')

        call set_srp_scale_error(-0.1_DP)
        call compute_solar_radiation_pressure(position, epoch, minus, &
                                              Cr=CR, SMR=SMR, RP=RP, velocity=velocity)
        call assert_vector_close(minus, 0.9_DP*nominal, label//' -10 percent')

        call clear_srp_scale_error()
        call compute_solar_radiation_pressure(position, epoch, restored, &
                                              Cr=CR, SMR=SMR, RP=RP, velocity=velocity)
        call assert_vector_close(restored, nominal, label//' clear restores nominal')
    end subroutine exercise_model

    subroutine configure_common_geometry()
        config%use_srp = .true.
        config%srp_mass_kg = 1200.0_DP
        config%srp_box_dimensions_m = [2.0_DP, 3.0_DP, 4.0_DP]
        config%srp_box_optical = [0.30_DP, 0.40_DP, 0.30_DP]
        config%srp_array_total_area_m2 = 24.0_DP
        config%srp_array_tracking_mode = 'single_axis'
        config%srp_array_hinge_axis_body = [0.0_DP, 1.0_DP, 0.0_DP]
        config%srp_array_reference_normal_body = [1.0_DP, 0.0_DP, 0.0_DP]
        config%srp_array_front_optical = [0.10_DP, 0.80_DP, 0.10_DP]
        config%srp_array_back_optical = [0.60_DP, 0.20_DP, 0.20_DP]
        config%srp_roll_reference = 'orbit_normal'
        config%srp_primary_axis_body = [0.0_DP, 0.0_DP, 1.0_DP]
        config%srp_secondary_axis_body = [0.0_DP, 1.0_DP, 0.0_DP]
        config%srp_pressure_1au_n_m2 = RP
        config%srp_geometry_tolerance = 1.0e-12_DP
        config%srp_scale_da_span = 0.0_DP
        config%srp_attitude_bias_span_arcsec = 0.0_DP
        config%srp_array_angle_span_deg = 0.0_DP
    end subroutine configure_common_geometry

    subroutine assert_vector_close(actual, expected_value, label)
        real(DP), intent(in) :: actual(3), expected_value(3)
        character(len=*), intent(in) :: label
        real(DP) :: error, scale

        error = maxval(abs(actual-expected_value))
        scale = max(maxval(abs(expected_value)), ATOL)
        if (error > ATOL + RTOL*scale) then
            write(*,'(a,2(1x,es24.16))') 'FAIL: '//trim(label)//' error/scale', error, scale
            error stop 1
        end if
    end subroutine assert_vector_close

    subroutine assert_true(condition, label)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: label
        if (.not. condition) then
            write(*,'(a)') 'FAIL: '//trim(label)
            error stop 1
        end if
    end subroutine assert_true

end program test_srp_real_scale_error

