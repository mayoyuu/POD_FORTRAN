program test_spacecraft_geometry_config
    use pod_global, only: DP
    use pod_config, only: config, set_default_config, load_config, validate_config
    implicit none

    character(len=*), parameter :: temp_file = 'test_box_wing_config.tmp'
    integer :: unit

    call set_default_config()
    call assert_true(trim(config%srp_model) == 'cannonball', 'default SRP model')
    call assert_close(config%srp_primary_axis_body, [0.0_DP, 0.0_DP, 1.0_DP], 'default primary axis')
    call assert_close(config%srp_secondary_axis_body, [0.0_DP, 1.0_DP, 0.0_DP], 'default secondary axis')
    call assert_true(validate_config(), 'cannonball does not require box-wing geometry')

    open(newunit=unit, file=temp_file, status='replace', action='write')
    write(unit, '(A)') 'srp_model = box_wing'
    write(unit, '(A)') 'srp_attitude_mode = moon'
    write(unit, '(A)') 'srp_roll_reference = orbit_normal'
    write(unit, '(A)') 'srp_primary_axis_body = 0 0 1'
    write(unit, '(A)') 'srp_secondary_axis_body = 0 1 0'
    write(unit, '(A)') 'srp_mass_kg = 1200'
    write(unit, '(A)') 'srp_box_dimensions_m = 2 3 4'
    write(unit, '(A)') 'srp_box_optical = 0.30 0.40 0.30'
    write(unit, '(A)') 'srp_array_total_area_m2 = 24'
    write(unit, '(A)') 'srp_array_tracking_mode = single_axis'
    write(unit, '(A)') 'srp_array_hinge_axis_body = 0 1 0'
    write(unit, '(A)') 'srp_array_reference_normal_body = 1 0 0'
    write(unit, '(A)') 'srp_array_front_optical = 0.10 0.80 0.10'
    write(unit, '(A)') 'srp_array_back_optical = 0.60 0.20 0.20'
    write(unit, '(A)') 'srp_pressure_1au_n_m2 = 4.56e-6'
    write(unit, '(A)') 'srp_geometry_tolerance = 1.0e-12'
    write(unit, '(A)') 'srp_scale_da_span = 0.20'
    write(unit, '(A)') 'srp_attitude_bias_span_arcsec = 10 20 30'
    write(unit, '(A)') 'srp_array_angle_span_deg = 2'
    close(unit)

    call load_config(temp_file)
    call assert_true(trim(config%srp_model) == 'box_wing', 'parsed SRP model')
    call assert_true(trim(config%srp_attitude_mode) == 'moon', 'parsed attitude mode')
    call assert_close(config%srp_box_dimensions_m, [2.0_DP, 3.0_DP, 4.0_DP], 'parsed box dimensions')
    call assert_close(config%srp_attitude_bias_span_arcsec, [10.0_DP, 20.0_DP, 30.0_DP], &
                      'parsed attitude spans')
    call assert_true(validate_config(), 'valid box-wing configuration')

    config%srp_mass_kg = -1.0_DP
    call assert_true(.not. validate_config(), 'negative mass rejected')
    config%srp_mass_kg = 1200.0_DP

    config%srp_secondary_axis_body = config%srp_primary_axis_body
    call assert_true(.not. validate_config(), 'parallel primary and secondary axes rejected')
    config%srp_secondary_axis_body = [0.0_DP, 1.0_DP, 0.0_DP]

    config%srp_box_optical = [0.3_DP, 0.4_DP, 0.4_DP]
    call assert_true(.not. validate_config(), 'optical coefficients must sum to one')

    open(newunit=unit, file=temp_file, status='old')
    close(unit, status='delete')
    write(*,*) 'Spacecraft geometry configuration tests passed.'

contains

    subroutine assert_true(condition, label)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: label
        if (.not. condition) then
            write(*,*) 'FAILED: ', trim(label)
            stop 1
        end if
    end subroutine assert_true

    subroutine assert_close(actual, expected, label)
        real(DP), intent(in) :: actual(3), expected(3)
        character(len=*), intent(in) :: label
        if (maxval(abs(actual - expected)) > 1.0e-14_DP) then
            write(*,*) 'FAILED: ', trim(label), actual, expected
            stop 1
        end if
    end subroutine assert_close

end program test_spacecraft_geometry_config
