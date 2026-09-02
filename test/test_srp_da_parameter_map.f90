program test_srp_da_parameter_map
    use pod_global, only: DP
    use pod_config, only: config, set_default_config
    use pod_dace_classes, only: dace_initialize
    use pod_spacecraft_geometry, only: srp_da_parameter_map_type
    use pod_da_force_model_module, only: set_srp_scale_uncertainty, clear_srp_scale_uncertainty, &
                                         build_srp_da_parameter_map
    implicit none

    type(srp_da_parameter_map_type) :: map
    real(DP), parameter :: arcsec_to_rad = acos(-1.0_DP)/(180.0_DP*3600.0_DP)
    real(DP), parameter :: deg_to_rad = acos(-1.0_DP)/180.0_DP

    call set_default_config()
    config%srp_attitude_bias_span_arcsec = [10.0_DP, 20.0_DP, 30.0_DP]
    config%srp_array_angle_span_deg = 2.0_DP

    call dace_initialize(2, 6)
    call clear_srp_scale_uncertainty()
    call build_srp_da_parameter_map(map)
    call assert_true(map%global_scale_index == 0, '6D has no scale variable')
    call assert_true(all(map%attitude_bias_index == 0), '6D has no attitude variables')
    call assert_true(map%array_angle_index == 0, '6D has no array variable')

    call dace_initialize(2, 7)
    call set_srp_scale_uncertainty(7, 0.1_DP, 0.2_DP)
    call build_srp_da_parameter_map(map)
    call assert_true(map%global_scale_index == 7, '7D scale index')
    call assert_close(map%global_scale_nominal, 0.1_DP, 'scale nominal')
    call assert_close(map%global_scale_span, 0.2_DP, 'scale span')

    call dace_initialize(2, 10)
    call build_srp_da_parameter_map(map)
    call assert_true(all(map%attitude_bias_index == [8,9,10]), '10D attitude indices')
    call assert_vector_close(map%attitude_bias_span_rad, &
                             [10.0_DP,20.0_DP,30.0_DP]*arcsec_to_rad, 'attitude spans')
    call assert_true(map%array_angle_index == 0, '10D has no array variable')

    call dace_initialize(2, 11)
    call build_srp_da_parameter_map(map)
    call assert_true(map%array_angle_index == 11, '11D array index')
    call assert_close(map%array_angle_span_rad, 2.0_DP*deg_to_rad, 'array span')
    call clear_srp_scale_uncertainty()
    write(*,*) 'SRP DA parameter-map tests passed.'

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
        real(DP), intent(in) :: actual, expected
        character(len=*), intent(in) :: label
        if (abs(actual-expected) > 1.0e-15_DP) then
            write(*,*) 'FAILED: ', trim(label), actual, expected
            stop 1
        end if
    end subroutine assert_close
    subroutine assert_vector_close(actual, expected, label)
        real(DP), intent(in) :: actual(3), expected(3)
        character(len=*), intent(in) :: label
        if (maxval(abs(actual-expected)) > 1.0e-15_DP) then
            write(*,*) 'FAILED: ', trim(label), actual, expected
            stop 1
        end if
    end subroutine assert_vector_close
end program test_srp_da_parameter_map
