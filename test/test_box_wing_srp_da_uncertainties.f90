program test_box_wing_srp_da_uncertainties
    use pod_global, only: DP
    use pod_dace_classes
    use pod_spacecraft_geometry, only: spacecraft_geometry_type, srp_da_parameter_map_type
    implicit none

    type(spacecraft_geometry_type) :: geometry
    type(srp_da_parameter_map_type) :: map
    type(AlgebraicVector) :: position_da, velocity_da, acceleration_da
    real(DP) :: position(3), velocity(3), sun_position(3), earth_position(3), moon_position(3)
    real(DP) :: a0(3), dscale(3), dattitude(3,3), darray(3)
    real(DP) :: c_i_b(3,3)
    integer :: i, j, status
    character(len=256) :: message

    call dace_initialize(2, 11)
    call configure_geometry(geometry)
    position = [1000.0_DP, 2000.0_DP, 3000.0_DP]
    velocity = [0.2_DP, 0.7_DP, -0.1_DP]
    sun_position = [149597870.7_DP, 2.0e7_DP, 1.0e7_DP]
    earth_position = 0.0_DP
    moon_position = [384400.0_DP, 0.0_DP, 0.0_DP]
    map%global_scale_index = 7
    map%global_scale_span = 0.2_DP
    map%attitude_bias_index = [8,9,10]
    map%attitude_bias_span_rad = [0.01_DP,0.01_DP,0.01_DP]
    map%array_angle_index = 11
    map%array_angle_span_rad = 0.02_DP

    call position_da%init(3); call velocity_da%init(3); call acceleration_da%init(3)
    do i = 1, 3
        position_da%elements(i) = position(i) + da_var(i)
        velocity_da%elements(i) = velocity(i) + da_var(i+3)
    end do
    call geometry%compute_srp_from_ephemerides_da(position_da, velocity_da, sun_position, &
                                                   earth_position, moon_position, map, acceleration_da, &
                                                   status, message)
    call assert_true(status >= 0, '11D DA status: '//trim(message))
    a0 = acceleration_da%cons()
    do i = 1, 3
        dscale(i) = acceleration_da%elements(i)%get_deriv_value(7)
        do j = 1, 3
            dattitude(i,j) = acceleration_da%elements(i)%get_deriv_value(7+j)
        end do
        darray(i) = acceleration_da%elements(i)%get_deriv_value(11)
    end do
    call assert_close_vector(dscale, 0.2_DP*a0, 2.0e-20_DP, 'global scale derivative')
    call assert_true(maxval(abs(dattitude(:,2:3))) > 1.0e-20_DP, &
                     'pitch/yaw attitude variables affect SRP')
    call assert_true(maxval(abs(darray)) > 1.0e-20_DP, 'array-angle variable affects SRP')

    c_i_b = 0.0_DP
    c_i_b(1,1)=1.0_DP; c_i_b(2,2)=1.0_DP; c_i_b(3,3)=1.0_DP
    call geometry%compute_srp_with_real_attitude_da(position_da, sun_position, c_i_b, map, &
                                                     acceleration_da, status, message)
    call assert_true(status >= 0, 'external real C_I_B enters DA path')
    call assert_true(maxval(abs([acceleration_da%elements(1)%get_deriv_value(8), &
                                acceleration_da%elements(2)%get_deriv_value(9), &
                                acceleration_da%elements(3)%get_deriv_value(10)])) > 1.0e-20_DP, &
                     'external attitude receives DA rotation-vector biases')
    c_i_b(1,1)=2.0_DP
    call geometry%compute_srp_with_real_attitude_da(position_da, sun_position, c_i_b, map, &
                                                     acceleration_da, status, message)
    call assert_true(status < 0, 'invalid external real C_I_B rejected by DA path')

    call position_da%destroy(); call velocity_da%destroy(); call acceleration_da%destroy()
    write(*,*) 'DA box-wing uncertainty sensitivity tests passed.'

contains
    subroutine configure_geometry(g)
        type(spacecraft_geometry_type), intent(out) :: g
        g%mass_kg = 1000.0_DP
        g%box_dimensions_m = [2.0_DP,3.0_DP,4.0_DP]
        g%box_optical%absorptivity = 0.3_DP
        g%box_optical%specular_reflectivity = 0.4_DP
        g%box_optical%diffuse_reflectivity = 0.3_DP
        g%array_total_area_m2 = 20.0_DP
        g%array_tracking_mode = 'single_axis'
        g%array_hinge_axis_body = [0.0_DP,1.0_DP,0.0_DP]
        g%array_reference_normal_body = [1.0_DP,0.0_DP,0.0_DP]
        g%array_front_optical = g%box_optical
        g%array_back_optical = g%box_optical
        g%attitude_mode = 'sun'
        g%roll_reference = 'inertial_z'
        g%primary_axis_body = [1.0_DP,0.0_DP,0.0_DP]
        g%secondary_axis_body = [0.0_DP,0.0_DP,1.0_DP]
        g%pressure_1au_n_m2 = 4.56e-6_DP
    end subroutine configure_geometry
    subroutine assert_close_vector(actual, expected, atol, label)
        real(DP), intent(in) :: actual(3), expected(3), atol
        character(len=*), intent(in) :: label
        if (maxval(abs(actual-expected)) > atol) then
            write(*,*) 'FAILED: ', trim(label), actual, expected
            stop 1
        end if
    end subroutine assert_close_vector
    subroutine assert_true(condition, label)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: label
        if (.not. condition) then
            write(*,*) 'FAILED: ', trim(label)
            stop 1
        end if
    end subroutine assert_true
end program test_box_wing_srp_da_uncertainties
