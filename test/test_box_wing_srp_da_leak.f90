program test_box_wing_srp_da_leak
    use pod_global, only: DP
    use pod_dace_classes
    use pod_spacecraft_geometry, only: spacecraft_geometry_type, srp_da_parameter_map_type
    implicit none

    type(spacecraft_geometry_type) :: geometry
    type(srp_da_parameter_map_type) :: map
    type(AlgebraicVector) :: position_da, velocity_da, acceleration_da
    real(DP) :: position(3), velocity(3), sun_position(3), earth_position(3), moon_position(3)
    integer :: i, iteration, status, baseline_count, warm_count, after_count, final_count
    character(len=256) :: message

    call dace_initialize(2, 11)
    baseline_count = active_da_count()
    call configure_geometry(geometry)
    position = [1000.0_DP,2000.0_DP,3000.0_DP]
    velocity = [0.2_DP,0.7_DP,-0.1_DP]
    sun_position = [149597870.7_DP,2.0e7_DP,1.0e7_DP]
    earth_position = 0.0_DP
    moon_position = [384400.0_DP,0.0_DP,0.0_DP]
    map%global_scale_index = 7
    map%global_scale_span = 0.2_DP
    map%attitude_bias_index = [8,9,10]
    map%attitude_bias_span_rad = 0.01_DP
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
    call assert_true(status >= 0, 'warm-up status: '//trim(message))
    warm_count = active_da_count()

    do iteration = 1, 200
        call geometry%compute_srp_from_ephemerides_da(position_da, velocity_da, sun_position, &
                                                       earth_position, moon_position, map, acceleration_da, &
                                                       status, message)
        call assert_true(status >= 0, 'loop status: '//trim(message))
    end do
    after_count = active_da_count()
    call assert_true(after_count == warm_count, 'active DA count stable over repeated calls')

    call position_da%destroy(); call velocity_da%destroy(); call acceleration_da%destroy()
    final_count = active_da_count()
    call assert_true(final_count == baseline_count, 'all test-owned DA handles released')
    write(*,*) 'DA box-wing leak test passed.'

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
    subroutine assert_true(condition, label)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: label
        if (.not. condition) then
            write(*,*) 'FAILED: ', trim(label)
            stop 1
        end if
    end subroutine assert_true
end program test_box_wing_srp_da_leak
