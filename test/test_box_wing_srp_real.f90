program test_box_wing_srp_real
    use pod_global, only: DP
    use pod_spacecraft_geometry, only: spacecraft_geometry_type
    implicit none

    real(DP), parameter :: au_km = 149597870.7_DP
    real(DP), parameter :: pressure = 4.0e-6_DP
    real(DP), parameter :: mass = 1000.0_DP
    real(DP), parameter :: tol = 1.0e-20_DP
    type(spacecraft_geometry_type) :: geometry
    real(DP) :: position(3), velocity(3), sun_position(3), earth_position(3), moon_position(3)
    real(DP) :: acceleration(3), expected(3)
    integer :: status
    character(len=256) :: message

    geometry%mass_kg = mass
    geometry%box_dimensions_m = [2.0_DP, 3.0_DP, 4.0_DP]
    geometry%box_optical%absorptivity = 1.0_DP
    geometry%box_optical%specular_reflectivity = 0.0_DP
    geometry%box_optical%diffuse_reflectivity = 0.0_DP
    geometry%array_total_area_m2 = 20.0_DP
    geometry%array_tracking_mode = 'single_axis'
    geometry%array_hinge_axis_body = [0.0_DP, 1.0_DP, 0.0_DP]
    geometry%array_reference_normal_body = [1.0_DP, 0.0_DP, 0.0_DP]
    geometry%array_front_optical = geometry%box_optical
    geometry%array_back_optical = geometry%box_optical
    geometry%attitude_mode = 'sun'
    geometry%roll_reference = 'inertial_z'
    geometry%primary_axis_body = [1.0_DP, 0.0_DP, 0.0_DP]
    geometry%secondary_axis_body = [0.0_DP, 0.0_DP, 1.0_DP]
    geometry%pressure_1au_n_m2 = pressure

    position = 0.0_DP
    velocity = [0.0_DP, 1.0_DP, 0.0_DP]
    sun_position = [au_km, 0.0_DP, 0.0_DP]
    earth_position = [-1000.0_DP, 0.0_DP, 0.0_DP]
    moon_position = [0.0_DP, 1000.0_DP, 0.0_DP]

    call geometry%compute_srp_from_ephemerides_real(position, velocity, sun_position, earth_position, &
                                                     moon_position, acceleration, status, message)
    ! +X box face area = Ly*Lz = 12 m^2; tracked array adds 20 m^2.
    expected = [-pressure * (12.0_DP + 20.0_DP) / mass * 1.0e-3_DP, 0.0_DP, 0.0_DP]
    call assert_true(status >= 0, 'aggregate status: '//trim(message))
    call assert_close(acceleration, expected, 'box plus double-sided array at 1 AU')

    sun_position = [2.0_DP*au_km, 0.0_DP, 0.0_DP]
    call geometry%compute_srp_from_ephemerides_real(position, velocity, sun_position, earth_position, &
                                                     moon_position, acceleration, status, message)
    call assert_close(acceleration, 0.25_DP*expected, 'inverse-square solar distance')

    write(*,*) 'Real box-wing aggregate tests passed.'

contains
    subroutine assert_close(actual, expected_value, label)
        real(DP), intent(in) :: actual(3), expected_value(3)
        character(len=*), intent(in) :: label
        if (maxval(abs(actual - expected_value)) > tol) then
            write(*,*) 'FAILED: ', trim(label), actual, expected_value
            stop 1
        end if
    end subroutine assert_close

    subroutine assert_true(condition, label)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: label
        if (.not. condition) then
            write(*,*) 'FAILED: ', trim(label)
            stop 1
        end if
    end subroutine assert_true
end program test_box_wing_srp_real
