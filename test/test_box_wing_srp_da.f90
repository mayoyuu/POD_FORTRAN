program test_box_wing_srp_da
    use pod_global, only: DP
    use pod_dace_classes
    use pod_spacecraft_geometry, only: spacecraft_geometry_type, srp_da_parameter_map_type
    implicit none

    integer, parameter :: nvar = 6
    real(DP), parameter :: h = 100.0_DP
    type(spacecraft_geometry_type) :: geometry
    type(srp_da_parameter_map_type) :: da_map
    type(AlgebraicVector) :: position_da, velocity_da, acceleration_da
    real(DP) :: position(3), velocity(3), sun_position(3), earth_position(3), moon_position(3)
    real(DP) :: acceleration_real(3), acceleration_da_constant(3)
    real(DP) :: jacobian_da(3,3), jacobian_fd(3,3), plus(3), minus(3), ap(3), am(3)
    integer :: i, j, status
    character(len=256) :: message

    call dace_initialize(2, nvar)
    call configure_geometry(geometry)
    position = [1000.0_DP, 2000.0_DP, 3000.0_DP]
    velocity = [0.2_DP, 0.7_DP, -0.1_DP]
    sun_position = [149597870.7_DP, 2.0e7_DP, 1.0e7_DP]
    earth_position = 0.0_DP
    moon_position = [384400.0_DP, 0.0_DP, 0.0_DP]

    call position_da%init(3)
    call velocity_da%init(3)
    call acceleration_da%init(3)
    do i = 1, 3
        position_da%elements(i) = position(i) + da_var(i)
        velocity_da%elements(i) = velocity(i) + da_var(i+3)
    end do

    call geometry%compute_srp_from_ephemerides_real(position, velocity, sun_position, earth_position, &
                                                     moon_position, acceleration_real, status, message)
    call assert_true(status >= 0, 'real status: '//trim(message))
    call geometry%compute_srp_from_ephemerides_da(position_da, velocity_da, sun_position, earth_position, &
                                                   moon_position, da_map, acceleration_da, status, message)
    call assert_true(status >= 0, 'DA status: '//trim(message))
    acceleration_da_constant = acceleration_da%cons()
    call assert_vector_close(acceleration_da_constant, acceleration_real, 2.0e-20_DP, &
                             'DA constant equals real')

    do i = 1, 3
        do j = 1, 3
            jacobian_da(i,j) = acceleration_da%elements(i)%get_deriv_value(j)
        end do
    end do
    do j = 1, 3
        plus = position
        minus = position
        plus(j) = plus(j) + h
        minus(j) = minus(j) - h
        call geometry%compute_srp_from_ephemerides_real(plus, velocity, sun_position, earth_position, &
                                                        moon_position, ap, status, message)
        call geometry%compute_srp_from_ephemerides_real(minus, velocity, sun_position, earth_position, &
                                                        moon_position, am, status, message)
        jacobian_fd(:,j) = (ap - am) / (2.0_DP*h)
    end do
    call assert_matrix_close(jacobian_da, jacobian_fd, 2.0e-22_DP, 2.0e-5_DP, &
                             'DA position Jacobian equals finite difference')

    call position_da%destroy()
    call velocity_da%destroy()
    call acceleration_da%destroy()
    write(*,*) 'DA box-wing SRP tests passed.'

contains
    subroutine configure_geometry(g)
        type(spacecraft_geometry_type), intent(out) :: g
        g%mass_kg = 1000.0_DP
        g%box_dimensions_m = [2.0_DP, 3.0_DP, 4.0_DP]
        g%box_optical%absorptivity = 0.3_DP
        g%box_optical%specular_reflectivity = 0.4_DP
        g%box_optical%diffuse_reflectivity = 0.3_DP
        g%array_total_area_m2 = 0.0_DP
        g%attitude_mode = 'sun'
        g%roll_reference = 'inertial_z'
        g%primary_axis_body = [1.0_DP, 0.0_DP, 0.0_DP]
        g%secondary_axis_body = [0.0_DP, 0.0_DP, 1.0_DP]
        g%pressure_1au_n_m2 = 4.56e-6_DP
    end subroutine configure_geometry

    subroutine assert_vector_close(actual, expected, atol, label)
        real(DP), intent(in) :: actual(3), expected(3), atol
        character(len=*), intent(in) :: label
        if (maxval(abs(actual-expected)) > atol) then
            write(*,*) 'FAILED: ', trim(label), actual, expected
            stop 1
        end if
    end subroutine assert_vector_close

    subroutine assert_matrix_close(actual, expected, atol, rtol, label)
        real(DP), intent(in) :: actual(3,3), expected(3,3), atol, rtol
        character(len=*), intent(in) :: label
        real(DP) :: limit
        limit = atol + rtol*maxval(abs(expected))
        if (maxval(abs(actual-expected)) > limit) then
            write(*,*) 'FAILED: ', trim(label), maxval(abs(actual-expected)), limit
            stop 1
        end if
    end subroutine assert_matrix_close

    subroutine assert_true(condition, label)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: label
        if (.not. condition) then
            write(*,*) 'FAILED: ', trim(label)
            stop 1
        end if
    end subroutine assert_true
end program test_box_wing_srp_da
