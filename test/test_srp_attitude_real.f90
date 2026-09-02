program test_srp_attitude_real
    use pod_global, only: DP
    use pod_spacecraft_geometry, only: compute_pointing_attitude_real
    implicit none

    real(DP), parameter :: tol = 2.0e-13_DP
    real(DP) :: position(3), velocity(3), sun_position(3), earth_position(3), moon_position(3)
    real(DP) :: primary_body(3), secondary_body(3), c_i_b(3,3), expected(3)
    integer :: status
    character(len=256) :: message

    position = [1.0_DP, 2.0_DP, 3.0_DP]
    velocity = [0.2_DP, 0.7_DP, -0.1_DP]
    sun_position = [11.0_DP, 2.0_DP, 3.0_DP]
    earth_position = [1.0_DP, 12.0_DP, 3.0_DP]
    moon_position = [1.0_DP, 2.0_DP, 13.0_DP]
    primary_body = [0.0_DP, 0.0_DP, 1.0_DP]
    secondary_body = [0.0_DP, 1.0_DP, 0.0_DP]

    expected = [1.0_DP, 0.0_DP, 0.0_DP]
    call check_mode('sun', expected)
    expected = [0.0_DP, 1.0_DP, 0.0_DP]
    call check_mode('earth', expected)
    expected = [0.0_DP, 0.0_DP, 1.0_DP]
    call check_mode('moon', expected)

    write(*,*) 'Real SRP attitude tests passed.'

contains

    subroutine check_mode(mode, target_direction)
        character(len=*), intent(in) :: mode
        real(DP), intent(in) :: target_direction(3)
        real(DP) :: identity_error(3,3)

        call compute_pointing_attitude_real(position, velocity, sun_position, earth_position, &
                                             moon_position, mode, 'orbit_normal', primary_body, &
                                             secondary_body, c_i_b, status, message)
        call assert_true(status >= 0, trim(mode)//' attitude status: '//trim(message))
        call assert_close(matmul(c_i_b, primary_body), target_direction, trim(mode)//' primary pointing')
        identity_error = matmul(transpose(c_i_b), c_i_b)
        identity_error(1,1) = identity_error(1,1) - 1.0_DP
        identity_error(2,2) = identity_error(2,2) - 1.0_DP
        identity_error(3,3) = identity_error(3,3) - 1.0_DP
        call assert_true(maxval(abs(identity_error)) < tol, trim(mode)//' orthonormal attitude')
        call assert_true(abs(det3(c_i_b) - 1.0_DP) < tol, trim(mode)//' right-handed attitude')
    end subroutine check_mode

    real(DP) function det3(a)
        real(DP), intent(in) :: a(3,3)
        det3 = a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2)) &
             - a(1,2)*(a(2,1)*a(3,3)-a(2,3)*a(3,1)) &
             + a(1,3)*(a(2,1)*a(3,2)-a(2,2)*a(3,1))
    end function det3

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
end program test_srp_attitude_real
