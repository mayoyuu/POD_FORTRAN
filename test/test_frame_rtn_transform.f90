!> @file test_frame_rtn_transform.f90
!! @brief Unit tests for reusable inertial-to-RTN vector and covariance transforms.
!!
!! The RTN convention tested here is:
!!   R = r/|r|, N = (r x v)/|r x v|, T = N x R
!! and rows of C_RTN_FROM_I are R^T, T^T, N^T, so u_RTN=C*u_I.
program test_frame_rtn_transform
    use pod_global, only: DP
    use pod_frame_module, only: build_rtn_rotation, transform_vector_to_rtn, &
                                transform_covariance3_to_rtn, transform_covariance6_to_rtn
    implicit none

    real(DP), parameter :: TOL = 5.0e-13_DP
    real(DP) :: r(3), v(3), c(3,3), identity(3,3), expected(3,3)
    real(DP) :: vector_i(3), vector_rtn(3), p3(3,3), p3_rtn(3,3)
    real(DP) :: p6(6,6), p6_rtn(6,6), b6(6,6)
    real(DP) :: rhat(3), that(3), nhat(3)
    integer :: status, i
    character(len=256) :: message

    identity = 0.0_DP
    do i = 1, 3
        identity(i,i) = 1.0_DP
    end do

    ! For this circular equatorial state, RTN and inertial axes coincide.
    r = [7000.0_DP, 0.0_DP, 0.0_DP]
    v = [0.0_DP, 7.5_DP, 0.0_DP]
    call build_rtn_rotation(r, v, c, status, message)
    call assert_true(status == 0, 'axis-aligned RTN status: '//trim(message))
    call assert_matrix_close(c, identity, TOL, 'axis-aligned RTN matrix')

    vector_i = [1.0_DP, 2.0_DP, 3.0_DP]
    call transform_vector_to_rtn(vector_i, c, vector_rtn)
    call assert_vector_close(vector_rtn, vector_i, TOL, 'identity vector transform')

    ! A general state must produce a proper, right-handed orthonormal matrix.
    r = [7000.0_DP, 1000.0_DP, 2000.0_DP]
    v = [-1.0_DP, 7.0_DP, 2.0_DP]
    call build_rtn_rotation(r, v, c, status, message)
    call assert_true(status == 0, 'general RTN status: '//trim(message))
    call assert_matrix_close(matmul(c, transpose(c)), identity, 2.0e-12_DP, &
                             'RTN orthonormality')
    call assert_close(determinant3(c), 1.0_DP, 2.0e-12_DP, 'RTN determinant')

    rhat = c(1,:)
    that = c(2,:)
    nhat = c(3,:)
    call assert_vector_close(cross3(nhat, rhat), that, 2.0e-12_DP, &
                             'T equals N cross R')

    ! The explicit 3x3 covariance helper must match C*P*C^T.
    p3 = reshape([4.0_DP, 0.5_DP, -0.2_DP, &
                  0.5_DP, 2.0_DP,  0.3_DP, &
                 -0.2_DP, 0.3_DP,  1.0_DP], [3,3])
    call transform_covariance3_to_rtn(p3, c, p3_rtn)
    expected = matmul(c, matmul(p3, transpose(c)))
    call assert_matrix_close(p3_rtn, expected, 2.0e-12_DP, '3x3 covariance transform')
    call assert_close(trace3(p3_rtn), trace3(p3), 2.0e-12_DP, &
                      '3x3 covariance trace invariance')

    ! The 6x6 helper rotates position, velocity and cross covariance blocks.
    p6 = 0.0_DP
    p6(1:3,1:3) = p3
    p6(4:6,4:6) = 0.25_DP*p3
    p6(1:3,4:6) = 0.1_DP*p3
    p6(4:6,1:3) = transpose(p6(1:3,4:6))
    b6 = 0.0_DP
    b6(1:3,1:3) = c
    b6(4:6,4:6) = c
    call transform_covariance6_to_rtn(p6, c, p6_rtn)
    call assert_matrix6_close(p6_rtn, matmul(b6,matmul(p6,transpose(b6))), &
                              3.0e-12_DP, '6x6 covariance transform')
    call assert_close(sum([(p6_rtn(i,i),i=1,6)]), sum([(p6(i,i),i=1,6)]), &
                      3.0e-12_DP, '6x6 covariance trace invariance')

    ! RTN is undefined for a zero radius or zero angular momentum.
    call build_rtn_rotation([0.0_DP,0.0_DP,0.0_DP], v, c, status, message)
    call assert_true(status < 0, 'zero position must be rejected')
    call build_rtn_rotation([1.0_DP,0.0_DP,0.0_DP], &
                            [2.0_DP,0.0_DP,0.0_DP], c, status, message)
    call assert_true(status < 0, 'radial-only velocity must be rejected')

    write(*,'(a)') 'PASS: inertial/RTN vector and covariance transforms.'

contains

    pure function cross3(a, b) result(cross)
        real(DP), intent(in) :: a(3), b(3)
        real(DP) :: cross(3)
        cross = [a(2)*b(3)-a(3)*b(2), a(3)*b(1)-a(1)*b(3), &
                 a(1)*b(2)-a(2)*b(1)]
    end function cross3

    pure real(DP) function determinant3(a)
        real(DP), intent(in) :: a(3,3)
        determinant3 = a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2)) - &
                       a(1,2)*(a(2,1)*a(3,3)-a(2,3)*a(3,1)) + &
                       a(1,3)*(a(2,1)*a(3,2)-a(2,2)*a(3,1))
    end function determinant3

    pure real(DP) function trace3(a)
        real(DP), intent(in) :: a(3,3)
        trace3 = a(1,1)+a(2,2)+a(3,3)
    end function trace3

    subroutine assert_true(condition, label)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: label
        if (.not. condition) then
            write(*,'(a)') 'FAIL: '//trim(label)
            error stop 1
        end if
    end subroutine assert_true

    subroutine assert_close(actual, expected_value, tolerance, label)
        real(DP), intent(in) :: actual, expected_value, tolerance
        character(len=*), intent(in) :: label
        if (abs(actual-expected_value) > tolerance) then
            write(*,'(a,3(1x,es24.16))') 'FAIL: '//trim(label), actual, expected_value, tolerance
            error stop 1
        end if
    end subroutine assert_close

    subroutine assert_vector_close(actual, expected_value, tolerance, label)
        real(DP), intent(in) :: actual(3), expected_value(3), tolerance
        character(len=*), intent(in) :: label
        if (maxval(abs(actual-expected_value)) > tolerance) then
            write(*,'(a,1x,es24.16)') 'FAIL: '//trim(label), maxval(abs(actual-expected_value))
            error stop 1
        end if
    end subroutine assert_vector_close

    subroutine assert_matrix_close(actual, expected_value, tolerance, label)
        real(DP), intent(in) :: actual(3,3), expected_value(3,3), tolerance
        character(len=*), intent(in) :: label
        if (maxval(abs(actual-expected_value)) > tolerance) then
            write(*,'(a,1x,es24.16)') 'FAIL: '//trim(label), maxval(abs(actual-expected_value))
            error stop 1
        end if
    end subroutine assert_matrix_close

    subroutine assert_matrix6_close(actual, expected_value, tolerance, label)
        real(DP), intent(in) :: actual(6,6), expected_value(6,6), tolerance
        character(len=*), intent(in) :: label
        if (maxval(abs(actual-expected_value)) > tolerance) then
            write(*,'(a,1x,es24.16)') 'FAIL: '//trim(label), maxval(abs(actual-expected_value))
            error stop 1
        end if
    end subroutine assert_matrix6_close

end program test_frame_rtn_transform

