program test_srp_array_real
    use pod_global, only: DP
    use pod_spacecraft_geometry, only: compute_array_normal_real
    implicit none

    real(DP), parameter :: tol = 2.0e-13_DP
    real(DP) :: c_i_b(3,3), e_s_i(3), hinge_b(3), reference_b(3), normal_b(3), expected(3)
    integer :: status
    character(len=256) :: message

    c_i_b = 0.0_DP
    c_i_b(1,1) = 1.0_DP
    c_i_b(2,2) = 1.0_DP
    c_i_b(3,3) = 1.0_DP
    hinge_b = [0.0_DP, 1.0_DP, 0.0_DP]
    reference_b = [1.0_DP, 0.0_DP, 0.0_DP]
    e_s_i = [1.0_DP, 1.0_DP, 1.0_DP] / sqrt(3.0_DP)
    expected = [1.0_DP, 0.0_DP, 1.0_DP] / sqrt(2.0_DP)

    call compute_array_normal_real(c_i_b, e_s_i, 'single_axis', hinge_b, reference_b, &
                                   0.0_DP, normal_b, status, message)
    call assert_true(status == 0, 'single-axis status: '//trim(message))
    call assert_close(normal_b, expected, 'single-axis projected Sun direction')
    call assert_true(abs(dot_product(normal_b, hinge_b)) < tol, 'normal perpendicular to hinge')

    call compute_array_normal_real(c_i_b, e_s_i, 'fixed', hinge_b, reference_b, &
                                   0.5_DP*acos(-1.0_DP), normal_b, status, message)
    call assert_true(status == 0, 'fixed status: '//trim(message))
    call assert_close(normal_b, [0.0_DP, 0.0_DP, -1.0_DP], 'fixed panel angle bias')

    write(*,*) 'Real SRP array tests passed.'

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
end program test_srp_array_real
