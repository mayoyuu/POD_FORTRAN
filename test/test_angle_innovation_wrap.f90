program test_angle_innovation_wrap
    use pod_global, only: DP
    use pod_basicmath_module, only: PI, wrap_angle_rad
    implicit none

    call test_ra_wrap_near_zero()
    call test_no_wrap_for_small_angles()
    call test_wrap_positive_wraparound()
    call test_wrap_negative_wraparound()
    call test_identity_for_small_delta()

    write(*,*) 'test_angle_innovation_wrap passed'

contains

    subroutine assert_scalar_close(actual, expected, tol, message)
        real(DP), intent(in) :: actual, expected, tol
        character(len=*), intent(in) :: message
        if (abs(actual - expected) > tol) then
            write(*,*) 'FAIL: ', trim(message)
            write(*,*) '  actual  = ', actual
            write(*,*) '  expected= ', expected
            error stop 1
        end if
    end subroutine assert_scalar_close

    subroutine test_ra_wrap_near_zero()
        real(DP) :: ra_meas, ra_pred, innov, innov_wrapped
        real(DP), parameter :: tol = 1.0e-12_DP

        ra_meas  = 359.999_DP * PI / 180.0_DP
        ra_pred  =   0.001_DP * PI / 180.0_DP
        innov = ra_meas - ra_pred
        innov_wrapped = wrap_angle_rad(innov)
        call assert_scalar_close(innov_wrapped, -0.002_DP * PI / 180.0_DP, tol, &
            'RA 359.999deg vs 0.001deg should wrap to ~ -0.002deg')
    end subroutine test_ra_wrap_near_zero

    subroutine test_no_wrap_for_small_angles()
        real(DP) :: innov, wrapped
        real(DP), parameter :: tol = 1.0e-12_DP

        innov = 0.1_DP * PI / 180.0_DP
        wrapped = wrap_angle_rad(innov)
        call assert_scalar_close(wrapped, innov, tol, &
            'small angle delta should pass through unchanged')
    end subroutine test_no_wrap_for_small_angles

    subroutine test_wrap_positive_wraparound()
        real(DP) :: delta, wrapped
        real(DP), parameter :: tol = 1.0e-12_DP

        delta = 350.0_DP * PI / 180.0_DP
        wrapped = wrap_angle_rad(delta)
        call assert_scalar_close(wrapped, -10.0_DP * PI / 180.0_DP, tol, &
            '+350 deg should wrap to -10 deg')
    end subroutine test_wrap_positive_wraparound

    subroutine test_wrap_negative_wraparound()
        real(DP) :: delta, wrapped
        real(DP), parameter :: tol = 1.0e-12_DP

        delta = -350.0_DP * PI / 180.0_DP
        wrapped = wrap_angle_rad(delta)
        call assert_scalar_close(wrapped, 10.0_DP * PI / 180.0_DP, tol, &
            '-350 deg should wrap to +10 deg')
    end subroutine test_wrap_negative_wraparound

    subroutine test_identity_for_small_delta()
        real(DP) :: delta, wrapped
        real(DP), parameter :: tol = 1.0e-15_DP
        integer :: i

        do i = -179, 179
            delta = real(i, DP) * PI / 180.0_DP
            wrapped = wrap_angle_rad(delta)
            call assert_scalar_close(wrapped, delta, tol, &
                'delta within [-179,179] deg should be identity')
        end do
    end subroutine test_identity_for_small_delta

end program test_angle_innovation_wrap
