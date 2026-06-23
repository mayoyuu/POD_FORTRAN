program test_orbit_error_mahalanobis
    use pod_global, only: DP
    use pod_error_analysis_module, only: compute_orbit_error
    implicit none

    real(DP) :: est_state(6), ref_state(6), cov(6,6)
    real(DP) :: pos_err(3), vel_err(3), pos_rms, vel_rms, mahal_d
    integer :: i, n_fail

    n_fail = 0
    est_state = 0.0_DP
    ref_state = 0.0_DP
    cov = 0.0_DP

    do i = 1, 6
        cov(i,i) = 4.0_DP
    end do

    est_state(1) = 2.0_DP
    call compute_orbit_error(est_state, cov, ref_state, pos_err, vel_err, pos_rms, vel_rms, mahal_d)
    call assert_close(mahal_d, 1.0_DP, 'one 1-sigma state error gives Mahalanobis distance 1', n_fail)
    call assert_close(pos_rms, 2.0_DP, 'position RMS is Euclidean norm in km', n_fail)

    est_state = 0.0_DP
    call compute_orbit_error(est_state, cov, ref_state, pos_err, vel_err, pos_rms, vel_rms, mahal_d)
    call assert_close(mahal_d, 0.0_DP, 'zero state error gives Mahalanobis distance 0', n_fail)

    if (n_fail /= 0) then
        write(*,*) 'test_orbit_error_mahalanobis failed: ', n_fail
        stop 1
    end if

    write(*,*) 'test_orbit_error_mahalanobis passed'

contains

    subroutine assert_close(actual, expected, label, n_fail)
        real(DP), intent(in) :: actual, expected
        character(len=*), intent(in) :: label
        integer, intent(inout) :: n_fail

        if (abs(actual - expected) > 1.0e-12_DP) then
            write(*,*) 'FAIL: ', trim(label)
            write(*,*) '  actual:   ', actual
            write(*,*) '  expected: ', expected
            n_fail = n_fail + 1
        end if
    end subroutine assert_close

end program test_orbit_error_mahalanobis
