!> @file test_da_rkf78_tiny_step_no_leak.f90
!! @brief Regression test for the RKF78 tiny-final-step DA handle leak.
!!
!! Adaptive propagation can leave a roundoff-sized remainder at the requested
!! final epoch. da_rkf78_step treats dt <= 1e-15 as a no-op. The no-op guard
!! must run before allocating the 13 RKF stage vectors and nine work vectors;
!! otherwise one such call leaks (13 + 9) * 6 = 132 DA handles.
program test_da_rkf78_tiny_step_no_leak
    use pod_global, only: DP
    use pod_dace_classes, only: AlgebraicVector, dace_initialize, &
        active_da_count, assignment(=)
    use pod_da_integrator_module, only: da_rkf78_step
    implicit none

    real(DP), parameter :: TINY_DT = 1.0e-16_DP
    real(DP), parameter :: TEST_STATE(6) = [ &
        0.10_DP, -0.20_DP, 0.30_DP, 0.40_DP, -0.50_DP, 0.60_DP ]
    real(DP), parameter :: TOL = 1.0e-15_DP
    type(AlgebraicVector) :: state, state_7th, state_8th, error_estimate
    integer :: i, baseline_count, count_before, count_after

    call dace_initialize(1, 6)
    baseline_count = active_da_count()

    ! Preallocate every caller-owned vector. Its 24 handles are expected to
    ! remain live across the call and are therefore included in count_before.
    call state%init(6)
    call state_7th%init(6)
    call state_8th%init(6)
    call error_estimate%init(6)
    do i = 1, 6
        state%elements(i) = TEST_STATE(i)
    end do

    count_before = active_da_count()
    call da_rkf78_step(state, TINY_DT, 0.0_DP, state_7th, state_8th, &
                       error_estimate)
    count_after = active_da_count()

    call assert_true(count_after == count_before, &
        'tiny RKF78 step retained internal DA work handles')
    call assert_true(maxval(abs(state_7th%cons() - TEST_STATE)) <= TOL, &
        'tiny RKF78 step changed the seventh-order state')
    call assert_true(maxval(abs(state_8th%cons() - TEST_STATE)) <= TOL, &
        'tiny RKF78 step changed the eighth-order state')
    call assert_true(maxval(abs(error_estimate%cons())) <= TOL, &
        'tiny RKF78 step returned a nonzero error estimate')

    call state%destroy()
    call state_7th%destroy()
    call state_8th%destroy()
    call error_estimate%destroy()
    call assert_true(active_da_count() == baseline_count, &
        'caller-owned DA handles were not fully released')

    write(*,'(a)') 'PASS: RKF78 tiny-step branch has no DA handle leak.'

contains

    subroutine assert_true(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (.not. condition) then
            write(*,'(a)') 'FAIL: '//trim(message)
            error stop 1
        end if
    end subroutine assert_true

end program test_da_rkf78_tiny_step_no_leak
