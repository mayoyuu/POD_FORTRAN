!> @file test_compiled_da_eval_into.f90
!! @brief Verify allocation-free evaluation of a compiled DA vector.
!!
!! The initial-state SRP uncertainty study evaluates thousands of points at
!! every output epoch. The historical CompiledDA%eval function allocates a
!! result vector for every point. This test specifies a caller-owned output
!! interface so those allocations can be avoided.
program test_compiled_da_eval_into
    use pod_global, only: DP
    use pod_dace_classes, only: DA, AlgebraicVector, CompiledDA, &
        dace_initialize, active_da_count
    implicit none

    type(DA) :: x1, x2
    type(AlgebraicVector) :: state
    type(CompiledDA) :: compiled, invalid
    real(DP), allocatable :: allocated_result(:)
    real(DP) :: result(2), wrong_size(1), point(2)
    integer :: status, i, handles_before_loop

    call dace_initialize(2, 2)
    call x1%init_var(1)
    call x2%init_var(2)
    call state%init(2)
    call state%set(1, x1)
    call state%set(2, x2)
    call x1%destroy()
    call x2%destroy()

    compiled = state%compile()
    point = [0.25_DP, -0.50_DP]
    allocated_result = compiled%eval(point)

    call compiled%eval_into(point, result, status)
    call assert_true(status == 0, 'valid eval_into status')
    call assert_close(maxval(abs(result-allocated_result)), 0.0_DP, &
                      1.0e-14_DP, 'eval_into matches allocating eval')

    call compiled%eval_into(point, wrong_size, status)
    call assert_true(status == -1, 'wrong result size is rejected')

    call invalid%eval_into(point, result, status)
    call assert_true(status == -2, 'invalid compiled handle is rejected')

    handles_before_loop = active_da_count()
    do i = 1, 1000
        call compiled%eval_into(point, result, status)
        call assert_true(status == 0, 'loop evaluation status')
    end do
    call assert_true(active_da_count() == handles_before_loop, &
                     'eval_into must not create DA handles')

    call compiled%destroy()
    call state%destroy()
    if (allocated(allocated_result)) deallocate(allocated_result)

    call assert_true(active_da_count() == 0, 'all DA handles are released')
    write(*,'(a)') 'PASS: allocation-free CompiledDA evaluation.'

contains

    subroutine assert_true(condition, label)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: label
        if (.not. condition) then
            write(*,'(a)') 'FAIL: '//trim(label)
            error stop 1
        end if
    end subroutine assert_true

    subroutine assert_close(actual, expected, tolerance, label)
        real(DP), intent(in) :: actual, expected, tolerance
        character(len=*), intent(in) :: label
        if (abs(actual-expected) > tolerance) then
            write(*,'(a,3(1x,es24.16))') 'FAIL: '//trim(label), &
                actual, expected, tolerance
            error stop 1
        end if
    end subroutine assert_close

end program test_compiled_da_eval_into
