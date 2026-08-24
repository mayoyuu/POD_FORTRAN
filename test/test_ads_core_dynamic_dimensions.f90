program test_ads_core_dynamic_dimensions
    use pod_global, only: DP
    use pod_dace_classes, only: AlgebraicVector, DA, dace_initialize, &
        da_exp_sub, active_da_count, assignment(=)
    use pod_ads_split_module, only: patch_type, patch_init, patch_destroy, &
        patch_get_trunc_err, patch_get_split_dir, patch_split, sh_center, sh_width, sh_contain, &
        sh_map_point
    implicit none

    type(AlgebraicVector) :: source
    type(DA) :: x7
    type(patch_type) :: parent, left, right
    real(DP), allocatable :: center(:), width(:), values(:)
    real(DP) :: point(7), errors(3)
    integer :: n_fail, direction, active_before

    n_fail = 0
    call dace_initialize(2, 7)
    active_before = active_da_count()

    ! Use the no-temporary DA API so finalizable function temporaries cannot
    ! interfere with this ADS ownership regression test.
    call source%init(3)
    call x7%init_var(7)
    call da_exp_sub(x7, source%elements(1))
    source%elements(2) = 2.0_DP
    source%elements(3) = -3.0_DP
    call patch_init(parent, source)
    call x7%destroy()

    call assert_equal_int(parent%da_vec%size, 3, &
        'Patch preserves its three output components', n_fail)
    call patch_get_trunc_err(parent, 2, errors)
    call assert_true(errors(1) > 0.0_DP, &
        'truncation error accepts the three-component output shape', n_fail)

    direction = patch_get_split_dir(parent, 1, 2)
    call assert_equal_int(direction, 7, &
        'split direction searches all seven DA variables', n_fail)

    call patch_split(parent, 7, left, right)
    call assert_equal_int(left%da_vec%size, 3, &
        'left child preserves output component count', n_fail)
    call assert_equal_int(right%da_vec%size, 3, &
        'right child preserves output component count', n_fail)

    center = sh_center(left%history)
    width = sh_width(left%history)
    call assert_equal_int(size(center), 7, &
        'history center follows the active DACE dimension', n_fail)
    call assert_equal_int(size(width), 7, &
        'history width follows the active DACE dimension', n_fail)
    call assert_close(center(7), -0.5_DP, 1.0e-14_DP, &
        'left child center in variable seven', n_fail)
    call assert_close(width(7), 1.0_DP, 1.0e-14_DP, &
        'left child width in variable seven', n_fail)

    point = 0.0_DP
    point(7) = -0.75_DP
    call assert_true(sh_contain(left%history, point), &
        'left child contains a global unit-domain point', n_fail)
    call assert_true(.not. sh_contain(right%history, point), &
        'right child rejects a left-domain point', n_fail)
    call sh_map_point(left%history, point)
    call assert_close(point(7), -0.5_DP, 1.0e-14_DP, &
        'global point maps to the left local coordinate', n_fail)

    values = left%da_vec%eval(point)
    call assert_close(values(1), 1.0_DP - 0.75_DP + 0.5_DP * 0.75_DP**2, 1.0e-13_DP, &
        'translated left polynomial retains the global value', n_fail)

    point = 0.0_DP
    point(7) = 0.75_DP
    call assert_true(sh_contain(right%history, point), &
        'right child contains a global unit-domain point', n_fail)
    call sh_map_point(right%history, point)
    call assert_close(point(7), 0.5_DP, 1.0e-14_DP, &
        'global point maps to the right local coordinate', n_fail)
    values = right%da_vec%eval(point)
    call assert_close(values(1), 1.0_DP + 0.75_DP + 0.5_DP * 0.75_DP**2, 1.0e-13_DP, &
        'translated right polynomial retains the global value', n_fail)

    if (allocated(center)) deallocate(center)
    if (allocated(width)) deallocate(width)
    if (allocated(values)) deallocate(values)
    call patch_destroy(parent)
    call patch_destroy(left)
    call patch_destroy(right)
    call source%destroy()
    call assert_equal_int(active_da_count(), active_before, &
        'ADS split releases every temporary and Patch DA handle', n_fail)

    if (n_fail /= 0) then
        write(*,'(A,I0)') 'FAIL: ADS Core dynamic-dimension checks failed: ', n_fail
        stop 1
    end if

    write(*,'(A)') 'PASS: ADS Core separates 7 DA variables from 3 Patch components'

contains

    subroutine assert_equal_int(actual, expected, label, failures)
        integer, intent(in) :: actual, expected
        character(len=*), intent(in) :: label
        integer, intent(inout) :: failures

        if (actual /= expected) then
            write(*,'(A,A,A,I0,A,I0)') 'FAIL: ', trim(label), &
                ', got ', actual, ', expected ', expected
            failures = failures + 1
        end if
    end subroutine assert_equal_int

    subroutine assert_close(actual, expected, tolerance, label, failures)
        real(DP), intent(in) :: actual, expected, tolerance
        character(len=*), intent(in) :: label
        integer, intent(inout) :: failures

        if (abs(actual - expected) > tolerance) then
            write(*,'(A,A,A,ES22.14,A,ES22.14)') 'FAIL: ', trim(label), &
                ', got ', actual, ', expected ', expected
            failures = failures + 1
        end if
    end subroutine assert_close

    subroutine assert_true(condition, label, failures)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: label
        integer, intent(inout) :: failures

        if (.not. condition) then
            write(*,'(A,A)') 'FAIL: ', trim(label)
            failures = failures + 1
        end if
    end subroutine assert_true

end program test_ads_core_dynamic_dimensions
