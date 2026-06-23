program test_reference_orbit_lookup
    use pod_global, only: DP
    use pod_obs_io_module, only: ref_orbit_record, find_reference_by_et
    implicit none

    type(ref_orbit_record) :: refs(3)
    real(DP) :: state(6)
    logical :: found
    integer :: n_fail

    n_fail = 0
    call init_refs(refs)

    call find_reference_by_et(20.2_DP, refs, state, found, tolerance=0.5_DP)
    call assert_true(found, 'finds a reference within tolerance', n_fail)
    call assert_close(state, refs(2)%state, 'returns nearest reference state', n_fail)

    call find_reference_by_et(28.0_DP, refs, state, found, tolerance=0.5_DP)
    call assert_false(found, 'rejects a reference outside tolerance', n_fail)

    if (n_fail /= 0) then
        write(*,*) 'test_reference_orbit_lookup failed: ', n_fail
        stop 1
    end if

    write(*,*) 'test_reference_orbit_lookup passed'

contains

    subroutine init_refs(refs)
        type(ref_orbit_record), intent(out) :: refs(:)
        integer :: i

        do i = 1, size(refs)
            refs(i)%et = 10.0_DP * real(i, DP)
            refs(i)%state = real(i, DP)
        end do
    end subroutine init_refs

    subroutine assert_true(value, label, n_fail)
        logical, intent(in) :: value
        character(len=*), intent(in) :: label
        integer, intent(inout) :: n_fail

        if (.not. value) then
            write(*,*) 'FAIL: ', trim(label)
            n_fail = n_fail + 1
        end if
    end subroutine assert_true

    subroutine assert_false(value, label, n_fail)
        logical, intent(in) :: value
        character(len=*), intent(in) :: label
        integer, intent(inout) :: n_fail

        if (value) then
            write(*,*) 'FAIL: ', trim(label)
            n_fail = n_fail + 1
        end if
    end subroutine assert_false

    subroutine assert_close(actual, expected, label, n_fail)
        real(DP), intent(in) :: actual(6), expected(6)
        character(len=*), intent(in) :: label
        integer, intent(inout) :: n_fail

        if (maxval(abs(actual - expected)) > 1.0e-12_DP) then
            write(*,*) 'FAIL: ', trim(label)
            write(*,*) '  actual:   ', actual
            write(*,*) '  expected: ', expected
            n_fail = n_fail + 1
        end if
    end subroutine assert_close

end program test_reference_orbit_lookup
