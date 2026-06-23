program test_json_real_format
    use pod_global, only: DP
    use pod_data_format_module, only: format_json_real
    implicit none

    integer :: n_fail
    character(len=64) :: text

    n_fail = 0

    text = format_json_real(1.0e-129_DP)
    call assert_contains(text, 'E-', 'small negative exponent keeps E marker', n_fail)
    call assert_not_contains(text, '*', 'small exponent does not overflow field', n_fail)

    text = format_json_real(1.0e+123_DP)
    call assert_contains(text, 'E+', 'large positive exponent keeps E marker', n_fail)
    call assert_not_contains(text, '*', 'large exponent does not overflow field', n_fail)

    if (n_fail /= 0) then
        write(*,*) 'test_json_real_format failed: ', n_fail
        stop 1
    end if

    write(*,*) 'test_json_real_format passed'

contains

    subroutine assert_contains(text, needle, label, n_fail)
        character(len=*), intent(in) :: text, needle, label
        integer, intent(inout) :: n_fail

        if (index(text, needle) <= 0) then
            write(*,*) 'FAIL: ', trim(label)
            write(*,*) '  text: ', trim(text)
            write(*,*) '  expected substring: ', trim(needle)
            n_fail = n_fail + 1
        end if
    end subroutine assert_contains

    subroutine assert_not_contains(text, needle, label, n_fail)
        character(len=*), intent(in) :: text, needle, label
        integer, intent(inout) :: n_fail

        if (index(text, needle) > 0) then
            write(*,*) 'FAIL: ', trim(label)
            write(*,*) '  text: ', trim(text)
            write(*,*) '  unexpected substring: ', trim(needle)
            n_fail = n_fail + 1
        end if
    end subroutine assert_not_contains

end program test_json_real_format
