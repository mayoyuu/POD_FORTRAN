program test_export_l1halo2_first_gap_case
    use iso_fortran_env, only: real64
    implicit none

    integer, parameter :: DP = real64
    character(len=*), parameter :: case_dir = &
        "OPM/supp_single_R91_floor1e-8_batch_20260625_233457/L1Halo-2"
    character(len=*), parameter :: out_dir = &
        "output/L1Halo-2_supp_single_R91_1h_floor_1p0em8_first_gap_case"
    character(len=*), parameter :: before_boundary = "2027-01-02T23:24:00.000"
    character(len=*), parameter :: after_boundary = "2027-01-13T12:44:00.384"
    character(len=*), parameter :: before_file = &
        out_dir // "/L1Halo-2_first_gap_before_errors.txt"
    character(len=*), parameter :: after_file = &
        out_dir // "/L1Halo-2_first_gap_after_errors.txt"

    integer :: before_count, after_count, n_fail

    n_fail = 0
    before_count = 0
    after_count = 0

    call export_l1halo2_first_gap_case(before_count, after_count)

    call assert_equal(before_count, 27, "before-gap error scatter row count", n_fail)
    call assert_equal(after_count, 795, "after-gap error scatter row count", n_fail)
    call assert_equal(count_rows(before_file), 27, "before-gap file row count", n_fail)
    call assert_equal(count_rows(after_file), 795, "after-gap file row count", n_fail)
    call assert_true(file_exists(before_file), "before-gap scatter file exists", n_fail)
    call assert_true(file_exists(after_file), "after-gap scatter file exists", n_fail)
    call assert_gmm_json(1, n_fail)
    call assert_gmm_json(3, n_fail)
    call assert_gmm_json(5, n_fail)

    if (n_fail /= 0) then
        write(*,*) "test_export_l1halo2_first_gap_case failed: ", n_fail
        stop 1
    end if

    write(*,*) "Exported scatter files:"
    write(*,*) "  ", before_file
    write(*,*) "  ", after_file
    write(*,*) "Exported OPM JSON files in: ", out_dir
    write(*,*) "test_export_l1halo2_first_gap_case passed"

contains

    subroutine export_l1halo2_first_gap_case(before_count, after_count)
        integer, intent(out) :: before_count, after_count
        integer :: unit_before, unit_after, exitstat, i
        integer, parameter :: n_components(3) = [1, 3, 5]

        before_count = 0
        after_count = 0

        call execute_command_line("mkdir -p " // out_dir, exitstat=exitstat)
        if (exitstat /= 0) error stop "failed to create first-gap export directory"

        open(newunit=unit_before, file=before_file, status="replace", action="write")
        open(newunit=unit_after, file=after_file, status="replace", action="write")

        do i = 1, size(n_components)
            call append_component_errors(n_components(i), unit_before, unit_after, &
                                         before_count, after_count)
            call copy_component_opm(n_components(i))
        end do

        close(unit_before)
        close(unit_after)
    end subroutine export_l1halo2_first_gap_case

    subroutine append_component_errors(n_component, unit_before, unit_after, &
                                       before_count, after_count)
        integer, intent(in) :: n_component, unit_before, unit_after
        integer, intent(inout) :: before_count, after_count
        character(len=512) :: err_path
        character(len=32) :: n_text, utc
        character(len=2048) :: line
        real(DP) :: et_seconds, state_error(6)
        integer :: unit_in, ios, k

        write(n_text, '(I0)') n_component
        err_path = case_dir // "/L1Halo-2_supp_single_R91_1h_floor_1p0em8_n" // &
                   trim(n_text) // "_p10000.err"

        open(newunit=unit_in, file=trim(err_path), status="old", action="read", iostat=ios)
        if (ios /= 0) error stop "failed to open source .err file"

        do
            read(unit_in, '(A)', iostat=ios) line
            if (ios /= 0) exit
            if (len_trim(line) == 0) cycle
            if (line(1:1) == "#") cycle

            utc = adjustl(line(1:24))
            call read_real_field(line, 25, 16, et_seconds)
            do k = 1, 6
                call read_real_field(line, 41 + (k - 1) * 14, 14, state_error(k))
            end do

            if (trim(utc) <= before_boundary) then
                call write_error_scatter(unit_before, state_error)
                before_count = before_count + 1
            else if (trim(utc) >= after_boundary) then
                call write_error_scatter(unit_after, state_error)
                after_count = after_count + 1
            end if
        end do

        close(unit_in)
    end subroutine append_component_errors


    subroutine read_real_field(line, start_col, field_width, value)
        character(len=*), intent(in) :: line
        integer, intent(in) :: start_col, field_width
        real(DP), intent(out) :: value
        character(len=64) :: field
        integer :: ios, end_col

        field = ""
        end_col = min(start_col + field_width - 1, len(line))
        if (start_col > len(line)) error stop "source .err row ended before expected field"
        field(1:end_col - start_col + 1) = line(start_col:end_col)

        read(field, *, iostat=ios) value
        if (ios /= 0) then
            write(*,*) "Failed parsing fixed-width .err field: ", trim(field)
            write(*,*) trim(line)
            error stop "failed to parse source .err fixed-width field"
        end if
    end subroutine read_real_field
    subroutine write_error_scatter(unit_out, state_error)
        integer, intent(in) :: unit_out
        real(DP), intent(in) :: state_error(6)

        write(unit_out, '(6(ES26.17E3,","))') state_error
    end subroutine write_error_scatter

    subroutine copy_component_opm(n_component)
        integer, intent(in) :: n_component
        character(len=512) :: src_path, dst_path
        character(len=32) :: n_text

        write(n_text, '(I0)') n_component
        src_path = case_dir // "/L1Halo-2_supp_single_R91_1h_floor_1p0em8_n" // &
                   trim(n_text) // "_p10000.opm.json"
        dst_path = out_dir // "/L1Halo-2_supp_single_R91_1h_floor_1p0em8_n" // &
                   trim(n_text) // "_p10000.opm.json"

        call copy_text_file(trim(src_path), trim(dst_path))
    end subroutine copy_component_opm

    subroutine copy_text_file(src_path, dst_path)
        character(len=*), intent(in) :: src_path, dst_path
        character(len=4096) :: line
        integer :: unit_src, unit_dst, ios

        open(newunit=unit_src, file=src_path, status="old", action="read", iostat=ios)
        if (ios /= 0) error stop "failed to open source JSON file"
        open(newunit=unit_dst, file=dst_path, status="replace", action="write", iostat=ios)
        if (ios /= 0) error stop "failed to open destination JSON file"

        do
            read(unit_src, '(A)', iostat=ios) line
            if (ios /= 0) exit
            write(unit_dst, '(A)') trim(line)
        end do

        close(unit_src)
        close(unit_dst)
    end subroutine copy_text_file

    subroutine assert_gmm_json(n_component, n_fail)
        integer, intent(in) :: n_component
        integer, intent(inout) :: n_fail
        character(len=512) :: json_path
        character(len=64) :: expected_count
        character(len=32) :: n_text

        write(n_text, '(I0)') n_component
        write(expected_count, '(A,I0)') '"GMM_N_COMPONENTS": ', n_component
        json_path = out_dir // "/L1Halo-2_supp_single_R91_1h_floor_1p0em8_n" // &
                    trim(n_text) // "_p10000.opm.json"

        call assert_true(file_exists(json_path), "OPM JSON exists for n" // trim(n_text), n_fail)
        call assert_file_contains(json_path, "GMM_COMPONENTS", &
                                  "OPM JSON has GMM_COMPONENTS for n" // trim(n_text), n_fail)
        call assert_file_contains(json_path, trim(expected_count), &
                                  "OPM JSON has expected component count for n" // trim(n_text), n_fail)
    end subroutine assert_gmm_json

    function file_exists(path) result(exists)
        character(len=*), intent(in) :: path
        logical :: exists

        inquire(file=trim(path), exist=exists)
    end function file_exists

    function count_rows(path) result(n_rows)
        character(len=*), intent(in) :: path
        integer :: n_rows
        character(len=2048) :: line
        integer :: unit_in, ios

        n_rows = 0
        open(newunit=unit_in, file=trim(path), status="old", action="read", iostat=ios)
        if (ios /= 0) return

        do
            read(unit_in, '(A)', iostat=ios) line
            if (ios /= 0) exit
            if (len_trim(line) > 0) n_rows = n_rows + 1
        end do

        close(unit_in)
    end function count_rows

    subroutine assert_file_contains(path, needle, label, n_fail)
        character(len=*), intent(in) :: path, needle, label
        integer, intent(inout) :: n_fail
        character(len=4096) :: line
        integer :: unit_in, ios
        logical :: found

        found = .false.
        open(newunit=unit_in, file=trim(path), status="old", action="read", iostat=ios)
        if (ios /= 0) then
            call assert_true(.false., label, n_fail)
            return
        end if

        do
            read(unit_in, '(A)', iostat=ios) line
            if (ios /= 0) exit
            if (index(line, needle) > 0) found = .true.
        end do

        close(unit_in)
        call assert_true(found, label, n_fail)
    end subroutine assert_file_contains

    subroutine assert_equal(actual, expected, label, n_fail)
        integer, intent(in) :: actual, expected
        character(len=*), intent(in) :: label
        integer, intent(inout) :: n_fail

        if (actual /= expected) then
            write(*,*) "FAIL: ", trim(label)
            write(*,*) "  actual:   ", actual
            write(*,*) "  expected: ", expected
            n_fail = n_fail + 1
        end if
    end subroutine assert_equal

    subroutine assert_true(value, label, n_fail)
        logical, intent(in) :: value
        character(len=*), intent(in) :: label
        integer, intent(inout) :: n_fail

        if (.not. value) then
            write(*,*) "FAIL: ", trim(label)
            n_fail = n_fail + 1
        end if
    end subroutine assert_true

end program test_export_l1halo2_first_gap_case