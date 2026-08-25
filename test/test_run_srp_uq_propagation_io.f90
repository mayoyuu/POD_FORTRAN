program test_run_srp_uq_propagation_io
    use pod_global, only: DP
    implicit none

    character(len=*), parameter :: OPM_FILE = 'OPM/L1Halo-1/L1Halo-1_init.opm.json'
    character(len=*), parameter :: PREFIX = '/tmp/pod_run_srp_uq_io_test'
    character(len=*), parameter :: CSV_FILE = PREFIX // '_particles.csv'
    character(len=*), parameter :: JSON_FILE = PREFIX // '_moments.json'
    character(len=*), parameter :: SUCCESS_LOG = PREFIX // '_success.log'
    integer, parameter :: N_PARTICLES = 4
    integer :: n_fail

    n_fail = 0

    call cleanup_outputs()
    call test_successful_run(n_fail)
    call test_invalid_particle_count_fails(n_fail)
    call test_dt_and_epoch_are_mutually_exclusive(n_fail)

    if (n_fail /= 0) then
        write(*,*) 'test_run_srp_uq_propagation_io failed: ', n_fail
        write(*,*) 'See log: ', SUCCESS_LOG
        stop 1
    end if

    write(*,*) 'test_run_srp_uq_propagation_io passed'

contains

    subroutine cleanup_outputs()
        integer :: exit_status, cmd_status

        call execute_command_line('rm -f ' // CSV_FILE // ' ' // JSON_FILE // ' ' // &
                                  SUCCESS_LOG // ' ' // PREFIX // '_bad_n.log ' // &
                                  PREFIX // '_dt_et.log', wait=.true., &
                                  exitstat=exit_status, cmdstat=cmd_status)
    end subroutine cleanup_outputs

    subroutine test_successful_run(n_fail)
        integer, intent(inout) :: n_fail
        character(len=2048) :: command
        integer :: exit_status, cmd_status

        command = 'fpm run run_srp_uq_propagation -- ' // &
                  '-opm ' // OPM_FILE // ' ' // &
                  '-dt 60 ' // &
                  '-o ' // PREFIX // ' ' // &
                  '-n 4 ' // &
                  '-da 2 ' // &
                  '-cr 1.33 ' // &
                  '-smr 8.1e-3 ' // &
                  '-rp 4.2e-6 ' // &
                  '-srp-mean 0.125 ' // &
                  '--zero-init-cov ' // &
                  '-srp-sigma 0.05 ' // &
                  '> ' // SUCCESS_LOG // ' 2>&1'

        call execute_command_line(trim(command), wait=.true., exitstat=exit_status, cmdstat=cmd_status)
        call assert_command_success(cmd_status, exit_status, 'valid run_srp_uq_propagation invocation', n_fail)

        call assert_file_exists(CSV_FILE, 'particles CSV output', n_fail)
        call assert_file_exists(JSON_FILE, 'moments JSON output', n_fail)
        call assert_particles_csv(CSV_FILE, n_fail)
        call assert_json_output(JSON_FILE, n_fail)
        call assert_file_contains(SUCCESS_LOG, 'Initial orbit covariance: ZERO', &
                                  'zero initial covariance log', n_fail)
        call assert_file_contains(JSON_FILE, '5.000000000000000E-002', &
                                  'nonzero eta sigma retained', n_fail)
    end subroutine test_successful_run

    subroutine test_invalid_particle_count_fails(n_fail)
        integer, intent(inout) :: n_fail
        character(len=2048) :: command
        integer :: exit_status, cmd_status

        command = 'fpm run run_srp_uq_propagation -- ' // &
                  '-opm ' // OPM_FILE // ' -dt 60 -o ' // PREFIX // '_bad_n -n 1 ' // &
                  '> ' // PREFIX // '_bad_n.log 2>&1'

        call execute_command_line(trim(command), wait=.true., exitstat=exit_status, cmdstat=cmd_status)
        call assert_command_failure(cmd_status, exit_status, '-n 1 validation', n_fail)
    end subroutine test_invalid_particle_count_fails

    subroutine test_dt_and_epoch_are_mutually_exclusive(n_fail)
        integer, intent(inout) :: n_fail
        character(len=2048) :: command
        integer :: exit_status, cmd_status

        command = 'fpm run run_srp_uq_propagation -- ' // &
                  '-opm ' // OPM_FILE // ' -dt 60 -et 2027-01-01T00:10:00 -o ' // &
                  PREFIX // '_dt_et -n 4 ' // &
                  '> ' // PREFIX // '_dt_et.log 2>&1'

        call execute_command_line(trim(command), wait=.true., exitstat=exit_status, cmdstat=cmd_status)
        call assert_command_failure(cmd_status, exit_status, '-dt and -et mutual exclusion', n_fail)
    end subroutine test_dt_and_epoch_are_mutually_exclusive

    subroutine assert_command_success(cmd_status, exit_status, label, n_fail)
        integer, intent(in) :: cmd_status, exit_status
        character(len=*), intent(in) :: label
        integer, intent(inout) :: n_fail

        if (cmd_status /= 0 .or. exit_status /= 0) then
            write(*,*) 'FAIL: ', trim(label)
            write(*,*) '  cmd_status=', cmd_status, ' exit_status=', exit_status
            n_fail = n_fail + 1
        end if
    end subroutine assert_command_success

    subroutine assert_command_failure(cmd_status, exit_status, label, n_fail)
        integer, intent(in) :: cmd_status, exit_status
        character(len=*), intent(in) :: label
        integer, intent(inout) :: n_fail

        if (cmd_status /= 0) then
            write(*,*) 'FAIL: ', trim(label)
            write(*,*) '  command processor failed, cmd_status=', cmd_status
            n_fail = n_fail + 1
        else if (exit_status == 0) then
            write(*,*) 'FAIL: ', trim(label)
            write(*,*) '  expected non-zero exit status'
            n_fail = n_fail + 1
        end if
    end subroutine assert_command_failure

    subroutine assert_file_exists(path, label, n_fail)
        character(len=*), intent(in) :: path, label
        integer, intent(inout) :: n_fail
        logical :: exists

        inquire(file=path, exist=exists)
        if (.not. exists) then
            write(*,*) 'FAIL: missing ', trim(label), ': ', trim(path)
            n_fail = n_fail + 1
        end if
    end subroutine assert_file_exists

    subroutine assert_particles_csv(path, n_fail)
        character(len=*), intent(in) :: path
        integer, intent(inout) :: n_fail
        character(len=1024) :: line
        real(DP) :: values(7), eta_min, eta_max
        integer :: unit, io_status, read_status, row_count

        open(newunit=unit, file=path, status='old', action='read', iostat=io_status)
        if (io_status /= 0) then
            write(*,*) 'FAIL: could not open CSV: ', trim(path)
            n_fail = n_fail + 1
            return
        end if

        read(unit, '(A)', iostat=io_status) line
        if (io_status /= 0 .or. trim(line) /= 'x,y,z,vx,vy,vz,eta_srp') then
            write(*,*) 'FAIL: unexpected CSV header: ', trim(line)
            n_fail = n_fail + 1
        end if

        row_count = 0
        eta_min = huge(1.0_DP)
        eta_max = -huge(1.0_DP)
        do
            read(unit, '(A)', iostat=io_status) line
            if (io_status /= 0) exit
            if (len_trim(line) == 0) cycle

            values = 0.0_DP
            read(line, *, iostat=read_status) values
            if (read_status /= 0) then
                write(*,*) 'FAIL: CSV row is not seven real values: ', trim(line)
                n_fail = n_fail + 1
            else
                eta_min = min(eta_min, values(7))
                eta_max = max(eta_max, values(7))
            end if
            row_count = row_count + 1
        end do

        close(unit)

        if (row_count /= N_PARTICLES) then
            write(*,*) 'FAIL: CSV row count mismatch, got ', row_count, ' expected ', N_PARTICLES
            n_fail = n_fail + 1
        end if

        if (eta_max <= eta_min) then
            write(*,*) 'FAIL: nonzero eta sigma did not produce sample variation'
            n_fail = n_fail + 1
        end if
    end subroutine assert_particles_csv

    subroutine assert_json_output(path, n_fail)
        character(len=*), intent(in) :: path
        integer, intent(inout) :: n_fail

        call assert_file_contains(path, '"method": "DA-MC-SRP"', 'JSON method', n_fail)
        call assert_file_contains(path, '"n_particles": 4', 'JSON n_particles', n_fail)
        call assert_file_contains(path, '"da_order": 2', 'JSON da_order', n_fail)
        call assert_file_contains(path, '"state_labels"', 'JSON state_labels key', n_fail)
        call assert_file_contains(path, 'eta_srp', 'JSON eta_srp label', n_fail)
        call assert_file_contains(path, '"srp"', 'JSON SRP block', n_fail)
        call assert_file_contains(path, '"eta_mean_input"', 'JSON eta mean input', n_fail)
        call assert_file_contains(path, '"eta_sigma_input"', 'JSON eta sigma input', n_fail)
        call assert_file_contains(path, '"mean"', 'JSON mean vector', n_fail)
        call assert_file_contains(path, '"covariance"', 'JSON covariance matrix', n_fail)
    end subroutine assert_json_output

    subroutine assert_file_contains(path, needle, label, n_fail)
        character(len=*), intent(in) :: path, needle, label
        integer, intent(inout) :: n_fail
        character(len=2048) :: line
        integer :: unit, io_status
        logical :: found

        open(newunit=unit, file=path, status='old', action='read', iostat=io_status)
        if (io_status /= 0) then
            write(*,*) 'FAIL: could not open file for substring check: ', trim(path)
            n_fail = n_fail + 1
            return
        end if

        found = .false.
        do
            read(unit, '(A)', iostat=io_status) line
            if (io_status /= 0) exit
            if (index(line, needle) > 0) then
                found = .true.
                exit
            end if
        end do
        close(unit)

        if (.not. found) then
            write(*,*) 'FAIL: missing ', trim(label)
            write(*,*) '  expected substring: ', trim(needle)
            n_fail = n_fail + 1
        end if
    end subroutine assert_file_contains

end program test_run_srp_uq_propagation_io
