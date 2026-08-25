program test_run_hfem_uprop_zero_cov_io
    implicit none

    character(len=*), parameter :: prefix = '/tmp/test_hfem_zero_cov'
    character(len=*), parameter :: log_file = '/tmp/test_hfem_zero_cov.log'
    character(len=2048) :: command
    integer :: cmd_status, exit_status, u, ios
    character(len=512) :: line
    logical :: found

    call cleanup_outputs()

    command = 'fpm run run_HFEM_uprop -- ' // &
              '-cfg config/config.txt ' // &
              '-opm OPM/L1Halo-1/L1Halo-1_init.opm.json ' // &
              '-m DA -dt 60 -o ' // prefix // ' -n 4 -da 2 ' // &
              '--zero-init-cov > ' // log_file // ' 2>&1'
    call execute_command_line(command, wait=.true., cmdstat=cmd_status, exitstat=exit_status)

    if (cmd_status /= 0 .or. exit_status /= 0) then
        write(*,*) 'FAIL: HFEM zero-cov command failed, cmd_status=', cmd_status, &
                   ' exit_status=', exit_status
        error stop 1
    end if

    inquire(file=prefix // '_moments.json', exist=found)
    if (.not. found) error stop 'HFEM zero-cov moments output missing'

    found = .false.
    open(newunit=u, file=log_file, status='old', action='read', iostat=ios)
    if (ios /= 0) error stop 'HFEM zero-cov log missing'
    do
        read(u, '(A)', iostat=ios) line
        if (ios /= 0) exit
        if (index(line, 'Initial orbit covariance: ZERO') > 0) found = .true.
    end do
    close(u)
    if (.not. found) error stop 'HFEM zero-cov log marker missing'

    call cleanup_outputs()
    write(*,*) 'test_run_hfem_uprop_zero_cov_io passed'

contains

    subroutine cleanup_outputs()
        call execute_command_line('rm -f ' // prefix // '_particles.csv ' // &
                                  prefix // '_moments.json ' // log_file, wait=.true.)
    end subroutine cleanup_outputs

end program test_run_hfem_uprop_zero_cov_io
