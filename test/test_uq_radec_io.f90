program test_uq_radec_io
    use pod_global, only: DP
    use pod_uq_state_module, only: uq_state_type
    use pod_uq_sample_io_module, only: write_uq_moments_json, read_uq_moments_json
    implicit none

    character(len=*), parameter :: MOM6 = 'test_uq_radec_mom6.json'
    character(len=*), parameter :: MOM7 = 'test_uq_radec_mom7.json'
    character(len=*), parameter :: MOM_NOEPOCH = 'test_uq_radec_noepoch.json'
    type(uq_state_type) :: state
    real(DP), allocatable :: mean(:), cov(:,:)
    real(DP) :: epoch_et, epoch_ref
    character(len=256) :: message
    integer :: status, failures

    failures = 0

    ! --- Test 1: 6D round-trip with epoch_et ---
    call state%allocate_memory(6, 3)
    state%samples(:,1) = [1.0_DP, 2.0_DP, 3.0_DP, 0.1_DP, 0.2_DP, 0.3_DP]
    state%samples(:,2) = [1.5_DP, 2.5_DP, 3.5_DP, 0.2_DP, 0.3_DP, 0.4_DP]
    state%samples(:,3) = [0.5_DP, 1.5_DP, 2.5_DP, 0.0_DP, 0.1_DP, 0.2_DP]
    call state%compute_moments()
    epoch_ref = 123456789.0_DP
    call write_uq_moments_json(MOM6, state, 'ADS', status, message, epoch_et=epoch_ref)
    call assert_true(status == 0, 'write 6D moments', failures)

    epoch_et = 0.0_DP
    call read_uq_moments_json(MOM6, mean, cov, status, message, epoch_et=epoch_et)
    call assert_true(status == 0, 'read 6D moments', failures)
    call assert_true(allocated(mean) .and. size(mean) == 6, 'mean size 6', failures)
    call assert_true(allocated(cov) .and. size(cov,1) == 6 .and. size(cov,2) == 6, &
        'cov size 6x6', failures)
    if (allocated(mean) .and. allocated(state%mean)) then
        call assert_true(maxval(abs(mean - state%mean)) < 1.0e-10_DP, 'mean 6D round-trip', failures)
    end if
    if (allocated(cov) .and. allocated(state%cov)) then
        call assert_true(maxval(abs(cov - state%cov)) < 1.0e-10_DP, 'cov 6D round-trip', failures)
    end if
    call assert_true(abs(epoch_et - epoch_ref) < 1.0e-6_DP, 'epoch round-trip', failures)
    call state%deallocate_memory()
    if (allocated(mean)) deallocate(mean)
    if (allocated(cov)) deallocate(cov)

    ! --- Test 2: 7D (ADS with eta_srp) -> reader extracts 6D position/velocity block ---
    call state%allocate_memory(7, 3)
    state%samples(:,1) = [1.0_DP, 2.0_DP, 3.0_DP, 0.1_DP, 0.2_DP, 0.3_DP, 0.10_DP]
    state%samples(:,2) = [1.5_DP, 2.5_DP, 3.5_DP, 0.2_DP, 0.3_DP, 0.4_DP, -0.05_DP]
    state%samples(:,3) = [0.5_DP, 1.5_DP, 2.5_DP, 0.0_DP, 0.1_DP, 0.2_DP, 0.02_DP]
    call state%compute_moments()
    call write_uq_moments_json(MOM7, state, 'ADS', status, message, epoch_et=epoch_ref)
    call assert_true(status == 0, 'write 7D moments', failures)

    epoch_et = 0.0_DP
    call read_uq_moments_json(MOM7, mean, cov, status, message, epoch_et=epoch_et)
    call assert_true(status == 0, 'read 7D moments', failures)
    if (allocated(mean) .and. allocated(state%mean)) then
        call assert_true(maxval(abs(mean - state%mean(1:6))) < 1.0e-10_DP, &
            'mean 7D->6D truncation', failures)
    end if
    if (allocated(cov) .and. allocated(state%cov)) then
        call assert_true(maxval(abs(cov - state%cov(1:6,1:6))) < 1.0e-10_DP, &
            'cov 7D->6D truncation', failures)
    end if
    call state%deallocate_memory()
    if (allocated(mean)) deallocate(mean)
    if (allocated(cov)) deallocate(cov)

    ! --- Test 3: missing epoch_et -> status -5 ---
    call state%allocate_memory(6, 2)
    state%samples(:,1) = [1.0_DP, 2.0_DP, 3.0_DP, 0.1_DP, 0.2_DP, 0.3_DP]
    state%samples(:,2) = [1.5_DP, 2.5_DP, 3.5_DP, 0.2_DP, 0.3_DP, 0.4_DP]
    call state%compute_moments()
    call write_uq_moments_json(MOM_NOEPOCH, state, 'ADS', status, message)
    call assert_true(status == 0, 'write moments without epoch', failures)

    epoch_et = 0.0_DP
    call read_uq_moments_json(MOM_NOEPOCH, mean, cov, status, message, epoch_et=epoch_et)
    call assert_true(status == -5, 'missing epoch_et yields -5', failures)
    call state%deallocate_memory()
    if (allocated(mean)) deallocate(mean)
    if (allocated(cov)) deallocate(cov)

    call delete_file(MOM6)
    call delete_file(MOM7)
    call delete_file(MOM_NOEPOCH)

    if (failures /= 0) then
        write(*,'(A,I0)') 'FAIL: uq radec I/O checks: ', failures
        stop 1
    end if
    write(*,'(A)') 'PASS: uq radec moments JSON round-trip'

contains

    subroutine delete_file(filename)
        character(len=*), intent(in) :: filename
        integer :: delete_unit, io
        open(newunit=delete_unit, file=filename, status='old', iostat=io)
        if (io == 0) close(delete_unit, status='delete')
    end subroutine delete_file

    subroutine assert_true(condition, label, n_fail)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: label
        integer, intent(inout) :: n_fail
        if (.not. condition) then
            write(*,'(A,A)') 'FAIL: ', trim(label)
            n_fail = n_fail + 1
        end if
    end subroutine assert_true

end program test_uq_radec_io
