program test_uq_sample_io
    use pod_global, only: DP
    use pod_uq_state_module, only: uq_state_type
    use pod_uq_hfem_ads_module, only: hfem_ads_stats_type
    use pod_uq_sample_io_module, only: read_uq_perturbations_csv, &
        write_uq_particles_csv, write_uq_moments_json, write_ads_stats_json
    implicit none

    character(len=*), parameter :: FILE6 = 'test_uq_samples_6.csv'
    character(len=*), parameter :: FILE7 = 'test_uq_samples_7.csv'
    character(len=*), parameter :: MIXED = 'test_uq_samples_mixed.csv'
    character(len=*), parameter :: PARTICLES = 'test_uq_particles.csv'
    character(len=*), parameter :: MOMENTS = 'test_uq_moments.json'
    character(len=*), parameter :: STATS = 'test_uq_ads_stats.json'
    type(uq_state_type) :: state
    type(hfem_ads_stats_type) :: ads_stats
    real(DP), allocatable :: samples(:,:)
    character(len=256) :: message, header
    integer :: unit, status, failures

    failures = 0
    open(newunit=unit,file=FILE6,status='replace',action='write')
    write(unit,'(A)') 'x,y,z,vx,vy,vz'
    write(unit,'(A)') '1,2,3,4,5,6'
    write(unit,'(A)') '7,8,9,10,11,12'
    close(unit)
    call read_uq_perturbations_csv(FILE6,samples,status,message)
    call assert_true(status == 0 .and. size(samples,1) == 6 .and. &
        size(samples,2) == 2, 'six-column input', failures)

    open(newunit=unit,file=FILE7,status='replace',action='write')
    write(unit,'(A)') 'x,y,z,vx,vy,vz,eta_srp'
    write(unit,'(A)') '1,2,3,4,5,6,0.1'
    write(unit,'(A)') '7,8,9,10,11,12,-0.1'
    close(unit)
    call read_uq_perturbations_csv(FILE7,samples,status,message)
    call assert_true(status == 0 .and. size(samples,1) == 7, &
        'seven-column input', failures)

    open(newunit=unit,file=MIXED,status='replace',action='write')
    write(unit,'(A)') '1,2,3,4,5,6'
    write(unit,'(A)') '1,2,3,4,5,6,7'
    close(unit)
    call read_uq_perturbations_csv(MIXED,samples,status,message)
    call assert_true(status < 0, 'mixed columns are rejected', failures)

    call state%allocate_memory(7,2)
    state%samples(:,1) = [1.0_DP,2.0_DP,3.0_DP,4.0_DP,5.0_DP,6.0_DP,-0.1_DP]
    state%samples(:,2) = [2.0_DP,3.0_DP,4.0_DP,5.0_DP,6.0_DP,7.0_DP,0.1_DP]
    call state%compute_moments()
    allocate(ads_stats%split_counts(7),source=0)
    ads_stats%n_variables = 7
    ads_stats%written_count = 2
    call write_uq_particles_csv(PARTICLES,state%samples,status,message)
    call assert_true(status == 0, 'dynamic particle writer', failures)
    open(newunit=unit,file=PARTICLES,status='old',action='read')
    read(unit,'(A)') header
    close(unit)
    call assert_true(trim(header) == 'x,y,z,vx,vy,vz,eta_srp', &
        'seven-dimensional CSV header', failures)
    call write_uq_moments_json(MOMENTS,state,'ADS',status,message)
    call assert_true(status == 0, 'dynamic moments JSON writer', failures)
    call write_ads_stats_json(STATS,ads_stats,status,message)
    call assert_true(status == 0, 'ADS diagnostics JSON writer', failures)

    call state%deallocate_memory()
    call delete_file(FILE6)
    call delete_file(FILE7)
    call delete_file(MIXED)
    call delete_file(PARTICLES)
    call delete_file(MOMENTS)
    call delete_file(STATS)
    if (failures /= 0) then
        write(*,'(A,I0)') 'FAIL: dynamic UQ sample I/O checks: ', failures
        stop 1
    end if
    write(*,'(A)') 'PASS: dynamic UQ sample and ADS diagnostic I/O'

contains

    subroutine delete_file(filename)
        character(len=*), intent(in) :: filename
        integer :: delete_unit, io
        open(newunit=delete_unit,file=filename,status='old',iostat=io)
        if (io == 0) close(delete_unit,status='delete')
    end subroutine delete_file

    subroutine assert_true(condition,label,n_fail)
        logical,intent(in) :: condition
        character(len=*),intent(in) :: label
        integer,intent(inout) :: n_fail
        if (.not. condition) then
            write(*,'(A,A)') 'FAIL: ',trim(label)
            n_fail=n_fail+1
        end if
    end subroutine assert_true

end program test_uq_sample_io
