!> Dynamic 6D/7D uncertainty-cloud and ADS diagnostic I/O.
module pod_uq_sample_io_module
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use pod_global, only: DP
    use pod_uq_state_module, only: uq_state_type
    use pod_uq_hfem_ads_module, only: hfem_ads_stats_type
    use pod_uq_ads_coordinates_module, only: ADS_COORD_COMPONENT, ADS_COORD_WHITENED
    implicit none
    private

    integer, parameter :: LINE_LENGTH = 4096
    public :: read_uq_perturbations_csv, write_uq_particles_csv
    public :: write_uq_moments_json, write_ads_stats_json

contains

    subroutine read_uq_perturbations_csv(filename, samples, status, message)
        character(len=*), intent(in) :: filename
        real(DP), allocatable, intent(out) :: samples(:,:)
        integer, intent(out) :: status
        character(len=*), intent(out) :: message
        character(len=LINE_LENGTH) :: line
        real(DP), allocatable :: values(:)
        integer :: unit, io, parse_status, n_columns, n_rows, row
        logical :: header_skipped

        status = 0
        message = ''
        n_columns = 0
        n_rows = 0
        header_skipped = .false.
        open(newunit=unit, file=trim(filename), status='old', action='read', iostat=io)
        if (io /= 0) then
            call fail(-1, 'cannot open perturbation CSV')
            return
        end if
        do
            read(unit,'(A)',iostat=io) line
            if (io < 0) exit
            if (io > 0) then
                close(unit)
                call fail(-2, 'cannot read perturbation CSV')
                return
            end if
            call parse_numeric_line(line, values, parse_status)
            if (parse_status == 1) cycle
            if (parse_status < 0) then
                if (n_rows == 0 .and. .not. header_skipped) then
                    header_skipped = .true.
                    cycle
                end if
                close(unit)
                call fail(-3, 'non-numeric row in perturbation CSV')
                return
            end if
            if (n_columns == 0) then
                n_columns = size(values)
                if (n_columns /= 6 .and. n_columns /= 7) then
                    close(unit)
                    call fail(-4, 'perturbation CSV must contain six or seven columns')
                    return
                end if
            else if (size(values) /= n_columns) then
                close(unit)
                call fail(-5, 'mixed column counts in perturbation CSV')
                return
            end if
            if (.not. all(ieee_is_finite(values))) then
                close(unit)
                call fail(-6, 'perturbation CSV contains a non-finite value')
                return
            end if
            n_rows = n_rows + 1
        end do
        close(unit)
        if (n_rows == 0) then
            call fail(-7, 'perturbation CSV contains no samples')
            return
        end if

        allocate(samples(n_columns,n_rows))
        open(newunit=unit, file=trim(filename), status='old', action='read', iostat=io)
        row = 0
        do
            read(unit,'(A)',iostat=io) line
            if (io /= 0) exit
            call parse_numeric_line(line, values, parse_status)
            if (parse_status /= 0) cycle
            row = row + 1
            samples(:,row) = values
        end do
        close(unit)

    contains
        subroutine fail(code, text)
            integer, intent(in) :: code
            character(len=*), intent(in) :: text
            status = code
            message = text
        end subroutine fail
    end subroutine read_uq_perturbations_csv

    subroutine write_uq_particles_csv(filename, samples, status, message)
        character(len=*), intent(in) :: filename
        real(DP), intent(in) :: samples(:,:)
        integer, intent(out) :: status
        character(len=*), intent(out) :: message
        integer :: unit, io, j

        status = 0
        message = ''
        if (size(samples,1) /= 6 .and. size(samples,1) /= 7) then
            status = -1
            message = 'particle output supports only six or seven rows'
            return
        end if
        open(newunit=unit, file=trim(filename), status='replace', action='write', iostat=io)
        if (io /= 0) then
            status = -2
            message = 'cannot open particle CSV for writing'
            return
        end if
        if (size(samples,1) == 6) then
            write(unit,'(A)',iostat=io) 'x,y,z,vx,vy,vz'
        else
            write(unit,'(A)',iostat=io) 'x,y,z,vx,vy,vz,eta_srp'
        end if
        do j = 1, size(samples,2)
            write(unit,'(*(ES24.16E3,:,","))',iostat=io) samples(:,j)
            if (io /= 0) exit
        end do
        close(unit)
        if (io /= 0) then
            status = -3
            message = 'failed while writing particle CSV'
        end if
    end subroutine write_uq_particles_csv

    subroutine write_uq_moments_json(filename, state, method, status, message)
        character(len=*), intent(in) :: filename, method
        type(uq_state_type), intent(in) :: state
        integer, intent(out) :: status
        character(len=*), intent(out) :: message
        real(DP), allocatable :: skewness(:), kurtosis(:)
        integer :: unit, io, particle_count

        status = 0
        message = ''
        if (.not. allocated(state%mean) .or. .not. allocated(state%cov)) then
            status = -1
            message = 'moments are not allocated'
            return
        end if
        particle_count = 0
        if (allocated(state%samples)) then
            particle_count = size(state%samples,2)
            call state%compute_higher_moments(skewness,kurtosis)
        end if
        open(newunit=unit, file=trim(filename), status='replace', action='write', iostat=io)
        if (io /= 0) then
            status = -2
            message = 'cannot open moments JSON for writing'
            return
        end if
        write(unit,'(A)') '{'
        write(unit,'(A,A,A)') '  "method": "', trim(method), '",'
        write(unit,'(A,I0,A)') '  "particle_count": ', particle_count, ','
        write(unit,'(A)') '  "covariance_convention": "sample_n_minus_1",'
        call write_real_vector(unit, 'mean', state%mean, .true.)
        call write_real_matrix(unit, 'covariance', state%cov, allocated(skewness))
        if (allocated(skewness)) then
            call write_real_vector(unit, 'marginal_skewness', skewness, .true.)
            call write_real_vector(unit, 'marginal_kurtosis', kurtosis, .false.)
        end if
        write(unit,'(A)') '}'
        close(unit)
    end subroutine write_uq_moments_json

    subroutine write_ads_stats_json(filename, stats, status, message)
        character(len=*), intent(in) :: filename
        type(hfem_ads_stats_type), intent(in) :: stats
        integer, intent(out) :: status
        character(len=*), intent(out) :: message
        integer :: unit, io
        character(len=16) :: coordinate_name

        status = 0
        message = ''
        select case (stats%coordinate_mode)
        case (ADS_COORD_COMPONENT)
            coordinate_name = 'component'
        case (ADS_COORD_WHITENED)
            coordinate_name = 'whitened'
        case default
            coordinate_name = 'unknown'
        end select
        open(newunit=unit, file=trim(filename), status='replace', action='write', iostat=io)
        if (io /= 0) then
            status = -1
            message = 'cannot open ADS stats JSON for writing'
            return
        end if
        write(unit,'(A)') '{'
        write(unit,'(A,I0,A)') '  "requested_count": ', stats%requested_count, ','
        write(unit,'(A,I0,A)') '  "sampled_count": ', stats%sampled_count, ','
        write(unit,'(A,I0,A)') '  "rejected_count": ', stats%rejected_count, ','
        write(unit,'(A,I0,A)') '  "input_count": ', stats%input_count, ','
        write(unit,'(A,I0,A)') '  "inside_count": ', stats%inside_count, ','
        write(unit,'(A,I0,A)') '  "outside_count": ', stats%outside_count, ','
        write(unit,'(A,I0,A)') '  "propagated_count": ', stats%propagated_count, ','
        write(unit,'(A,I0,A)') '  "written_count": ', stats%written_count, ','
        write(unit,'(A,I0,A)') '  "patch_count": ', stats%n_patches, ','
        write(unit,'(A,I0,A)') '  "bfs_iterations": ', stats%bfs_iterations, ','
        write(unit,'(A,I0,A)') '  "max_queue_size": ', stats%max_queue_size, ','
        write(unit,'(A,I0,A)') '  "depth_limited_patches": ', &
            stats%depth_limited_patches, ','
        write(unit,'(A,A,A)') '  "coordinate_mode": "', trim(coordinate_name), '",'
        write(unit,'(A,I0,A)') '  "coordinate_mode_code": ', stats%coordinate_mode, ','
        write(unit,'(A,I0,A)') '  "n_variables": ', stats%n_variables, ','
        write(unit,'(A,ES24.16E3,A)') '  "domain_sigma": ', stats%domain_sigma, ','
        write(unit,'(A,ES24.16E3,A)') '  "srp_sigma": ', stats%srp_sigma, ','
        write(unit,'(A,ES24.16E3,A)') '  "elapsed_seconds": ', stats%elapsed_seconds, ','
        call write_real_matrix(unit, 'basis6', stats%basis6, .true.)
        if (allocated(stats%split_counts)) then
            call write_integer_vector(unit, 'split_counts', stats%split_counts)
        else
            write(unit,'(A)') '  "split_counts": []'
        end if
        write(unit,'(A)') '}'
        close(unit)
    end subroutine write_ads_stats_json

    subroutine parse_numeric_line(line, values, status)
        character(len=*), intent(in) :: line
        real(DP), allocatable, intent(out) :: values(:)
        integer, intent(out) :: status
        character(len=LINE_LENGTH) :: cleaned
        integer :: i, n_columns, io

        cleaned = adjustl(line)
        if (len_trim(cleaned) == 0 .or. cleaned(1:1) == '#') then
            status = 1
            return
        end if
        n_columns = 1
        do i = 1, len_trim(cleaned)
            if (cleaned(i:i) == ',') then
                cleaned(i:i) = ' '
                n_columns = n_columns + 1
            end if
        end do
        allocate(values(n_columns))
        read(cleaned,*,iostat=io) values
        status = merge(0,-1,io == 0)
    end subroutine parse_numeric_line

    subroutine write_real_vector(unit, name, values, trailing_comma)
        integer, intent(in) :: unit
        character(len=*), intent(in) :: name
        real(DP), intent(in) :: values(:)
        logical, intent(in) :: trailing_comma
        integer :: i

        write(unit,'(A)',advance='no') '  "'//trim(name)//'": ['
        do i = 1, size(values)
            if (i > 1) write(unit,'(A)',advance='no') ', '
            write(unit,'(ES24.16E3)',advance='no') values(i)
        end do
        if (trailing_comma) then
            write(unit,'(A)') '],'
        else
            write(unit,'(A)') ']'
        end if
    end subroutine write_real_vector

    subroutine write_integer_vector(unit, name, values)
        integer, intent(in) :: unit
        character(len=*), intent(in) :: name
        integer, intent(in) :: values(:)
        integer :: i

        write(unit,'(A)',advance='no') '  "'//trim(name)//'": ['
        do i = 1, size(values)
            if (i > 1) write(unit,'(A)',advance='no') ', '
            write(unit,'(I0)',advance='no') values(i)
        end do
        write(unit,'(A)') ']'
    end subroutine write_integer_vector

    subroutine write_real_matrix(unit, name, values, trailing_comma)
        integer, intent(in) :: unit
        character(len=*), intent(in) :: name
        real(DP), intent(in) :: values(:,:)
        logical, intent(in) :: trailing_comma
        integer :: i, j

        write(unit,'(A)') '  "'//trim(name)//'": ['
        do i = 1, size(values,1)
            write(unit,'(A)',advance='no') '    ['
            do j = 1, size(values,2)
                if (j > 1) write(unit,'(A)',advance='no') ', '
                write(unit,'(ES24.16E3)',advance='no') values(i,j)
            end do
            if (i < size(values,1)) then
                write(unit,'(A)') '],'
            else
                write(unit,'(A)') ']'
            end if
        end do
        if (trailing_comma) then
            write(unit,'(A)') '  ],'
        else
            write(unit,'(A)') '  ]'
        end if
    end subroutine write_real_matrix

end module pod_uq_sample_io_module
