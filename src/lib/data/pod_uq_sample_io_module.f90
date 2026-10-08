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
    public :: write_uq_moments_json, read_uq_moments_json, write_ads_stats_json

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

    subroutine write_uq_moments_json(filename, state, method, status, message, epoch_et)
        character(len=*), intent(in) :: filename, method
        type(uq_state_type), intent(in) :: state
        integer, intent(out) :: status
        character(len=*), intent(out) :: message
        real(DP), optional, intent(in) :: epoch_et
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
        if (present(epoch_et)) then
            write(unit,'(A,ES24.16E3,A)') '  "epoch_et": ', epoch_et, ','
        end if
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

    !> Read a moments JSON (any of the three writer shapes) back into mean/cov.
    !> Only the common "mean" and "covariance" fields are extracted, plus the
    !> optional "epoch_et" scalar when present.
    subroutine read_uq_moments_json(filename, mean, cov, status, message, epoch_et)
        character(len=*), intent(in) :: filename
        real(DP), allocatable, intent(out) :: mean(:)
        real(DP), allocatable, intent(out) :: cov(:,:)
        integer, intent(out) :: status
        character(len=*), intent(out) :: message
        real(DP), optional, intent(out) :: epoch_et

        character(len=LINE_LENGTH) :: line
        real(DP) :: row_vals(6), tmp
        integer :: unit, io, row
        logical :: found_mean, in_cov, found_epoch

        status = 0
        message = ''
        found_mean = .false.
        found_epoch = .false.
        in_cov = .false.
        row = 0

        allocate(mean(6), cov(6,6))
        mean = 0.0_DP
        cov = 0.0_DP

        open(newunit=unit, file=trim(filename), status='old', action='read', iostat=io)
        if (io /= 0) then
            status = -1
            message = 'cannot open moments JSON'
            return
        end if

        do
            read(unit, '(A)', iostat=io) line
            if (io < 0) exit
            if (io > 0) then
                close(unit)
                status = -2
                message = 'cannot read moments JSON'
                return
            end if

            if (in_cov) then
                if (parse_bracket_6(line, row_vals)) then
                    row = row + 1
                    if (row <= 6) cov(row, :) = row_vals
                    if (row == 6) in_cov = .false.
                end if
                cycle
            end if

            if (index(line, '"mean"') > 0) then
                if (parse_bracket_6(line, row_vals)) then
                    mean = row_vals
                    found_mean = .true.
                end if
            else if (index(line, '"covariance"') > 0) then
                in_cov = .true.
                row = 0
            else if (present(epoch_et) .and. index(line, '"epoch_et"') > 0) then
                if (parse_json_scalar(line, '"epoch_et"', tmp)) then
                    epoch_et = tmp
                    found_epoch = .true.
                end if
            end if
        end do
        close(unit)

        if (.not. found_mean) then
            status = -3
            message = 'mean field not found in moments JSON'
            return
        end if
        if (row < 6) then
            status = -4
            message = 'covariance matrix incomplete in moments JSON'
            return
        end if
        if (present(epoch_et) .and. .not. found_epoch) then
            status = -5
            message = 'epoch_et field not found in moments JSON'
        end if
    end subroutine read_uq_moments_json

    !> Extract six real numbers from the [ ... ] bracket group on a line.
    function parse_bracket_6(line, vals) result(ok)
        character(len=*), intent(in) :: line
        real(DP), intent(out) :: vals(6)
        logical :: ok
        character(len=LINE_LENGTH) :: buf
        integer :: b1, b2, i, io

        ok = .false.
        b1 = index(line, '[')
        b2 = index(line, ']')
        if (b1 == 0 .or. b2 <= b1) return
        buf = line(b1+1:b2-1)
        do i = 1, len_trim(buf)
            if (buf(i:i) == ',') buf(i:i) = ' '
        end do
        read(buf, *, iostat=io) vals
        ok = (io == 0)
    end function parse_bracket_6

    !> Extract a scalar value following a "key": token on a JSON line.
    function parse_json_scalar(line, key, val) result(ok)
        character(len=*), intent(in) :: line, key
        real(DP), intent(out) :: val
        logical :: ok
        character(len=LINE_LENGTH) :: buf
        integer :: ik, ic, ie

        ok = .false.
        ik = index(line, trim(key))
        if (ik == 0) return
        ic = index(line(ik:), ':') + ik - 1
        if (ic < ik) return
        ie = index(line(ic+1:), ',')
        if (ie == 0) then
            buf = line(ic+1:)
        else
            buf = line(ic+1:ic+ie-1)
        end if
        read(buf, *, iostat=ie) val
        ok = (ie == 0)
    end function parse_json_scalar

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
