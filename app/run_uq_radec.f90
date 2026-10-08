!> @file run_uq_radec.f90
!> @brief 将 HFEM UQ 传播结果(粒子 CSV / 矩 JSON)转换到站心 RA/DEC 平面
!> @author Song Yu
!> @date 2026-09-29
!>
!> Usage:
!>   fpm run run_uq_radec -- -i <particles.csv> -j <moments.json> -s <station_id> -et <epoch> \
!>       [-o <out_prefix>] [--site <site.json>] [--deg|--rad] [-cfg <config.txt>]
!>
!>   -i 粒子 CSV(可选); -j 矩 JSON(可选); 至少给一个。
!>   -s 测站 ID(必填); -et 观测历元(可选, 若 -j 的 JSON 含 epoch_et 则可省略)。
!>   默认输出角秒(arcsec); --deg 输出度, --rad 输出弧度。
program run_uq_radec
    use pod_global, only: DP, MAX_STRING_LEN
    use pod_engine_module, only: pod_engine_init
    use pod_spice, only: str2et
    use pod_measurement_model_module, only: compute_measurement
    use pod_obs_io_module, only: station_record, preload_stations, find_station_by_id
    use pod_measurement_base_module, only: observation_station
    use pod_basicmath_module, only: PI
    use pod_uq_sample_io_module, only: read_uq_moments_json
    use pod_uq_transform_ut_module, only: generate_sigma_points, reconstruct_ut_moments

    implicit none

    character(len=MAX_STRING_LEN) :: input_file, moments_file, station_id, epoch_str
    character(len=MAX_STRING_LEN) :: output_csv, output_json, out_prefix
    character(len=MAX_STRING_LEN) :: site_json, config_file
    character(len=8) :: unit_mode
    character(len=MAX_STRING_LEN) :: csv_header
    character(len=4096) :: line
    real(DP) :: et, ra, dec, angle_factor
    logical :: has_input, has_moments, has_station, has_epoch, has_prefix
    integer :: num_args, i, u_in, u_out, u_json, ios, line_count, skipped, ncol, nsig, k
    character(len=MAX_STRING_LEN) :: arg_str
    type(station_record), allocatable :: station_list(:)
    type(observation_station) :: station
    real(DP) :: state(6)
    real(DP), allocatable :: measurement(:)
    real(DP), allocatable :: mean6(:), cov6(:,:)
    real(DP), allocatable :: sigmas_x(:,:), weights(:), sigmas_y(:,:)
    real(DP), allocatable :: mean_y(:), P_yy(:,:), P_xy(:,:)
    logical :: is_angle(2)
    real(DP) :: mean_radec(2), std_radec(2), cov_radec(2,2)
    integer :: io_status
    character(len=256) :: io_message

    ! Defaults
    site_json   = 'config/site-used.json'
    config_file = 'config/config.txt'
    unit_mode   = 'arcsec'
    has_input   = .false.
    has_moments = .false.
    has_station = .false.
    has_epoch   = .false.
    has_prefix  = .false.
    out_prefix  = ''

    ! Parse CLI
    num_args = command_argument_count()
    i = 1
    do while (i <= num_args)
        call get_command_argument(i, arg_str)
        select case (trim(arg_str))
        case ('-i')
            call get_command_argument(i+1, arg_str)
            input_file = trim(arg_str); has_input = .true.; i = i + 1
        case ('-j')
            call get_command_argument(i+1, arg_str)
            moments_file = trim(arg_str); has_moments = .true.; i = i + 1
        case ('-s')
            call get_command_argument(i+1, arg_str)
            station_id = trim(arg_str); has_station = .true.; i = i + 1
        case ('-et')
            call get_command_argument(i+1, arg_str)
            epoch_str = trim(arg_str); has_epoch = .true.; i = i + 1
        case ('-o', '--output')
            call get_command_argument(i+1, arg_str)
            out_prefix = trim(arg_str); has_prefix = .true.; i = i + 1
        case ('--site')
            call get_command_argument(i+1, arg_str)
            site_json = trim(arg_str); i = i + 1
        case ('--deg')
            unit_mode = 'deg'
        case ('--rad')
            unit_mode = 'rad'
        case ('-cfg', '--config')
            call get_command_argument(i+1, arg_str)
            config_file = trim(arg_str); i = i + 1
        case default
            write(*,*) 'Warning: ignoring unknown argument: ', trim(arg_str)
        end select
        i = i + 1
    end do

    ! Validate
    if (.not. has_input .and. .not. has_moments) then
        write(*,*) 'Error: at least one of -i <particles.csv> or -j <moments.json> is required.'
        stop 1
    end if
    if (.not. has_station) then
        write(*,*) 'Error: -s <station_id> is required.'
        stop 1
    end if

    ! Angle scale factor and CSV header
    select case (trim(unit_mode))
    case ('deg')
        angle_factor = 180.0_DP / PI
        csv_header = 'ra_deg,dec_deg'
    case ('rad')
        angle_factor = 1.0_DP
        csv_header = 'ra_rad,dec_rad'
    case default
        angle_factor = 180.0_DP * 3600.0_DP / PI
        csv_header = 'ra_arcsec,dec_arcsec'
    end select

    ! Init engine + stations
    write(*,*) '>>> Initializing POD engine...'
    call pod_engine_init(trim(config_file))
    write(*,*) '>>> Loading stations from: ', trim(site_json)
    call preload_stations(site_json, station_list)
    station = find_station_by_id(station_id, station_list)

    ! Resolve epoch
    if (has_epoch) then
        call str2et(trim(epoch_str), et)
    else if (has_moments) then
        call read_uq_moments_json(trim(moments_file), mean6, cov6, io_status, io_message, epoch_et=et)
        if (io_status == -5) then
            write(*,*) 'Error: ', trim(moments_file), ' has no epoch_et; provide -et <epoch>.'
            stop 1
        else if (io_status /= 0) then
            write(*,*) 'Error: ', trim(io_message)
            stop 1
        end if
    else
        write(*,*) 'Error: -et <epoch> or -j <moments.json> is required to determine epoch.'
        stop 1
    end if

    allocate(measurement(2))

    write(*,*) '========================================'
    write(*,*) '  UQ -> RA/DEC Conversion'
    write(*,*) '========================================'
    write(*,*) 'Station       : ', trim(station_id)
    write(*,*) 'ET (TDB sec)  : ', et
    write(*,*) 'Output units  : ', trim(unit_mode)
    write(*,*) '----------------------------------------'

    ! ---------- CSV path ----------
    if (has_input) then
        if (has_prefix) then
            output_csv = trim(out_prefix) // '_radec.csv'
        else
            output_csv = trim(strip_extension(input_file)) // '_radec.csv'
        end if

        open(newunit=u_in, file=trim(input_file), status='old', action='read', iostat=ios)
        if (ios /= 0) then
            write(*,*) 'Error: cannot open input CSV: ', trim(input_file)
            stop 1
        end if
        open(newunit=u_out, file=trim(output_csv), status='replace', action='write')
        write(u_out,'(A)') trim(csv_header)

        ! Skip header line
        read(u_in,'(A)',iostat=ios) line
        if (ios /= 0) then
            write(*,*) 'Error: input CSV is empty.'
            stop 1
        end if

        line_count = 0
        skipped = 0
        do
            read(u_in,'(A)',iostat=ios) line
            if (ios /= 0) exit
            call count_columns(line, ncol)
            if (ncol /= 6 .and. ncol /= 7) then
                skipped = skipped + 1
                cycle
            end if
            read(line, *, iostat=ios) state
            if (ios /= 0) then
                skipped = skipped + 1
                cycle
            end if
            call compute_measurement(state, et, station, 'OPTICAL', measurement)
            ra  = measurement(1) * angle_factor
            dec = measurement(2) * angle_factor
            write(u_out,'(ES22.14,",",ES22.14)') ra, dec
            line_count = line_count + 1
        end do
        close(u_in); close(u_out)

        if (line_count == 0) then
            write(*,*) 'Error: no valid data rows in input CSV.'
            stop 1
        end if
        if (skipped > 0) write(*,*) '  Skipped ', skipped, ' non-6/7-column rows.'
        write(*,*) '  CSV: ', line_count, ' particles -> ', trim(output_csv)
    end if

    ! ---------- JSON path ----------
    if (has_moments) then
        if (.not. allocated(mean6)) then
            call read_uq_moments_json(trim(moments_file), mean6, cov6, io_status, io_message)
            if (io_status /= 0) then
                write(*,*) 'Error: ', trim(io_message)
                stop 1
            end if
        end if

        if (has_prefix) then
            output_json = trim(out_prefix) // '_radec_moments.json'
        else
            output_json = trim(strip_extension(moments_file)) // '_radec_moments.json'
        end if

        if (maxval(abs(cov6)) <= 0.0_DP) then
            write(*,*) '  Warning: zero covariance; using mean-only RA/DEC with zero covariance.'
            call compute_measurement(mean6, et, station, 'OPTICAL', measurement)
            mean_radec(1) = measurement(1)
            mean_radec(2) = measurement(2)
            cov_radec = 0.0_DP
        else
            call generate_sigma_points(mean6, cov6, 0.5_DP, sigmas_x, weights)
            nsig = size(weights)
            allocate(sigmas_y(2, nsig))
            do k = 1, nsig
                state = sigmas_x(:, k)
                call compute_measurement(state, et, station, 'OPTICAL', measurement)
                sigmas_y(1, k) = measurement(1)
                sigmas_y(2, k) = measurement(2)
            end do
            is_angle = [.true., .false.]
            call reconstruct_ut_moments(sigmas_x, sigmas_y, weights, mean_y, P_yy, P_xy, is_angle)
            mean_radec = mean_y
            cov_radec = P_yy
        end if

        ! Standard deviation and unit conversion (common to both branches)
        std_radec(1) = sqrt(max(cov_radec(1,1), 0.0_DP))
        std_radec(2) = sqrt(max(cov_radec(2,2), 0.0_DP))
        mean_radec = mean_radec * angle_factor
        std_radec = std_radec * angle_factor
        cov_radec = cov_radec * angle_factor**2

        open(newunit=u_json, file=trim(output_json), status='replace', action='write', iostat=ios)
        if (ios /= 0) then
            write(*,*) 'Error: cannot open output JSON: ', trim(output_json)
            stop 1
        end if
        write(u_json,'(A)') '{'
        write(u_json,'(A)') '  "method": "radec_ut",'
        write(u_json,'(A,ES24.16E3,A)') '  "epoch_et": ', et, ','
        write(u_json,'(A)',advance='no') '  "mean": ['
        write(u_json,'(ES24.16E3)',advance='no') mean_radec(1)
        write(u_json,'(A,ES24.16E3)',advance='no') ', ', mean_radec(2)
        write(u_json,'(A)') '],'
        write(u_json,'(A)',advance='no') '  "std": ['
        write(u_json,'(ES24.16E3)',advance='no') std_radec(1)
        write(u_json,'(A,ES24.16E3)',advance='no') ', ', std_radec(2)
        write(u_json,'(A)') '],'
        write(u_json,'(A)') '  "covariance": ['
        write(u_json,'(A,ES24.16E3,A,ES24.16E3,A)') '    [', cov_radec(1,1), ', ', cov_radec(1,2), '],'
        write(u_json,'(A,ES24.16E3,A,ES24.16E3,A)') '    [', cov_radec(2,1), ', ', cov_radec(2,2), ']'
        write(u_json,'(A)') '  ]'
        write(u_json,'(A)') '}'
        close(u_json)
        write(*,*) '  JSON: mean/std/cov -> ', trim(output_json)
        write(*,*) '  RA  sigma = ', std_radec(1), ' ', trim(unit_mode)
        write(*,*) '  DEC sigma = ', std_radec(2), ' ', trim(unit_mode)
    end if

    deallocate(station_list)
    write(*,*) '========================================'
    write(*,*) '  Done.'
    write(*,*) '========================================'

contains

    function strip_extension(name) result(stem)
        character(len=*), intent(in) :: name
        character(len=MAX_STRING_LEN) :: stem
        integer :: p
        p = index(name, '.', back=.true.)
        if (p > 0) then
            stem = name(1:p-1)
        else
            stem = name
        end if
    end function strip_extension

    subroutine count_columns(line, n)
        character(len=*), intent(in) :: line
        integer, intent(out) :: n
        integer :: j
        logical :: in_token
        n = 0
        in_token = .false.
        do j = 1, len_trim(line)
            if (line(j:j) == ',' .or. line(j:j) == ' ' .or. line(j:j) == achar(9)) then
                in_token = .false.
            else if (.not. in_token) then
                n = n + 1
                in_token = .true.
            end if
        end do
    end subroutine count_columns

end program run_uq_radec
