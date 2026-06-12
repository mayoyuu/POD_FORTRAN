!> @file run_radec_convert.f90
!> @brief 将 HFEM UQ 粒子传播结果转换到站心 RADEC 坐标系
!> @author Song Yu
!> @date 2026-06-12
!>
!> Usage:
!>   fpm run run_radec_convert -- -i <particles.csv> -s <station_id> -et <epoch> [-o <output>] [--site <site.json>] [--deg]
!>
!> Input CSV format:  x,y,z,vx,vy,vz  (J2000, same as run_HFEM_uprop output)
!> Output CSV format: ra,dec  (radians by default, degrees with --deg)
program run_radec_convert
    use pod_global, only: DP, MAX_STRING_LEN
    use pod_engine_module, only: pod_engine_init
    use pod_spice, only: str2et
    use pod_measurement_model_module, only: set_station_from_geodetic, compute_measurement
    use pod_obs_io_module, only: station_record, preload_stations, find_station_by_id
    use pod_measurement_base_module, only: observation_station
    use pod_basicmath_module, only: PI

    implicit none

    character(len=MAX_STRING_LEN) :: input_file, output_file, station_id, epoch_str
    character(len=MAX_STRING_LEN) :: site_json, config_file
    real(DP) :: et
    logical :: use_degrees, has_input, has_station, has_epoch
    integer :: num_args, i, u_in, u_out, ios, line_count
    character(len=MAX_STRING_LEN) :: arg_str, line
    type(station_record), allocatable :: station_list(:)
    type(observation_station) :: station
    real(DP) :: pos_j2000(3), state(6), measurement(2)
    real(DP) :: ra, dec

    ! Defaults
    site_json    = 'config/site-used.json'
    config_file  = 'config/config.txt'
    use_degrees  = .false.
    has_input    = .false.
    has_station  = .false.
    has_epoch    = .false.

    ! Parse CLI arguments
    num_args = command_argument_count()
    i = 1
    do while (i <= num_args)
        call get_command_argument(i, arg_str)
        select case (trim(arg_str))

            case ('-i')
                call get_command_argument(i+1, arg_str)
                input_file = trim(arg_str)
                has_input = .true.
                i = i + 1

            case ('-s')
                call get_command_argument(i+1, arg_str)
                station_id = trim(arg_str)
                has_station = .true.
                i = i + 1

            case ('-et')
                call get_command_argument(i+1, arg_str)
                epoch_str = trim(arg_str)
                has_epoch = .true.
                i = i + 1

            case ('-o', '--output')
                call get_command_argument(i+1, arg_str)
                output_file = trim(arg_str)
                i = i + 1

            case ('--site')
                call get_command_argument(i+1, arg_str)
                site_json = trim(arg_str)
                i = i + 1

            case ('--deg')
                use_degrees = .true.

            case ('-cfg', '--config')
                call get_command_argument(i+1, arg_str)
                config_file = trim(arg_str)
                i = i + 1

            case default
                write(*,*) 'Warning: ignoring unknown argument: ', trim(arg_str)
        end select
        i = i + 1
    end do

    ! Validate
    if (.not. has_input) then
        write(*,*) 'Error: -i <particles.csv> is required.'
        stop 1
    end if
    if (.not. has_station) then
        write(*,*) 'Error: -s <station_id> is required.'
        stop 1
    end if
    if (.not. has_epoch) then
        write(*,*) 'Error: -et <epoch> is required.'
        stop 1
    end if

    ! Auto-generate output filename
    if (output_file == '') then
        i = index(input_file, '.', back=.true.)
        if (i > 0) then
            output_file = input_file(1:i-1) // '_radec.csv'
        else
            output_file = trim(input_file) // '_radec.csv'
        end if
    end if

    ! Init engine + SPICE
    write(*,*) '>>> Initializing POD engine...'
    call pod_engine_init(trim(config_file))
    write(*,*) '>>> Engine initialized.'

    ! Load stations
    write(*,*) '>>> Loading stations from: ', trim(site_json)
    call preload_stations(site_json, station_list)
    station = find_station_by_id(station_id, station_list)
    write(*,*) '>>> Station ', trim(station_id), ' loaded.'

    ! Parse epoch
    call str2et(trim(epoch_str), et)

    write(*,*) '========================================'
    write(*,*) '  RADEC Conversion Configuration'
    write(*,*) '========================================'
    write(*,*) 'Input         : ', trim(input_file)
    write(*,*) 'Station       : ', trim(station_id)
    write(*,*) 'Epoch (UTC)   : ', trim(epoch_str)
    write(*,*) 'ET (TDB sec)  : ', et
    write(*,*) 'Output        : ', trim(output_file)
    write(*,*) 'Output units  : ', merge('degrees', 'radians', use_degrees)
    write(*,*) '----------------------------------------'

    ! Process CSV
    open(newunit=u_in, file=trim(input_file), status='old', action='read', iostat=ios)
    if (ios /= 0) then
        write(*,*) 'Error: cannot open input file: ', trim(input_file)
        stop 1
    end if

    open(newunit=u_out, file=trim(output_file), status='replace', action='write')
    if (use_degrees) then
        write(u_out, '(A)') 'ra_deg,dec_deg'
    else
        write(u_out, '(A)') 'ra_rad,dec_rad'
    end if

    ! Skip header line
    read(u_in, '(A)', iostat=ios) line
    if (ios /= 0) then
        write(*,*) 'Error: input file is empty.'
        stop 1
    end if

    line_count = 0
    do
        read(u_in, '(A)', iostat=ios) line
        if (ios /= 0) exit

        read(line, *, iostat=ios) state
        if (ios /= 0) cycle

        pos_j2000 = state(1:3)

        call compute_measurement(state, et, station, 'OPTICAL', measurement)

        ra  = measurement(1)
        dec = measurement(2)

        if (use_degrees) then
            ra  = ra  * 180.0_DP / PI
            dec = dec * 180.0_DP / PI
        end if

        write(u_out, '(ES22.14, ",", ES22.14)') ra, dec
        line_count = line_count + 1

        if (mod(line_count, 100000) == 0) then
            write(*,*) '  Processed ', line_count, ' particles...'
        end if
    end do

    close(u_in)
    close(u_out)

    deallocate(station_list)

    write(*,*) '========================================'
    write(*,*) '  Done. ', line_count, ' particles written to: ', trim(output_file)
    write(*,*) '========================================'

end program run_radec_convert
