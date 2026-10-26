!> Diagnose the first non-finite value behind the DRO-1/C3 box-wing failure.
!!
!! This program intentionally does not change the production integrator or
!! force model.  It has two modes:
!!
!!   classifier
!!     Demonstrate how the current attitude routine classifies a quiet NaN.
!!
!!   propagate OPM CONFIG POINTS KERNEL_DIR [SECONDS] [NPROBES] [ONLY_RANK]
!!     Re-run the nominal state and the same evenly-spaced Real validation
!!     probes selected by run_srp_unified. ONLY_RANK=0 runs nominal plus every
!!     selected probe; ONLY_RANK=-1 runs nominal only; a positive value runs
!!     only that selected probe rank.
!!
!! Compile propagate mode with GNU Fortran IEEE traps.  The last BEGIN record
!! printed before the backtrace identifies the failing input probe.
program test_box_wing_nan_origin
    use, intrinsic :: iso_fortran_env, only: int64, output_unit, error_unit
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_quiet_nan, ieee_value
    use pod_global, only: DP
    use pod_config, only: config
    use pod_engine_module, only: pod_engine_init
    use pod_spice, only: set_spice_kernel_path
    use pod_data_format_module, only: load_initial_opm
    use pod_integrator_module, only: adaptive_step_integrate, METHOD_RKF78
    use pod_force_model_module, only: set_propagation_epoch, set_srp_scale_error, &
        clear_srp_scale_error
    use pod_spacecraft_geometry, only: compute_pointing_attitude_real
    implicit none

    integer, parameter :: PATH_LEN = 4096
    character(len=PATH_LEN) :: mode

    call get_command_argument(1, mode)
    select case (trim(mode))
    case ('classifier')
        call run_classifier_test()
    case ('propagate')
        call run_propagation_test()
    case default
        call print_usage()
        error stop 2
    end select

contains

    subroutine run_classifier_test()
        real(DP) :: position(3), velocity(3), sun_position(3)
        real(DP) :: earth_position(3), moon_position(3), c_i_b(3,3)
        real(DP) :: primary_axis(3), secondary_axis(3)
        integer :: status
        character(len=256) :: message

        position = [0.0_DP, 0.0_DP, 0.0_DP]
        velocity = [0.0_DP, 1.0_DP, 0.0_DP]
        sun_position = [149597870.7_DP, 0.0_DP, 0.0_DP]
        earth_position = [1000.0_DP, 0.0_DP, 0.0_DP]
        moon_position = [0.0_DP, 384400.0_DP, 0.0_DP]
        primary_axis = [0.0_DP, 0.0_DP, 1.0_DP]
        secondary_axis = [0.0_DP, 1.0_DP, 0.0_DP]

        call compute_pointing_attitude_real(position, velocity, sun_position, &
            earth_position, moon_position, 'earth', 'inertial_z', primary_axis, &
            secondary_axis, c_i_b, status, message)
        call require(status >= 0, 'finite Earth-pointing control failed: '//trim(message))
        call require(all(ieee_is_finite(c_i_b)), &
            'finite Earth-pointing control returned a non-finite matrix')

        position(1) = ieee_value(0.0_DP, ieee_quiet_nan)
        call compute_pointing_attitude_real(position, velocity, sun_position, &
            earth_position, moon_position, 'earth', 'inertial_z', primary_axis, &
            secondary_axis, c_i_b, status, message)

        call require(status < 0, 'NaN attitude input was not rejected')
        write(output_unit,'(a,i0)') 'DIAG classifier_status=', status
        write(output_unit,'(a,a)') 'DIAG classifier_message=', trim(message)
        if (status == -2 .and. index(message, 'coincides with attitude target') > 0) then
            write(output_unit,'(a)') &
                'PASS: current code maps a NaN target direction to the target-coincidence error.'
        else
            write(output_unit,'(a)') &
                'PASS: NaN input was rejected with a different diagnostic; inspect the message above.'
        end if
    end subroutine run_classifier_test

    subroutine run_propagation_test()
        character(len=PATH_LEN) :: opm_path, config_path, points_path, kernel_dir, text
        real(DP) :: epoch0, nominal(6), covariance(6,6), duration_seconds
        real(DP), allocatable :: samples(:,:), eta(:)
        integer(int64), allocatable :: sample_ids(:)
        integer, allocatable :: selected(:)
        integer :: count, requested_probes, only_rank, rank, source_index, io

        if (command_argument_count() < 5) then
            call print_usage()
            error stop 2
        end if
        call get_command_argument(2, opm_path)
        call get_command_argument(3, config_path)
        call get_command_argument(4, points_path)
        call get_command_argument(5, kernel_dir)

        duration_seconds = 3600.0_DP
        requested_probes = 64
        only_rank = 0
        if (command_argument_count() >= 6) then
            call get_command_argument(6, text)
            read(text, *, iostat=io) duration_seconds
            call require(io == 0 .and. ieee_is_finite(duration_seconds) .and. &
                duration_seconds > 0.0_DP, 'invalid diagnostic duration')
        end if
        if (command_argument_count() >= 7) then
            call get_command_argument(7, text)
            read(text, *, iostat=io) requested_probes
            call require(io == 0 .and. requested_probes > 0, 'invalid probe count')
        end if
        if (command_argument_count() >= 8) then
            call get_command_argument(8, text)
            read(text, *, iostat=io) only_rank
            call require(io == 0, 'invalid ONLY_RANK')
        end if

        call set_spice_kernel_path(with_trailing_slash(trim(kernel_dir)))
        call pod_engine_init(trim(config_path))
        call require(config%use_srp, 'resolved configuration has SRP disabled')
        call require(trim(config%srp_model) == 'box_wing', &
            'resolved configuration is not box_wing')
        call require(trim(config%srp_attitude_mode) == 'earth', &
            'resolved configuration is not Earth-pointing')

        call load_initial_opm(trim(opm_path), epoch0, nominal, covariance)
        call require(ieee_is_finite(epoch0) .and. all(ieee_is_finite(nominal)) .and. &
            all(ieee_is_finite(covariance)), 'OPM contains a non-finite value')
        call read_points_file(trim(points_path), samples, eta, sample_ids, count)
        call require(count >= requested_probes, 'fewer input samples than requested probes')
        call require(all(ieee_is_finite(samples)) .and. all(ieee_is_finite(eta)), &
            'points file contains a non-finite value')
        call require(maxval(abs(eta)) <= 1000.0_DP*epsilon(1.0_DP), &
            'C3 diagnostic requires zero eta')

        allocate(selected(requested_probes))
        call evenly_spaced(count, selected)
        call require(only_rank >= -1 .and. only_rank <= requested_probes, &
            'ONLY_RANK is outside [-1,NPROBES]')

        write(output_unit,'(a,es24.16e3)') 'DIAG epoch_tdb=', epoch0
        write(output_unit,'(a,es24.16e3)') 'DIAG duration_seconds=', duration_seconds
        write(output_unit,'(a,i0)') 'DIAG population_count=', count
        write(output_unit,'(a,i0)') 'DIAG selected_probe_count=', requested_probes
        write(output_unit,'(a,3(1x,a))') 'DIAG config=', trim(config%srp_model), &
            trim(config%srp_attitude_mode), trim(config%srp_roll_reference)
        call flush(output_unit)

        call set_propagation_epoch(epoch0)
        call clear_srp_scale_error()
        if (only_rank <= 0) then
            call propagate_one('nominal', 0, 0_int64, nominal, 0.0_DP, duration_seconds)
        end if
        if (only_rank == -1) then
            write(output_unit,'(a)') 'PASS: nominal diagnostic propagation completed.'
            return
        end if

        do rank = 1, requested_probes
            if (only_rank > 0 .and. rank /= only_rank) cycle
            source_index = selected(rank)
            call propagate_one('probe', rank, sample_ids(source_index), &
                samples(:,source_index), eta(source_index), duration_seconds)
        end do
        write(output_unit,'(a)') 'PASS: all requested diagnostic propagations reached the endpoint.'
    end subroutine run_propagation_test

    subroutine propagate_one(kind, rank, sample_id, physical_state, eta_value, duration_seconds)
        character(len=*), intent(in) :: kind
        integer, intent(in) :: rank
        integer(int64), intent(in) :: sample_id
        real(DP), intent(in) :: physical_state(6), eta_value, duration_seconds
        real(DP) :: nondim(6), expected_end, endpoint_tolerance
        real(DP), allocatable :: times(:), integrated(:,:)
        integer :: n_steps, row, component

        call require(all(ieee_is_finite(physical_state)) .and. ieee_is_finite(eta_value), &
            'attempted to propagate a non-finite input state')
        nondim(1:3) = physical_state(1:3)/config%LU
        nondim(4:6) = physical_state(4:6)/config%VU
        expected_end = duration_seconds/config%TU
        endpoint_tolerance = 1000.0_DP*epsilon(1.0_DP)*max(1.0_DP, abs(expected_end))

        write(output_unit,'(a,a,a,i0,a,i0)') 'BEGIN kind=', trim(kind), &
            ' rank=', rank, ' sample_id=', sample_id
        write(output_unit,'(a,6(1x,es24.16e3))') 'BEGIN state=', physical_state
        write(output_unit,'(a,es24.16e3)') 'BEGIN eta=', eta_value
        call flush(output_unit)

        call set_srp_scale_error(eta_value)
        n_steps = 0
        call adaptive_step_integrate(nondim, 0.0_DP, expected_end, METHOD_RKF78, &
            times, integrated, n_steps, config%max_propagation_steps, &
            config%rkf78_rel_tol, config%rkf78_abs_tol, &
            config%rkf78_min_step, config%rkf78_max_step)
        call clear_srp_scale_error()

        call require(n_steps >= 1, 'integrator returned no state for '//trim(kind))
        call require(allocated(times) .and. allocated(integrated), &
            'integrator returned unallocated arrays for '//trim(kind))
        call require(n_steps <= size(times) .and. n_steps <= size(integrated,1), &
            'integrator returned an invalid step count for '//trim(kind))

        do row = 1, n_steps
            do component = 1, 6
                if (.not. ieee_is_finite(integrated(row,component))) then
                    write(error_unit,'(a,a,a,i0,a,i0,a,i0)') 'FAIL nonfinite kind=', &
                        trim(kind), ' rank=', rank, ' row=', row, ' component=', component
                    write(error_unit,'(a,es24.16e3)') 'FAIL time_nondim=', times(row)
                    error stop 1
                end if
            end do
        end do
        if (.not. ieee_is_finite(times(n_steps)) .or. &
            abs(times(n_steps)-expected_end) > endpoint_tolerance) then
            write(error_unit,'(a,a,a,i0)') 'FAIL incomplete kind=', trim(kind), ' rank=', rank
            write(error_unit,'(a,2(1x,es24.16e3))') &
                'FAIL actual_expected_time_nondim=', times(n_steps), expected_end
            error stop 1
        end if

        write(output_unit,'(a,a,a,i0,a,i0,a,es24.16e3)') 'PASS kind=', trim(kind), &
            ' rank=', rank, ' steps=', n_steps, ' end_seconds=', times(n_steps)*config%TU
        call flush(output_unit)
        deallocate(times, integrated)
    end subroutine propagate_one

    subroutine read_points_file(path, samples, eta, sample_ids, count)
        character(len=*), intent(in) :: path
        real(DP), allocatable, intent(out) :: samples(:,:), eta(:)
        integer(int64), allocatable, intent(out) :: sample_ids(:)
        integer, intent(out) :: count
        character(len=8192) :: line
        character(len=*), parameter :: expected_header = &
            'sample_id,x,y,z,vx,vy,vz,eta,s,u,wm,wc'
        real(DP) :: values(11)
        integer(int64) :: sample_id
        integer :: unit, io, row

        open(newunit=unit, file=path, status='old', action='read', iostat=io)
        call require(io == 0, 'cannot open points file: '//trim(path))
        read(unit,'(a)',iostat=io) line
        call require(io == 0 .and. trim(line) == expected_header, &
            'unexpected points header')
        count = 0
        do
            read(unit,'(a)',iostat=io) line
            if (io < 0) exit
            call require(io == 0 .and. len_trim(line) > 0, &
                'cannot count points records')
            count = count + 1
        end do
        call require(count > 0, 'points file contains no records')

        allocate(samples(6,count), eta(count), sample_ids(count))
        rewind(unit)
        read(unit,'(a)',iostat=io) line
        do row = 1, count
            read(unit,'(a)',iostat=io) line
            call require(io == 0, 'cannot read points record')
            read(line,*,iostat=io) sample_id, values
            call require(io == 0, 'invalid points record at row')
            call require(all(ieee_is_finite(values)), 'non-finite points record')
            sample_ids(row) = sample_id
            samples(:,row) = values(1:6)
            eta(row) = values(7)
        end do
        close(unit)
    end subroutine read_points_file

    subroutine evenly_spaced(count, ids)
        integer, intent(in) :: count
        integer, intent(out) :: ids(:)
        integer :: index_id, selected_count

        selected_count = size(ids)
        if (selected_count == 1) then
            ids(1) = 1
        else
            do index_id = 1, selected_count
                ids(index_id) = 1 + nint(real(index_id-1,DP)*real(count-1,DP) / &
                    real(selected_count-1,DP))
            end do
        end if
    end subroutine evenly_spaced

    function with_trailing_slash(path) result(value)
        character(len=*), intent(in) :: path
        character(len=PATH_LEN) :: value
        integer :: length

        value = trim(path)
        length = len_trim(value)
        if (length > 0 .and. value(length:length) /= '/') value = trim(value)//'/'
    end function with_trailing_slash

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (.not. condition) then
            write(error_unit,'(a)') 'FAIL: '//trim(message)
            error stop 1
        end if
    end subroutine require

    subroutine print_usage()
        write(error_unit,'(a)') 'Usage:'
        write(error_unit,'(a)') '  test_box_wing_nan_origin classifier'
        write(error_unit,'(a)') &
            '  test_box_wing_nan_origin propagate OPM CONFIG POINTS KERNEL_DIR [SECONDS] [NPROBES] [ONLY_RANK]'
    end subroutine print_usage

end program test_box_wing_nan_origin
