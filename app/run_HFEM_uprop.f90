!> @file run_uq_propagation.f90
!> @brief 统一的不确定性传播 CLI 入口，支持 MC/DA/UT 三种方法
!> @author Song Yu
!> @date 2026-05-26
!>
!> Usage:
!>   fpm run run_HFEM_uprop -- -opm <file> -m <MC|DA|UT> -dt <seconds> -o <prefix>
!>   fpm run run_HFEM_uprop -- -opm <file> -m <MC|DA|UT> -et <epoch> -n 5000 -o <prefix>
!> fpm run run_HFEM_uprop -- -opm input/TD1_2604_2_times_100.opm -m DA -et 2026-06-12T12:00:00 -o output/TD1_2604_2_TO_260612_times_100
!> Output:
!>   MC/DA: <prefix>_particles.csv + <prefix>_moments.json (mean/cov/skewness/kurtosis)
!>   UT:    <prefix>_moments.json (mean/cov)
program run_HFEM_uprop
    use pod_global, only: DP, MAX_STRING_LEN
    use pod_engine_module, only: pod_engine_init
    use pod_dace_classes, only: dace_initialize
    use pod_spice, only: str2et
    use pod_data_format_module, only: load_initial_opm
    use pod_uq_propagation, only: run_uq_propagation, METHOD_MC, METHOD_DA, METHOD_UT
    use pod_uq_state_module, only: uq_state_type

    implicit none

    character(len=MAX_STRING_LEN) :: opm_file, method_str, output_prefix, config_file
    character(len=MAX_STRING_LEN) :: epoch_str, arg_str
    character(len=MAX_STRING_LEN) :: json_path, csv_path
    real(DP) :: dt_seconds, t_end_et, epoch0, dt
    integer  :: method_switch, n_particles, da_order, i, num_args, ext_pos
    logical  :: has_dt, has_et, has_opm, has_method, has_output

    type(uq_state_type) :: initial_state, final_state
    real(DP), allocatable :: skewness(:), kurtosis(:)
    real(DP) :: state(6), cov(6,6)
    integer :: u_csv, j

    ! Defaults
    config_file  = 'config/config.txt'
    n_particles  = 1000000
    da_order     = 4
    has_dt       = .false.
    has_et       = .false.
    has_opm      = .false.
    has_method   = .false.
    has_output   = .false.
    dt_seconds   = 0.0_DP
    t_end_et     = 0.0_DP
    epoch_str    = ''
    output_prefix = ''

    ! Parse CLI arguments
    num_args = command_argument_count()
    i = 1
    do while (i <= num_args)
        call get_command_argument(i, arg_str)
        select case (trim(arg_str))

            case ('-opm')
                call get_command_argument(i+1, arg_str)
                opm_file = trim(arg_str)
                has_opm = .true.
                i = i + 1

            case ('-m', '--method')
                call get_command_argument(i+1, arg_str)
                method_str = trim(arg_str)
                has_method = .true.
                i = i + 1

            case ('-dt')
                call get_command_argument(i+1, arg_str)
                read(arg_str, *) dt_seconds
                has_dt = .true.
                i = i + 1

            case ('-et')
                call get_command_argument(i+1, arg_str)
                epoch_str = trim(arg_str)
                has_et = .true.
                i = i + 1

            case ('-n', '--n-particles')
                call get_command_argument(i+1, arg_str)
                read(arg_str, *) n_particles
                i = i + 1

            case ('-da', '--da-order')
                call get_command_argument(i+1, arg_str)
                read(arg_str, *) da_order
                i = i + 1

            case ('-o', '--output')
                call get_command_argument(i+1, arg_str)
                output_prefix = trim(arg_str)
                has_output = .true.
                i = i + 1

            case ('-cfg', '--config')
                call get_command_argument(i+1, arg_str)
                config_file = trim(arg_str)
                i = i + 1

            case default
                write(*,*) 'Warning: ignoring unknown argument: ', trim(arg_str)
        end select
        i = i + 1
    end do

    ! Validate arguments
    if (.not. has_opm) then
        write(*,*) 'Error: -opm <file> is required.'
        stop 1
    end if
    if (.not. has_method) then
        write(*,*) 'Error: -m <MC|DA|UT> is required.'
        stop 1
    end if
    if (.not. has_dt .and. .not. has_et) then
        write(*,*) 'Error: either -dt <seconds> or -et <epoch> is required.'
        stop 1
    end if
    if (has_dt .and. has_et) then
        write(*,*) 'Error: -dt and -et are mutually exclusive. Provide only one.'
        stop 1
    end if

    ! Map method
    select case (trim(method_str))
        case ('MC', 'mc')
            method_switch = METHOD_MC
        case ('DA', 'da')
            method_switch = METHOD_DA
        case ('UT', 'ut')
            method_switch = METHOD_UT
        case default
            write(*,*) 'Error: unknown method "', trim(method_str), '". Use MC, DA, or UT.'
            stop 1
    end select

    ! Auto-generate output prefix from OPM filename if not specified
    if (.not. has_output) then
        ext_pos = index(opm_file, '.', back=.true.)
        if (ext_pos > 1) then
            output_prefix = opm_file(1:ext_pos-1) // '_' // trim(method_str)
        else
            output_prefix = trim(opm_file) // '_' // trim(method_str)
        end if
    end if

    ! Init engine
    write(*,*) '>>> Initializing POD physics engine...'
    call pod_engine_init(trim(config_file))
    write(*,*) '>>> Engine initialized.'

    ! Init DA if needed
    if (method_switch == METHOD_DA) then
        write(*,*) '>>> Initializing DACE for DA propagation...'
        call dace_initialize(da_order, 6)
        write(*,*) '>>> DACE initialized with order=', da_order, ' and nvars=6.'
    end if

    ! Load OPM
    write(*,*) '>>> Loading OPM file: ', trim(opm_file)
    call load_initial_opm(trim(opm_file), epoch0, state, cov)

    ! Compute dt
    if (has_dt) then
        dt = dt_seconds
    else
        call str2et(trim(epoch_str), t_end_et)
        dt = t_end_et - epoch0
        if (dt <= 0.0_DP) then
            write(*,*) 'Error: final epoch must be later than OPM epoch.'
            write(*,*) '  OPM epoch (TDB): ', epoch0
            write(*,*) '  Target epoch : ', t_end_et
            stop 1
        end if
    end if

    write(*,*) '========================================'
    write(*,*) '  UQ Propagation Configuration'
    write(*,*) '========================================'
    write(*,*) 'OPM file      : ', trim(opm_file)
    write(*,*) 'Method        : ', trim(method_str)
    write(*,*) 'Epoch0 (TDB)  : ', epoch0
    write(*,*) 'dt (seconds)  : ', dt
    write(*,*) 'n_particles   : ', n_particles
    if (method_switch == METHOD_DA) then
        write(*,*) 'DA order      : ', da_order
    end if
    write(*,*) 'Output prefix : ', trim(output_prefix)
    write(*,*) '----------------------------------------'

    ! Run propagation
    call run_uq_propagation( &
        nominal_state = state, &
        initial_cov   = cov, &
        epoch0        = epoch0, &
        t_start       = 0.0_DP, &
        t_end         = dt, &
        method_switch = method_switch, &
        n_particles   = n_particles, &
        save_results_to_file = .false., &
        da_order      = da_order, &
        initial_state_out = initial_state, &
        final_state_out   = final_state)

    write(*,*) '>>> Propagation complete.'

    ! Write outputs
    if (method_switch == METHOD_UT) then
        json_path = trim(output_prefix) // '_moments.json'
        write(*,*) '>>> Writing UT moments to: ', trim(json_path)
        call write_ut_json(json_path, final_state%mean, final_state%cov, method_str)
    else
        csv_path  = trim(output_prefix) // '_particles.csv'
        json_path = trim(output_prefix) // '_moments.json'

        write(*,*) '>>> Writing particles to: ', trim(csv_path)
        open(newunit=u_csv, file=trim(csv_path), status='replace', action='write')
        write(u_csv, '(A)') 'x,y,z,vx,vy,vz'
        do j = 1, size(final_state%samples, 2)
            write(u_csv, '(*(ES22.14, :, ","))') final_state%samples(:, j)
        end do
        close(u_csv)

        call final_state%compute_higher_moments(skewness, kurtosis)

        write(*,*) '>>> Writing moments to: ', trim(json_path)
        call write_mc_da_json(json_path, final_state%mean, final_state%cov, &
                               skewness, kurtosis, method_str)
    end if

    write(*,*) '=== UQ propagation finished successfully ==='

    ! Cleanup
    if (allocated(initial_state%samples)) deallocate(initial_state%samples)
    if (allocated(final_state%samples))   deallocate(final_state%samples)
    if (allocated(skewness))              deallocate(skewness)
    if (allocated(kurtosis))              deallocate(kurtosis)

contains

    subroutine write_ut_json(filename, mean_vec, cov_mat, method_label)
        character(len=*), intent(in) :: filename, method_label
        real(DP), intent(in) :: mean_vec(6), cov_mat(6,6)
        integer :: u, r

        open(newunit=u, file=trim(filename), status='replace', action='write')
        write(u, '(A)') '{'
        write(u, '(A,A,A)') '  "method": "', trim(method_label), '",'
        write(u, '(A,5(ES22.15,", "),ES22.15,A)') '  "mean": [', mean_vec, '],'
        write(u, '(A)') '  "covariance": ['
        do r = 1, 5
            write(u, '(A,5(ES22.15,", "),ES22.15,A)') '    [', cov_mat(r,:), '],'
        end do
        write(u, '(A,5(ES22.15,", "),ES22.15,A)') '    [', cov_mat(6,:), ']'
        write(u, '(A)') '  ]'
        write(u, '(A)') '}'
        close(u)
    end subroutine write_ut_json

    subroutine write_mc_da_json(filename, mean_vec, cov_mat, skew, kurt, method_label)
        character(len=*), intent(in) :: filename, method_label
        real(DP), intent(in) :: mean_vec(6), cov_mat(6,6), skew(6), kurt(6)
        integer :: u, r

        open(newunit=u, file=trim(filename), status='replace', action='write')
        write(u, '(A)') '{'
        write(u, '(A,A,A)') '  "method": "', trim(method_label), '",'
        write(u, '(A,5(ES22.15,", "),ES22.15,A)') '  "mean": [', mean_vec, '],'
        write(u, '(A)') '  "covariance": ['
        do r = 1, 5
            write(u, '(A,5(ES22.15,", "),ES22.15,A)') '    [', cov_mat(r,:), '],'
        end do
        write(u, '(A,5(ES22.15,", "),ES22.15,A)') '    [', cov_mat(6,:), ']'
        write(u, '(A)') '  ],'
        write(u, '(A,5(ES22.15,", "),ES22.15,A)') '  "marginal_skewness": [', skew, '],'
        write(u, '(A,5(ES22.15,", "),ES22.15,A)') '  "marginal_kurtosis": [', kurt, ']'
        write(u, '(A)') '}'
        close(u)
    end subroutine write_mc_da_json

end program run_HFEM_uprop
