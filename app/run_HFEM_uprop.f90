!> @file run_uq_propagation.f90
!> @brief 统一的不确定性传播 CLI 入口，支持 MC/DA/UT/ADS
!> @author Song Yu
!> @date 2026-05-26
!>
!> Usage:
!>   fpm run run_HFEM_uprop -- -opm <file> -m <MC|DA|UT> -dt <seconds> -o <prefix>
!>   fpm run run_HFEM_uprop -- -opm <file> -m <MC|DA|UT> -et <epoch> -n 5000 -o <prefix>
!>   Covariance override (isotropic diagonal):
!>     -pr <km>   3D position RMS in km  -> cov_ii = PR^2/3 for i=1..3
!>     -vr <m/s>  3D velocity RMS in m/s -> cov_ii = (VR/1000)^2/3 for i=4..6
!>     Must be used together (-pr and -vr).
!>   --zero-init-cov sets the 6D initial orbit covariance to exactly zero.
!>   ADS: --ads-coordinates component|whitened --ads-domain-sigma 3
!>        --srp-sigma <1-sigma> --perturb-file <dx,dy,dz,dvx,dvy,dvz[,eta] CSV>
!> Example: fpm run run_HFEM_uprop -- -opm input/init.opm -m ADS &
!>   -dt 86400 -o output/ads_run --srp-sigma 0.02
!> Output:
!>   MC/DA: <prefix>_particles.csv + <prefix>_moments.json (mean/cov/skewness/kurtosis)
!>   UT:    <prefix>_moments.json (mean/cov)
!>   ADS:   the particle/moments files plus <prefix>_ads_stats.json
program run_HFEM_uprop
    use pod_global, only: DP, MAX_STRING_LEN
    use pod_engine_module, only: pod_engine_init
    use pod_dace_classes, only: dace_initialize
    use pod_spice, only: str2et
    use pod_data_format_module, only: load_initial_opm
    use pod_uq_propagation, only: run_uq_propagation, METHOD_MC, METHOD_DA, &
        METHOD_UT, METHOD_ADS
    use pod_uq_state_module, only: uq_state_type
    use pod_uq_ads_coordinates_module, only: ADS_COORD_COMPONENT, ADS_COORD_WHITENED
    use pod_uq_hfem_ads_module, only: hfem_ads_options_type, hfem_ads_stats_type
    use pod_uq_sample_io_module, only: read_uq_perturbations_csv, &
        write_uq_particles_csv, write_uq_moments_json, write_ads_stats_json

    implicit none

    character(len=MAX_STRING_LEN) :: opm_file, method_str, output_prefix, config_file
    character(len=MAX_STRING_LEN) :: epoch_str, arg_str
    character(len=MAX_STRING_LEN) :: json_path, csv_path, ads_stats_path
    character(len=MAX_STRING_LEN) :: ads_coordinate_str, io_message
    character(len=MAX_STRING_LEN) :: perturb_file
    real(DP) :: pr_km, vr_ms, dt_seconds, t_end_et, epoch0, dt
    real(DP) :: ads_pos_tol, ads_vel_tol
    integer  :: method_switch, n_particles, da_order, i, num_args, ext_pos
    integer  :: io_status
    logical  :: has_dt, has_et, has_opm, has_method, has_output, has_pr, has_vr
    logical  :: zero_init_cov, has_perturb_file

    type(uq_state_type) :: initial_state, final_state
    type(hfem_ads_options_type) :: ads_options
    type(hfem_ads_stats_type) :: ads_stats
    real(DP), allocatable :: skewness(:), kurtosis(:)
    real(DP), allocatable :: initial_perturbations(:,:)
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
    has_pr       = .false.
    has_vr       = .false.
    zero_init_cov = .false.
    has_perturb_file = .false.
    dt_seconds   = 0.0_DP
    t_end_et     = 0.0_DP
    pr_km        = 0.0_DP
    vr_ms        = 0.0_DP
    epoch_str    = ''
    output_prefix = ''
    perturb_file = ''
    ads_coordinate_str = 'component'
    ads_pos_tol = ads_options%error_tolerance(1)
    ads_vel_tol = ads_options%error_tolerance(4)

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

            case ('-pr')
                call get_command_argument(i+1, arg_str)
                read(arg_str, *) pr_km
                has_pr = .true.
                i = i + 1

            case ('-vr')
                call get_command_argument(i+1, arg_str)
                read(arg_str, *) vr_ms
                has_vr = .true.
                i = i + 1

            case ('--zero-init-cov')
                zero_init_cov = .true.

            case ('-cfg', '--config')
                call get_command_argument(i+1, arg_str)
                config_file = trim(arg_str)
                i = i + 1

            case ('--ads-coordinates')
                call get_command_argument(i+1, arg_str)
                ads_coordinate_str = trim(arg_str)
                i = i + 1

            case ('--ads-domain-sigma')
                call get_command_argument(i+1, arg_str)
                read(arg_str, *) ads_options%domain_sigma
                i = i + 1

            case ('--ads-max-depth')
                call get_command_argument(i+1, arg_str)
                read(arg_str, *) ads_options%max_split_depth
                i = i + 1

            case ('--ads-pos-tol')
                call get_command_argument(i+1, arg_str)
                read(arg_str, *) ads_pos_tol
                i = i + 1

            case ('--ads-vel-tol')
                call get_command_argument(i+1, arg_str)
                read(arg_str, *) ads_vel_tol
                i = i + 1

            case ('--srp-sigma')
                call get_command_argument(i+1, arg_str)
                read(arg_str, *) ads_options%srp_sigma
                i = i + 1

            case ('--perturb-file')
                call get_command_argument(i+1, arg_str)
                perturb_file = trim(arg_str)
                has_perturb_file = .true.
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
        case ('ADS', 'ads')
            method_switch = METHOD_ADS
        case default
            write(*,*) 'Error: unknown method "', trim(method_str), &
                '". Use MC, DA, UT, or ADS.'
            stop 1
    end select

    select case (trim(ads_coordinate_str))
        case ('component', 'COMPONENT')
            ads_options%coordinate_mode = ADS_COORD_COMPONENT
        case ('whitened', 'WHITENED')
            ads_options%coordinate_mode = ADS_COORD_WHITENED
        case default
            write(*,*) 'Error: --ads-coordinates must be component or whitened.'
            stop 1
    end select
    ads_options%da_order = da_order
    ads_options%error_tolerance(1:3) = ads_pos_tol
    ads_options%error_tolerance(4:6) = ads_vel_tol

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
    if (method_switch == METHOD_DA .or. method_switch == METHOD_ADS) then
        write(*,*) '>>> Initializing DACE for DA propagation...'
        call dace_initialize(da_order, merge(7,6, &
            method_switch == METHOD_ADS .and. ads_options%srp_sigma > 0.0_DP))
        write(*,*) '>>> DACE initialized for DA/ADS propagation.'
    end if

    ! Load OPM
    write(*,*) '>>> Loading OPM file: ', trim(opm_file)
    call load_initial_opm(trim(opm_file), epoch0, state, cov)

    ! Override covariance from -pr / -vr if provided
    if (has_pr .and. has_vr) then
        call build_isotropic_cov(pr_km, vr_ms, cov)
        write(*,*) '>>> Covariance overridden by -pr / -vr'
    else if (has_pr .neqv. has_vr) then
        write(*,*) 'Error: -pr and -vr must be used together.'
        stop 1
    end if

    if (zero_init_cov) then
        cov = 0.0_DP
        write(*,*) 'Initial orbit covariance: ZERO'
    end if

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
    if (method_switch == METHOD_DA .or. method_switch == METHOD_ADS) then
        write(*,*) 'DA order      : ', da_order
    end if
    if (has_perturb_file) then
        if (method_switch /= METHOD_ADS) then
            write(*,*) 'Error: --perturb-file is only supported by ADS.'
            stop 1
        end if
        call read_uq_perturbations_csv(perturb_file, initial_perturbations, &
            io_status, io_message)
        if (io_status /= 0) error stop trim(io_message)
        n_particles = size(initial_perturbations,2)
    end if
    if (method_switch == METHOD_ADS) then
        write(*,*) 'ADS coordinates: ', trim(ads_coordinate_str)
        write(*,*) 'ADS domain sigma: ', ads_options%domain_sigma
        write(*,*) 'SRP sigma       : ', ads_options%srp_sigma
    end if
    if (has_pr .and. has_vr) then
        write(*,*) 'PR (3D RMS km): ', pr_km
        write(*,*) 'VR (3D RMS m/s): ', vr_ms
    end if
    write(*,*) 'Output prefix : ', trim(output_prefix)
    write(*,*) '----------------------------------------'

    ! Run propagation
    if (has_perturb_file) then
        call run_uq_propagation(state, cov, epoch0, 0.0_DP, dt, method_switch, &
            n_particles, .false., initial_state, final_state, da_order=da_order, &
            ads_options=ads_options, ads_stats=ads_stats, &
            initial_perturbations=initial_perturbations)
    else
        call run_uq_propagation(state, cov, epoch0, 0.0_DP, dt, method_switch, &
            n_particles, .false., initial_state, final_state, da_order=da_order, &
            ads_options=ads_options, ads_stats=ads_stats)
    end if

    write(*,*) '>>> Propagation complete.'

    ! Write outputs
    if (method_switch == METHOD_UT) then
        json_path = trim(output_prefix) // '_moments.json'
        write(*,*) '>>> Writing UT moments to: ', trim(json_path)
        call write_ut_json(json_path, final_state%mean, final_state%cov, method_str)
    else if (method_switch == METHOD_ADS) then
        csv_path = trim(output_prefix)//'_particles.csv'
        json_path = trim(output_prefix)//'_moments.json'
        ads_stats_path = trim(output_prefix)//'_ads_stats.json'
        call write_uq_particles_csv(csv_path, final_state%samples, io_status, io_message)
        if (io_status /= 0) error stop trim(io_message)
        call write_uq_moments_json(json_path, final_state, method_str, &
            io_status, io_message)
        if (io_status /= 0) error stop trim(io_message)
        ads_stats%written_count = size(final_state%samples,2)
        call write_ads_stats_json(ads_stats_path, ads_stats, io_status, io_message)
        if (io_status /= 0) error stop trim(io_message)
        write(*,*) '>>> ADS particles: ', size(final_state%samples,2)
        write(*,*) '>>> ADS patch count: ', ads_stats%n_patches
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

    subroutine build_isotropic_cov(pr_km, vr_ms, cov)
        real(DP), intent(in)  :: pr_km, vr_ms
        real(DP), intent(out) :: cov(6,6)
        real(DP) :: vr_kms, pos_var, vel_var
        integer  :: k

        vr_kms = vr_ms / 1000.0_DP
        pos_var = (pr_km ** 2) / 3.0_DP
        vel_var = (vr_kms ** 2) / 3.0_DP

        cov = 0.0_DP
        do k = 1, 3
            cov(k, k) = pos_var
            cov(k+3, k+3) = vel_var
        end do
    end subroutine build_isotropic_cov

end program run_HFEM_uprop
