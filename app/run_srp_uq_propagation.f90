program run_srp_uq_propagation
    use pod_global, only: DP, MAX_STRING_LEN
    use pod_engine_module, only: pod_engine_init
    use pod_config, only: config
    use pod_spice, only: str2et
    use pod_data_format_module, only: load_initial_opm
    use pod_random_module, only: init_random_seed, generate_multivariate_normal, randn
    use pod_uq_propagation, only: run_particle_propagation, METHOD_DA
    use pod_uq_state_module, only: uq_state_type
    use pod_da_force_model_module, only: set_srp_ballistic_parameters
    implicit none

    character(len=MAX_STRING_LEN) :: opm_file, output_prefix, config_file
    character(len=MAX_STRING_LEN) :: epoch_str, arg_str, csv_path, json_path
    real(DP) :: epoch0, target_et, dt_seconds, dt
    real(DP) :: state6(6), cov6(6,6), ref_out(6)
    real(DP) :: srp_cr, srp_smr, srp_rp, srp_eta_mean, srp_eta_sigma
    integer :: n_particles, da_order, num_args, i, j, ext_pos, unit_csv
    logical :: has_opm, has_dt, has_et, has_output
    type(uq_state_type) :: initial_state, final_state
    real(DP), allocatable :: orbit_samples(:,:)

    config_file = 'config/config.txt'
    output_prefix = ''
    opm_file = ''
    epoch_str = ''
    n_particles = 100000
    da_order = 4
    dt_seconds = 0.0_DP
    target_et = 0.0_DP
    srp_cr = 1.25_DP
    srp_smr = 7.5e-3_DP
    srp_rp = 1367.0_DP / 299792458.0_DP
    srp_eta_mean = 0.0_DP
    srp_eta_sigma = 0.0_DP
    has_opm = .false.
    has_dt = .false.
    has_et = .false.
    has_output = .false.

    num_args = command_argument_count()
    i = 1
    do while (i <= num_args)
        call get_command_argument(i, arg_str)
        select case (trim(arg_str))
        case ('-opm', '--opm')
            call get_command_argument(i+1, opm_file)
            has_opm = .true.
            i = i + 1
        case ('-dt', '--dt')
            call get_command_argument(i+1, arg_str)
            read(arg_str, *) dt_seconds
            has_dt = .true.
            i = i + 1
        case ('-et', '--epoch')
            call get_command_argument(i+1, epoch_str)
            has_et = .true.
            i = i + 1
        case ('-o', '--output')
            call get_command_argument(i+1, output_prefix)
            has_output = .true.
            i = i + 1
        case ('-n', '--n-particles')
            call get_command_argument(i+1, arg_str)
            read(arg_str, *) n_particles
            i = i + 1
        case ('-da', '--da-order')
            call get_command_argument(i+1, arg_str)
            read(arg_str, *) da_order
            i = i + 1
        case ('-cr', '--srp-cr')
            call get_command_argument(i+1, arg_str)
            read(arg_str, *) srp_cr
            i = i + 1
        case ('-smr', '--srp-smr')
            call get_command_argument(i+1, arg_str)
            read(arg_str, *) srp_smr
            i = i + 1
        case ('-rp', '--srp-rp')
            call get_command_argument(i+1, arg_str)
            read(arg_str, *) srp_rp
            i = i + 1
        case ('-srp-mean', '--eta-mean')
            call get_command_argument(i+1, arg_str)
            read(arg_str, *) srp_eta_mean
            i = i + 1
        case ('-srp-sigma', '--eta-sigma')
            call get_command_argument(i+1, arg_str)
            read(arg_str, *) srp_eta_sigma
            i = i + 1
        case ('-cfg', '--config')
            call get_command_argument(i+1, config_file)
            i = i + 1
        case ('-h', '--help')
            call print_usage()
            stop
        case default
            write(*,*) 'Warning: ignoring unknown argument: ', trim(arg_str)
        end select
        i = i + 1
    end do

    if (.not. has_opm) then
        write(*,*) 'Error: -opm <file.opm.json> is required.'
        call print_usage()
        stop 1
    end if
    if (.not. has_dt .and. .not. has_et) then
        write(*,*) 'Error: provide either -dt <seconds> or -et <UTC epoch>.'
        call print_usage()
        stop 1
    end if
    if (has_dt .and. has_et) then
        write(*,*) 'Error: -dt and -et are mutually exclusive.'
        stop 1
    end if
    if (n_particles < 2) then
        write(*,*) 'Error: -n must be at least 2.'
        stop 1
    end if
    if (srp_eta_sigma < 0.0_DP) then
        write(*,*) 'Error: -srp-sigma must be non-negative.'
        stop 1
    end if

    if (.not. has_output) then
        ext_pos = index(opm_file, '.', back=.true.)
        if (ext_pos > 1) then
            output_prefix = opm_file(1:ext_pos-1) // '_srp_uq'
        else
            output_prefix = trim(opm_file) // '_srp_uq'
        end if
    end if

    write(*,*) '>>> Initializing POD engine...'
    call pod_engine_init(trim(config_file))
    config%use_srp = .true.
    call set_srp_ballistic_parameters(srp_cr, srp_smr, srp_rp)

    write(*,*) '>>> Loading OPM: ', trim(opm_file)
    call load_initial_opm(trim(opm_file), epoch0, state6, cov6)

    if (has_dt) then
        dt = dt_seconds
    else
        call str2et(trim(epoch_str), target_et)
        dt = target_et - epoch0
    end if
    if (dt <= 0.0_DP) then
        write(*,*) 'Error: propagation duration must be positive.'
        stop 1
    end if

    write(*,*) '========================================'
    write(*,*) '  SRP-Parameter DA-MC Propagation'
    write(*,*) '========================================'
    write(*,*) 'OPM file      : ', trim(opm_file)
    write(*,*) 'Config file   : ', trim(config_file)
    write(*,*) 'Epoch0 (TDB)  : ', epoch0
    write(*,*) 'dt (seconds)  : ', dt
    write(*,*) 'n_particles   : ', n_particles
    write(*,*) 'DA order      : ', da_order
    write(*,*) 'SRP Cr        : ', srp_cr
    write(*,*) 'SRP SMR       : ', srp_smr
    write(*,*) 'SRP RP        : ', srp_rp
    write(*,*) 'eta_srp mean  : ', srp_eta_mean
    write(*,*) 'eta_srp sigma : ', srp_eta_sigma
    write(*,*) 'Output prefix : ', trim(output_prefix)
    write(*,*) '----------------------------------------'

    allocate(orbit_samples(6, n_particles))
    call init_random_seed(.true.)
    call generate_multivariate_normal(state6, cov6, orbit_samples)

    call initial_state%allocate_memory(7, n_particles)
    initial_state%samples(1:6, :) = orbit_samples
    do j = 1, n_particles
        if (srp_eta_sigma > 0.0_DP) then
            initial_state%samples(7, j) = srp_eta_mean + srp_eta_sigma * randn()
        else
            initial_state%samples(7, j) = srp_eta_mean
        end if
    end do

    call run_particle_propagation(initial_state, state6, epoch0, 0.0_DP, dt, &
                                  METHOD_DA, final_state, da_order=da_order, &
                                  reference_orbit_out=ref_out)
    call final_state%compute_moments()

    csv_path = trim(output_prefix) // '_particles.csv'
    json_path = trim(output_prefix) // '_moments.json'

    open(newunit=unit_csv, file=trim(csv_path), status='replace', action='write')
    write(unit_csv, '(A)') 'x,y,z,vx,vy,vz,eta_srp'
    do j = 1, n_particles
        write(unit_csv, '(*(ES22.14, :, ","))') final_state%samples(:, j)
    end do
    close(unit_csv)

    call write_moments_json(json_path, final_state%mean, final_state%cov, &
                            srp_cr, srp_smr, srp_rp, srp_eta_mean, srp_eta_sigma, &
                            epoch0, epoch0 + dt, n_particles, da_order)

    write(*,*) '>>> Results saved:'
    write(*,*) '  ', trim(csv_path)
    write(*,*) '  ', trim(json_path)

    call initial_state%deallocate_memory()
    call final_state%deallocate_memory()
    if (allocated(orbit_samples)) deallocate(orbit_samples)

contains

    subroutine print_usage()
        write(*,*) 'Usage:'
        write(*,*) '  fpm run run_srp_uq_propagation -- -opm <file.opm.json> -dt <seconds> -o <prefix> [options]'
        write(*,*) '  fpm run run_srp_uq_propagation -- -opm <file.opm.json> -et <UTC epoch> -o <prefix> [options]'
        write(*,*) 'Options:'
        write(*,*) '  -cfg <file>          Config file, default config/config.txt'
        write(*,*) '  -n <N>               Particle count, default 100000'
        write(*,*) '  -da <order>          DA order, default 4'
        write(*,*) '  -cr <Cr>             SRP reflectivity coefficient, default 1.25'
        write(*,*) '  -smr <A_over_m>      Area-to-mass ratio in m^2/kg, default 7.5e-3'
        write(*,*) '  -rp <pressure>       Solar radiation pressure at 1 AU in N/m^2'
        write(*,*) '  -srp-mean <eta>      Mean of eta_srp, default 0'
        write(*,*) '  -srp-sigma <sigma>   1-sigma of eta_srp, default 0'
    end subroutine print_usage

    subroutine write_moments_json(filename, mean_vec, cov_mat, cr, smr, rp, eta_mean, eta_sigma, &
                                  et0, etf, npart, order)
        character(len=*), intent(in) :: filename
        real(DP), intent(in) :: mean_vec(:), cov_mat(:,:)
        real(DP), intent(in) :: cr, smr, rp, eta_mean, eta_sigma, et0, etf
        integer, intent(in) :: npart, order
        integer :: u

        open(newunit=u, file=trim(filename), status='replace', action='write')
        write(u, '(A)') '{'
        write(u, '(A)') '  "method": "DA-MC-SRP",'
        write(u, '(A,I0,A)') '  "n_particles": ', npart, ','
        write(u, '(A,I0,A)') '  "da_order": ', order, ','
        write(u, '(A,ES24.15E3,A)') '  "epoch0_et": ', et0, ','
        write(u, '(A,ES24.15E3,A)') '  "epoch_final_et": ', etf, ','
        write(u, '(A)') '  "state_labels": ["x", "y", "z", "vx", "vy", "vz", "eta_srp"],'
        write(u, '(A)') '  "srp": {'
        write(u, '(A,ES24.15E3,A)') '    "Cr": ', cr, ','
        write(u, '(A,ES24.15E3,A)') '    "SMR": ', smr, ','
        write(u, '(A,ES24.15E3,A)') '    "RP": ', rp, ','
        write(u, '(A,ES24.15E3,A)') '    "eta_mean_input": ', eta_mean, ','
        write(u, '(A,ES24.15E3)')  '    "eta_sigma_input": ', eta_sigma
        write(u, '(A)') '  },'
        call write_json_vector(u, 'mean', mean_vec, .true.)
        call write_json_matrix(u, 'covariance', cov_mat)
        write(u, '(A)') '}'
        close(u)
    end subroutine write_moments_json

    subroutine write_json_vector(u, name, vec, trailing_comma)
        integer, intent(in) :: u
        character(len=*), intent(in) :: name
        real(DP), intent(in) :: vec(:)
        logical, intent(in) :: trailing_comma
        integer :: k

        write(u, '(A,A,A)', advance='no') '  "', trim(name), '": ['
        do k = 1, size(vec)
            if (k < size(vec)) then
                write(u, '(ES24.15E3,A)', advance='no') vec(k), ', '
            else
                write(u, '(ES24.15E3)', advance='no') vec(k)
            end if
        end do
        if (trailing_comma) then
            write(u, '(A)') '],'
        else
            write(u, '(A)') ']'
        end if
    end subroutine write_json_vector

    subroutine write_json_matrix(u, name, mat)
        integer, intent(in) :: u
        character(len=*), intent(in) :: name
        real(DP), intent(in) :: mat(:,:)
        integer :: r, c

        write(u, '(A,A,A)') '  "', trim(name), '": ['
        do r = 1, size(mat, 1)
            write(u, '(A)', advance='no') '    ['
            do c = 1, size(mat, 2)
                if (c < size(mat, 2)) then
                    write(u, '(ES24.15E3,A)', advance='no') mat(r, c), ', '
                else
                    write(u, '(ES24.15E3)', advance='no') mat(r, c)
                end if
            end do
            if (r < size(mat, 1)) then
                write(u, '(A)') '],'
            else
                write(u, '(A)') ']'
            end if
        end do
        write(u, '(A)') '  ]'
    end subroutine write_json_matrix

end program run_srp_uq_propagation
