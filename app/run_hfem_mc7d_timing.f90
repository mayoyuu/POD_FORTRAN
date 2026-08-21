!> @file run_hfem_mc7d_timing.f90
!> @brief Time one HFEM Monte-Carlo particle propagation and extrapolate batch cost.
program run_hfem_mc7d_timing
    use pod_global, only: DP, MAX_STRING_LEN
    use pod_engine_module, only: pod_engine_init
    use pod_data_format_module, only: load_initial_opm
    use pod_force_model_module, only: set_propagation_epoch
    use pod_config, only: config
    use pod_integrator_module, only: adaptive_step_integrate, METHOD_RKF45, METHOD_RKF78

    implicit none

    character(len=MAX_STRING_LEN) :: config_file, opm_file, arg
    real(DP) :: dt_seconds, epoch0, t_end_nondim
    real(DP) :: state(6), cov(6,6), current_state(6)
    integer :: repeats, particles_per_case, mc_cases, active_parallel
    integer :: integrator_type, i, n_steps, count_rate, count_start, count_end
    real(DP) :: cpu_start, cpu_end, wall_elapsed, cpu_elapsed
    real(DP) :: total_wall, total_cpu, per_wall, per_cpu
    real(DP) :: case_wall, case_cpu, waves, batch_wall
    real(DP), allocatable :: temp_times(:), temp_states(:,:)

    config_file = 'config/config.txt'
    opm_file = 'OPM/L1Halo-1/L1Halo-1_init.opm.json'
    dt_seconds = 604800.0_DP
    repeats = 1
    particles_per_case = 100000
    mc_cases = 24
    active_parallel = 6
    integrator_type = METHOD_RKF78

    call parse_args()

    call pod_engine_init(trim(config_file))
    call load_initial_opm(trim(opm_file), epoch0, state, cov)
    call set_propagation_epoch(epoch0)

    t_end_nondim = dt_seconds / config%TU
    total_wall = 0.0_DP
    total_cpu = 0.0_DP

    write(*,'(A)') '========================================'
    write(*,'(A)') '  HFEM MC 7-day single-particle timing'
    write(*,'(A)') '========================================'
    write(*,'(A,A)') 'OPM file          : ', trim(opm_file)
    write(*,'(A,A)') 'Config file       : ', trim(config_file)
    write(*,'(A,F12.3)') 'dt seconds        : ', dt_seconds
    write(*,'(A,I0)') 'repeats           : ', repeats
    write(*,'(A,I0)') 'particles / case  : ', particles_per_case
    write(*,'(A,I0)') 'MC cases          : ', mc_cases
    write(*,'(A,I0)') 'active parallel   : ', active_parallel
    write(*,'(A,A)') 'integrator        : ', integrator_name(integrator_type)
    write(*,'(A)') '----------------------------------------'

    do i = 1, repeats
        current_state = state
        current_state(1:3) = current_state(1:3) / config%LU
        current_state(4:6) = current_state(4:6) / config%VU

        call system_clock(count_start, count_rate)
        call cpu_time(cpu_start)
        call adaptive_step_integrate(current_state, 0.0_DP, t_end_nondim, integrator_type, &
                                     temp_times, temp_states, n_steps)
        call cpu_time(cpu_end)
        call system_clock(count_end)

        wall_elapsed = real(count_end - count_start, DP) / real(count_rate, DP)
        cpu_elapsed = cpu_end - cpu_start
        total_wall = total_wall + wall_elapsed
        total_cpu = total_cpu + cpu_elapsed

        current_state = temp_states(n_steps, :)
        current_state(1:3) = current_state(1:3) * config%LU
        current_state(4:6) = current_state(4:6) * config%VU

        write(*,'(A,I0,A,F12.6,A,F12.6,A,I0)') 'repeat ', i, ': wall=', wall_elapsed, &
            ' s, cpu=', cpu_elapsed, ' s, steps=', n_steps
        write(*,'(A,6(ES14.5,1X))') '  final state: ', current_state

        if (allocated(temp_times)) deallocate(temp_times)
        if (allocated(temp_states)) deallocate(temp_states)
    end do

    per_wall = total_wall / real(repeats, DP)
    per_cpu = total_cpu / real(repeats, DP)
    case_wall = per_wall * real(particles_per_case, DP)
    case_cpu = per_cpu * real(particles_per_case, DP)
    waves = ceiling(real(mc_cases, DP) / real(max(active_parallel, 1), DP))
    batch_wall = case_wall * waves

    write(*,'(A)') '----------------------------------------'
    write(*,'(A,F12.6,A)') 'per particle wall : ', per_wall, ' s'
    write(*,'(A,F12.6,A)') 'per particle cpu  : ', per_cpu, ' s'
    write(*,'(A,F12.3,A,F10.3,A)') 'one MC case wall  : ', case_wall, ' s = ', case_wall/3600.0_DP, ' h'
    write(*,'(A,F12.3,A,F10.3,A)') 'one MC case cpu   : ', case_cpu, ' s = ', case_cpu/3600.0_DP, ' h'
    write(*,'(A,F12.3,A,F10.3,A)') 'batch wall est.   : ', batch_wall, ' s = ', batch_wall/3600.0_DP, ' h'
    write(*,'(A)') '========================================'

contains

    subroutine parse_args()
        integer :: n, k

        n = command_argument_count()
        k = 1
        do while (k <= n)
            call get_command_argument(k, arg)
            select case (trim(arg))
                case ('-h', '--help')
                    call print_usage()
                    stop 0
                case ('-cfg', '--config')
                    call get_command_argument(k + 1, config_file)
                    k = k + 1
                case ('-opm')
                    call get_command_argument(k + 1, opm_file)
                    k = k + 1
                case ('-dt')
                    call get_command_argument(k + 1, arg)
                    read(arg, *) dt_seconds
                    k = k + 1
                case ('-r', '--repeats')
                    call get_command_argument(k + 1, arg)
                    read(arg, *) repeats
                    k = k + 1
                case ('-n', '--particles')
                    call get_command_argument(k + 1, arg)
                    read(arg, *) particles_per_case
                    k = k + 1
                case ('--mc-cases')
                    call get_command_argument(k + 1, arg)
                    read(arg, *) mc_cases
                    k = k + 1
                case ('--parallel')
                    call get_command_argument(k + 1, arg)
                    read(arg, *) active_parallel
                    k = k + 1
                case ('--rkf45')
                    integrator_type = METHOD_RKF45
                case ('--rkf78')
                    integrator_type = METHOD_RKF78
                case default
                    write(*,'(A,A)') 'Warning: ignoring unknown argument: ', trim(arg)
            end select
            k = k + 1
        end do

        if (repeats < 1) repeats = 1
        if (particles_per_case < 1) particles_per_case = 1
        if (mc_cases < 1) mc_cases = 1
        if (active_parallel < 1) active_parallel = 1
    end subroutine parse_args

    subroutine print_usage()
        write(*,'(A)') 'Usage: fpm run run_hfem_mc7d_timing -- [options]'
        write(*,'(A)') 'Options:'
        write(*,'(A)') '  -cfg <file>        Config file, default config/config.txt'
        write(*,'(A)') '  -opm <file>        Initial OPM, default OPM/L1Halo-1/L1Halo-1_init.opm.json'
        write(*,'(A)') '  -dt <seconds>      Propagation duration, default 604800'
        write(*,'(A)') '  -r <count>         Repeated single-particle timings, default 1'
        write(*,'(A)') '  -n <particles>     Particles per MC case, default 100000'
        write(*,'(A)') '  --mc-cases <count> MC cases in batch, default 24'
        write(*,'(A)') '  --parallel <count> Concurrent MC cases, default 6'
        write(*,'(A)') '  --rkf45|--rkf78    Integrator, default RKF78'
    end subroutine print_usage

    function integrator_name(integ) result(name)
        integer, intent(in) :: integ
        character(len=8) :: name

        select case (integ)
            case (METHOD_RKF45)
                name = 'RKF45'
            case default
                name = 'RKF78'
        end select
    end function integrator_name

end program run_hfem_mc7d_timing
