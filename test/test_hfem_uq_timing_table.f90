!> @file test_hfem_uq_timing_table.f90
!> @brief Timing table for one 7-day HFEM MC estimate and DAMC runs.
program test_hfem_uq_timing_table
    use pod_global, only: DP, MAX_STRING_LEN
    use pod_engine_module, only: pod_engine_init
    use pod_data_format_module, only: load_initial_opm
    use pod_force_model_module, only: set_propagation_epoch
    use pod_config, only: config
    use pod_integrator_module, only: adaptive_step_integrate, METHOD_RKF78
    use pod_uq_propagation, only: run_uq_propagation, METHOD_DA
    use pod_uq_state_module, only: uq_state_type

    implicit none

    integer, parameter :: MAX_ORDERS = 16

    character(len=MAX_STRING_LEN) :: config_file, opm_file, arg
    real(DP) :: dt_seconds, epoch0, t_end_nondim
    real(DP) :: state(6), cov(6,6), current_state(6)
    real(DP) :: pr_km, vr_ms
    integer :: n_sample, n_orders, orders(MAX_ORDERS)
    integer :: i, n_steps, count_rate, count_start, count_end
    real(DP) :: cpu_start, cpu_end, wall_s, cpu_s, mc_case_wall_s, mc_case_cpu_s
    logical :: has_pr, has_vr
    real(DP), allocatable :: temp_times(:), temp_states(:,:)
    type(uq_state_type) :: initial_state, final_state

    config_file = 'config/config.txt'
    opm_file = 'OPM/L1Halo-1/L1Halo-1_init.opm.json'
    dt_seconds = 604800.0_DP
    n_sample = 100000
    orders = 0
    orders(1:3) = [3, 4, 6]
    n_orders = 3
    pr_km = 100.0_DP
    vr_ms = 0.3_DP
    has_pr = .true.
    has_vr = .true.

    call parse_args()

    call pod_engine_init(trim(config_file))
    call load_initial_opm(trim(opm_file), epoch0, state, cov)
    if (has_pr .and. has_vr) call build_isotropic_cov(pr_km, vr_ms, cov)
    call set_propagation_epoch(epoch0)

    write(*,'(A)') '========================================'
    write(*,'(A)') '  HFEM UQ 7-day timing table'
    write(*,'(A)') '========================================'
    write(*,'(A,A)') 'OPM file     : ', trim(opm_file)
    write(*,'(A,A)') 'Config file  : ', trim(config_file)
    write(*,'(A,F12.3)') 'dt seconds   : ', dt_seconds
    write(*,'(A,I0)') 'n_sample     : ', n_sample
    if (has_pr .and. has_vr) then
        write(*,'(A,F10.3,A,F10.6)') 'cov override : PR(km)=', pr_km, ' VR(m/s)=', vr_ms
    else
        write(*,'(A)') 'cov override : none, using OPM covariance'
    end if
    write(*,'(A)', advance='no') 'DA orders    : '
    do i = 1, n_orders
        if (i > 1) write(*,'(A)', advance='no') ','
        write(*,'(I0)', advance='no') orders(i)
    end do
    write(*,*)
    write(*,'(A)') '----------------------------------------'
    write(*,'(A)') 'method order wall_s cpu_s estimated_case_s estimated_case_h note'

    call time_mc_single_point(wall_s, cpu_s)
    mc_case_wall_s = wall_s * real(n_sample, DP)
    mc_case_cpu_s = cpu_s * real(n_sample, DP)
    write(*,'(A,1X,I0,1X,F12.6,1X,F12.6,1X,F14.3,1X,F12.6,1X,A)') &
        'MC', 0, wall_s, cpu_s, mc_case_wall_s, mc_case_wall_s / 3600.0_DP, &
        'single_point_times_n_sample'
    if (mc_case_cpu_s < 0.0_DP) error stop 'unreachable'

    do i = 1, n_orders
        call time_damc_order(orders(i), wall_s, cpu_s)
        write(*,'(A,1X,I0,1X,F12.6,1X,F12.6,1X,F14.3,1X,F12.6,1X,A)') &
            'DAMC', orders(i), wall_s, cpu_s, wall_s, wall_s / 3600.0_DP, 'measured_full_damc'
    end do

    write(*,'(A)') '========================================'
    write(*,*) 'test_hfem_uq_timing_table completed'

contains

    subroutine time_mc_single_point(wall_out, cpu_out)
        real(DP), intent(out) :: wall_out, cpu_out

        current_state = state
        current_state(1:3) = current_state(1:3) / config%LU
        current_state(4:6) = current_state(4:6) / config%VU
        t_end_nondim = dt_seconds / config%TU

        call system_clock(count_start, count_rate)
        call cpu_time(cpu_start)
        call adaptive_step_integrate(current_state, 0.0_DP, t_end_nondim, METHOD_RKF78, &
                                     temp_times, temp_states, n_steps)
        call cpu_time(cpu_end)
        call system_clock(count_end)

        wall_out = real(count_end - count_start, DP) / real(count_rate, DP)
        cpu_out = cpu_end - cpu_start

        if (allocated(temp_times)) deallocate(temp_times)
        if (allocated(temp_states)) deallocate(temp_states)
    end subroutine time_mc_single_point

    subroutine time_damc_order(order, wall_out, cpu_out)
        integer, intent(in) :: order
        real(DP), intent(out) :: wall_out, cpu_out

        call system_clock(count_start, count_rate)
        call cpu_time(cpu_start)
        call run_uq_propagation( &
            nominal_state = state, &
            initial_cov   = cov, &
            epoch0        = epoch0, &
            t_start       = 0.0_DP, &
            t_end         = dt_seconds, &
            method_switch = METHOD_DA, &
            n_particles   = n_sample, &
            save_results_to_file = .false., &
            da_order      = order, &
            initial_state_out = initial_state, &
            final_state_out   = final_state)
        call cpu_time(cpu_end)
        call system_clock(count_end)

        wall_out = real(count_end - count_start, DP) / real(count_rate, DP)
        cpu_out = cpu_end - cpu_start

        call initial_state%deallocate_memory()
        call final_state%deallocate_memory()
    end subroutine time_damc_order

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
                case ('-n', '--n-sample')
                    call get_command_argument(k + 1, arg)
                    read(arg, *) n_sample
                    k = k + 1
                case ('--orders')
                    call get_command_argument(k + 1, arg)
                    call parse_orders(arg)
                    k = k + 1
                case ('-pr')
                    call get_command_argument(k + 1, arg)
                    read(arg, *) pr_km
                    has_pr = .true.
                    k = k + 1
                case ('-vr')
                    call get_command_argument(k + 1, arg)
                    read(arg, *) vr_ms
                    has_vr = .true.
                    k = k + 1
                case ('--opm-cov')
                    has_pr = .false.
                    has_vr = .false.
                case default
                    write(*,'(A,A)') 'Warning: ignoring unknown argument: ', trim(arg)
            end select
            k = k + 1
        end do

        if (n_sample < 1) n_sample = 1
        if (n_orders < 1) then
            n_orders = 3
            orders(1:3) = [3, 4, 6]
        end if
        if (has_pr .neqv. has_vr) error stop '-pr and -vr must be used together'
    end subroutine parse_args

    subroutine parse_orders(text)
        character(len=*), intent(in) :: text
        character(len=MAX_STRING_LEN) :: token
        integer :: pos, start, len_text, value

        orders = 0
        n_orders = 0
        len_text = len_trim(text)
        start = 1
        do
            pos = index(text(start:len_text), ',')
            if (pos == 0) then
                token = adjustl(text(start:len_text))
                if (len_trim(token) > 0) then
                    read(token, *) value
                    call add_order(value)
                end if
                exit
            else
                token = adjustl(text(start:start + pos - 2))
                if (len_trim(token) > 0) then
                    read(token, *) value
                    call add_order(value)
                end if
                start = start + pos
            end if
            if (start > len_text) exit
        end do
    end subroutine parse_orders

    subroutine add_order(value)
        integer, intent(in) :: value

        if (value <= 0) error stop 'DA order must be positive'
        if (n_orders >= MAX_ORDERS) error stop 'too many DA orders'
        n_orders = n_orders + 1
        orders(n_orders) = value
    end subroutine add_order

    subroutine print_usage()
        write(*,'(A)') 'Usage: fpm test test_hfem_uq_timing_table -- [options]'
        write(*,'(A)') 'Options:'
        write(*,'(A)') '  -cfg <file>        Config file, default config/config.txt'
        write(*,'(A)') '  -opm <file>        Initial OPM, default OPM/L1Halo-1/L1Halo-1_init.opm.json'
        write(*,'(A)') '  -dt <seconds>      Propagation duration, default 604800'
        write(*,'(A)') '  -n <count>         n_sample for MC estimate and DAMC samples, default 100000'
        write(*,'(A)') '  --orders <csv>     DA orders, default 3,4,6; use 2,4,6 if desired'
        write(*,'(A)') '  -pr <km> -vr <m/s> Covariance override, default PR=100 km VR=0.3 m/s'
        write(*,'(A)') '  --opm-cov          Use covariance from OPM instead of -pr/-vr override'
    end subroutine print_usage

    subroutine build_isotropic_cov(pr_km_in, vr_ms_in, cov_out)
        real(DP), intent(in)  :: pr_km_in, vr_ms_in
        real(DP), intent(out) :: cov_out(6,6)
        real(DP) :: vr_kms, pos_var, vel_var
        integer :: j

        vr_kms = vr_ms_in / 1000.0_DP
        pos_var = (pr_km_in ** 2) / 3.0_DP
        vel_var = (vr_kms ** 2) / 3.0_DP

        cov_out = 0.0_DP
        do j = 1, 3
            cov_out(j, j) = pos_var
            cov_out(j + 3, j + 3) = vel_var
        end do
    end subroutine build_isotropic_cov

end program test_hfem_uq_timing_table