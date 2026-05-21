program test_da_convergence
    use pod_global, only: DP
    use pod_engine_module, only: pod_engine_init
    use pod_config, only: config
    use pod_spice, only: str2et
    use pod_da_force_model_module, only: init_gravity_network
    use pod_dace_classes, only: dace_initialize
    use pod_data_format_module, only: load_initial_opm
    use pod_uq_state_module, only: uq_state_type
    use pod_random_module, only: init_random_seed, generate_multivariate_normal
    use pod_uq_propagation, only: run_particle_propagation, METHOD_MC, METHOD_DA

    implicit none

    character(len=*), parameter :: CONFIG_FILE = 'config/dummy_test_config.txt'
    character(len=*), parameter :: OPM_FILE    = 'OPM/L1Halo/L1Halo_init.opm.json'
    integer, parameter :: DA_MAX_ORDER = 5
    integer, parameter :: DA_NVARS     = 6
    integer, parameter :: N_PARTICLES  = 1000
    real(DP), parameter :: DURATION    = 14.0_DP * 86400.0_DP  ! 14 days in seconds

    ! Local variables
    real(DP) :: nominal_state(6), initial_cov(6,6)
    real(DP) :: epoch0
    type(uq_state_type) :: init_state
    real(DP), allocatable :: mc_golden(:,:)
    real(DP), allocatable :: da_best(:,:)

    write(*,*) '============================================================'
    write(*,*) '  DA Order Convergence Test - L1Halo / 14-day / RKF45'
    write(*,*) '============================================================'
    write(*,*) ''

    ! Phase 0: System init
    write(*,*) '>>> Phase 0: System initialization...'
    call pod_engine_init(CONFIG_FILE)
    call dace_initialize(DA_MAX_ORDER, DA_NVARS)
    write(*,*) '  DACE initialized with max_order=', DA_MAX_ORDER, ', nvars=', DA_NVARS
    write(*,*) ''

    ! Phase 1: Read OPM
    write(*,*) '>>> Phase 1: Reading OPM from ', OPM_FILE
    call load_initial_opm(OPM_FILE, epoch0, nominal_state, initial_cov)
    write(*,*) '  Epoch (TDB sec):  ', epoch0
    write(*,*) '  Nominal state:    ', nominal_state
    write(*,'(A,6ES12.4)') '  Position sigma (km):  ', &
        sqrt(initial_cov(1,1)), sqrt(initial_cov(2,2)), sqrt(initial_cov(3,3))
    write(*,'(A,6ES12.4)') '  Velocity sigma (km/s):', &
        sqrt(initial_cov(4,4)), sqrt(initial_cov(5,5)), sqrt(initial_cov(6,6))
    write(*,*) ''

    ! Phase 2: Generate 1000 particles
    write(*,*) '>>> Phase 2: Generating ', N_PARTICLES, ' initial particles...'
    call init_state%allocate_memory(6, N_PARTICLES)
    call init_random_seed(.true.)
    call generate_multivariate_normal(nominal_state, initial_cov, init_state%samples)
    init_state%mean = nominal_state
    write(*,*) '  Particles generated. Sample (first):', init_state%samples(:, 1)
    write(*,*) ''

    ! Phase 3: MC golden reference
    write(*,*) '>>> Phase 3: MC propagation (golden reference)...'
    write(*,*) '  Propagating ', N_PARTICLES, ' particles for 14 days with RKF45...'

    block
        type(uq_state_type) :: init_mc, final_mc
        integer :: integrator_rkf45

        integrator_rkf45 = 1  ! METHOD_RKF45

        ! Copy initial particles for MC
        call init_mc%allocate_memory(6, N_PARTICLES)
        init_mc%samples = init_state%samples
        init_mc%mean = nominal_state

        call run_particle_propagation( &
            initial_state   = init_mc, &
            reference_orbit = nominal_state, &
            epoch0          = epoch0, &
            t_start         = 0.0_DP, &
            t_end           = DURATION, &
            method_switch   = METHOD_MC, &
            final_state     = final_mc, &
            integrator_switch = integrator_rkf45)

        ! Save MC golden reference
        allocate(mc_golden(6, N_PARTICLES))
        mc_golden = final_mc%samples

        write(*,*) '  MC propagation complete.'
        write(*,*) '  MC final mean: ', final_mc%mean
        write(*,*) ''
    end block

    ! Phase 4: DA scan over orders 3,4,5,6
    write(*,*) '>>> Phase 4: DA order convergence scan...'
    write(*,*) ''

    block
        type(uq_state_type) :: init_da, final_da
        integer :: order, integrator_rkf45
        real(DP) :: max_rel_err(6), rel_err, abs_err
        integer :: i_p, comp, n_exceed
        logical :: all_pass
        character(len=20) :: status_str

        integrator_rkf45 = 1  ! METHOD_RKF45

        write(*,'(A)') ' Order |   X(%)  |   Y(%)  |   Z(%)  |  VX(%)  |  VY(%)  |  VZ(%)  | Exceed | Status'
        write(*,'(A)') '-------+---------+---------+---------+---------+---------+---------+--------+--------'

        do order = 2, DA_MAX_ORDER
            call init_da%allocate_memory(6, N_PARTICLES)
            init_da%samples = init_state%samples
            init_da%mean = nominal_state

            call run_particle_propagation( &
                initial_state   = init_da, &
                reference_orbit = nominal_state, &
                epoch0          = epoch0, &
                t_start         = 0.0_DP, &
                t_end           = DURATION, &
                method_switch   = METHOD_DA, &
                final_state     = final_da, &
                integrator_switch = integrator_rkf45, &
                da_order        = order)

            ! Compute max relative error per component
            max_rel_err = 0.0_DP
            n_exceed = 0
            do i_p = 1, N_PARTICLES
                do comp = 1, 6
                    abs_err = abs(final_da%samples(comp, i_p) - mc_golden(comp, i_p))
                    if (abs(mc_golden(comp, i_p)) > 1.0e-30_DP) then
                        rel_err = abs_err / abs(mc_golden(comp, i_p))
                        max_rel_err(comp) = max(max_rel_err(comp), rel_err)
                        if (rel_err > 0.01_DP) n_exceed = n_exceed + 1
                    end if
                end do
            end do

            all_pass = all(max_rel_err < 0.01_DP)
            if (all_pass) then
                status_str = 'PASS'
            else
                status_str = 'FAIL'
            end if

            write(*,'(I4,A,6(F7.3,A),I6,A,A)') &
                order, '   |', &
                max_rel_err(1)*100, ' |', &
                max_rel_err(2)*100, ' |', &
                max_rel_err(3)*100, ' |', &
                max_rel_err(4)*100, ' |', &
                max_rel_err(5)*100, ' |', &
                max_rel_err(6)*100, ' |', &
                n_exceed, ' | ', trim(status_str)
        end do

        ! Save highest DA order results
        allocate(da_best(6, N_PARTICLES))
        da_best = final_da%samples
    end block

    ! Save MC and DA best results to CSV
    call write_csv('mc_golden.csv', mc_golden, N_PARTICLES)
    call write_csv('da_best.csv', da_best, N_PARTICLES)
    write(*,*) '  Results saved to mc_golden.csv and da_best.csv'
    write(*,*) ''

    ! Phase 5: Summary
    write(*,*) ''
    write(*,*) '============================================================'
    write(*,*) '  Test complete.'
    write(*,*) '============================================================'

contains

    subroutine write_csv(filename, data, n_particles)
        character(len=*), intent(in) :: filename
        real(DP), intent(in) :: data(6, n_particles)
        integer, intent(in) :: n_particles
        integer :: u, i_p

        open(newunit=u, file=filename, status='replace', action='write')
        write(u, '(A)') 'x,y,z,vx,vy,vz'
        do i_p = 1, n_particles
            write(u, '(*(ES22.14, :, ","))') data(:, i_p)
        end do
        close(u)
    end subroutine write_csv

end program test_da_convergence
