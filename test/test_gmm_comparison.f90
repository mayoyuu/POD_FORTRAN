program test_gmm_comparison
    use pod_global, only: DP
    use pod_engine_module, only: pod_engine_init
    use pod_dace_classes, only: dace_initialize
    use pod_data_format_module, only: load_initial_opm
    use pod_uq_state_module, only: uq_state_type
    use pod_random_module, only: init_random_seed, generate_multivariate_normal
    use pod_uq_propagation, only: run_particle_propagation, METHOD_MC, METHOD_DA
    use pod_uq_gmm_state_module, only: uq_gmm_state_type
    use pod_gmm_math_module, only: fit_gmm_to_particles

    implicit none

    character(len=*), parameter :: CONFIG_FILE = 'config/dummy_test_config.txt'
    character(len=*), parameter :: OPM_FILE    = 'OPM/L1Halo_0noise/L1Halo_0noise_init.opm.json'
    integer, parameter :: DA_ORDER     = 5
    integer, parameter :: N_PARTICLES  = 5000
    integer, parameter :: EM_MAX_ITER  = 50
    real(DP), parameter :: EM_TOL      = 1.0e-4_DP
    integer, parameter :: COMP_LIST(3) = [1, 3, 5]
    character(len=*), parameter :: OUT_DIR = 'OPM/L1Halo_0noise/gmm_analysis/'

    real(DP) :: nominal_state(6), initial_cov(6,6), epoch0
    type(uq_state_type) :: init_state, final_mc, final_da
    type(uq_gmm_state_type) :: gmm_mc, gmm_da
    real(DP), allocatable :: mc_particles(:,:), da_particles(:,:)
    real(DP) :: global_mean(6), global_cov(6,6)
    integer :: ic, n_comp

    write(*,*) '============================================================'
    write(*,*) '  GMM Component Comparison Test - L1Halo / DA5 / MC ref'
    write(*,*) '============================================================'

    ! Phase 0: System init
    write(*,*) '>>> Phase 0: System initialization...'
    call pod_engine_init(CONFIG_FILE)
    call dace_initialize(DA_ORDER, 6)

    ! Phase 1: Read OPM
    write(*,*) '>>> Phase 1: Reading OPM from ', OPM_FILE
    call load_initial_opm(OPM_FILE, epoch0, nominal_state, initial_cov)
    write(*,*) '  Epoch (TDB): ', epoch0
    write(*,*) '  Nominal:     ', nominal_state

    ! Phase 2: Generate initial particles
    write(*,*) '>>> Phase 2: Generating ', N_PARTICLES, ' initial particles...'
    call init_state%allocate_memory(6, N_PARTICLES)
    call init_random_seed(.true.)
    call generate_multivariate_normal(nominal_state, initial_cov, init_state%samples)
    init_state%mean = nominal_state

    ! Compute propagation duration: from OPM epoch to first OBS
    block
        real(DP) :: duration
        duration = 821361789.2802_DP - epoch0  ! first OBS ET - OPM epoch ET
        write(*,*) '  Propagation duration: ', duration / 86400.0_DP, ' days'

        ! ============================================================
        ! Phase 3: MC propagation (golden reference)
        ! ============================================================
        write(*,*) '>>> Phase 3: MC propagation (', N_PARTICLES, ' particles)...'
        block
            type(uq_state_type) :: init_mc
            call init_mc%allocate_memory(6, N_PARTICLES)
            init_mc%samples = init_state%samples
            init_mc%mean = nominal_state

            call run_particle_propagation( &
                initial_state   = init_mc, &
                reference_orbit = nominal_state, &
                epoch0          = epoch0, &
                t_start         = 0.0_DP, &
                t_end           = duration, &
                method_switch   = METHOD_MC, &
                final_state     = final_mc)

            allocate(mc_particles(6, N_PARTICLES))
            mc_particles = final_mc%samples
            write(*,*) '  MC final mean: ', final_mc%mean
        end block

        ! ============================================================
        ! Phase 4: DA5 propagation
        ! ============================================================
        write(*,*) '>>> Phase 4: DA', DA_ORDER, ' propagation (', N_PARTICLES, ' particles)...'
        block
            type(uq_state_type) :: init_da
            call init_da%allocate_memory(6, N_PARTICLES)
            init_da%samples = init_state%samples
            init_da%mean = nominal_state

            call run_particle_propagation( &
                initial_state   = init_da, &
                reference_orbit = nominal_state, &
                epoch0          = epoch0, &
                t_start         = 0.0_DP, &
                t_end           = duration, &
                method_switch   = METHOD_DA, &
                final_state     = final_da, &
                da_order        = DA_ORDER)

            allocate(da_particles(6, N_PARTICLES))
            da_particles = final_da%samples
            write(*,*) '  DA final mean: ', final_da%mean
        end block
    end block

    ! ============================================================
    ! Phase 5: Save particles to CSV
    ! ============================================================
    write(*,*) '>>> Phase 5: Saving particle CSVs...'
    call write_csv(OUT_DIR // 'mc_particles.csv', mc_particles, N_PARTICLES)
    call write_csv(OUT_DIR // 'da5_particles.csv', da_particles, N_PARTICLES)

    ! ============================================================
    ! Phase 6: GMM fitting and JSON output
    ! ============================================================
    write(*,*) '>>> Phase 6: GMM fitting (1/3/5 components)...'

    do ic = 1, size(COMP_LIST)
        n_comp = COMP_LIST(ic)
        write(*,*) '  --- n_comp = ', n_comp, ' ---'

        ! --- MC particles → GMM ---
        write(*,*) '    MC GMM...'
        call gmm_mc%allocate_components(n_comp, 6)
        call fit_gmm_to_particles(mc_particles, gmm_mc, EM_MAX_ITER, EM_TOL)
        call compute_gmm_global(gmm_mc, global_mean, global_cov)
        write(*,*) '      MC GMM global mean: ', global_mean(1:3)
        call write_gmm_json(OUT_DIR // 'mc_gmm_comp' // int2str(n_comp) // '.json', &
            'MC', gmm_mc, global_mean, global_cov)

        ! --- DA5 particles → GMM ---
        write(*,*) '    DA5 GMM...'
        call gmm_da%allocate_components(n_comp, 6)
        call fit_gmm_to_particles(da_particles, gmm_da, EM_MAX_ITER, EM_TOL)
        call compute_gmm_global(gmm_da, global_mean, global_cov)
        write(*,*) '      DA5 GMM global mean: ', global_mean(1:3)
        call write_gmm_json(OUT_DIR // 'da5_gmm_comp' // int2str(n_comp) // '.json', &
            'DA5', gmm_da, global_mean, global_cov)
    end do

    write(*,*) ''
    write(*,*) '============================================================'
    write(*,*) '  Test complete. Output in ', OUT_DIR
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
        write(*,*) '    Saved: ', trim(filename)
    end subroutine write_csv

    subroutine compute_gmm_global(gmm, global_mean, global_cov)
        type(uq_gmm_state_type), intent(in) :: gmm
        real(DP), intent(out) :: global_mean(6), global_cov(6,6)
        integer :: i
        real(DP) :: diff(6)

        global_mean = 0.0_DP
        do i = 1, gmm%n_components
            global_mean = global_mean + gmm%components(i)%weight * gmm%components(i)%mean
        end do

        global_cov = 0.0_DP
        do i = 1, gmm%n_components
            diff = gmm%components(i)%mean - global_mean
            global_cov = global_cov + gmm%components(i)%weight * &
                (gmm%components(i)%cov + &
                 matmul(reshape(diff, (/6, 1/)), reshape(diff, (/1, 6/))))
        end do
    end subroutine compute_gmm_global

    subroutine write_gmm_json(filename, method, gmm, global_mean, global_cov)
        character(len=*), intent(in) :: filename, method
        type(uq_gmm_state_type), intent(in) :: gmm
        real(DP), intent(in) :: global_mean(6), global_cov(6,6)
        integer :: u, i, j, k

        open(newunit=u, file=filename, status='replace', action='write')
        write(u, '(A)') '{'
        write(u, '(A,A,A)') '  "method": "', trim(method), '",'
        write(u, '(A,I0,A)') '  "n_components": ', gmm%n_components, ','

        ! global_mean
        write(u, '(A)', advance='no') '  "global_mean": ['
        do j = 1, 6
            if (j < 6) then
                write(u, '(ES22.15,A)', advance='no') global_mean(j), ', '
            else
                write(u, '(ES22.15,A)') global_mean(j), '],'
            end if
        end do

        ! global_cov as 6x6 2D row-major array
        write(u, '(A)') '  "global_cov": ['
        do j = 1, 6
            write(u, '(A)', advance='no') '    ['
            do k = 1, 6
                if (k < 6) then
                    write(u, '(ES22.15,A)', advance='no') global_cov(j,k), ', '
                else
                    write(u, '(ES22.15)', advance='no') global_cov(j,k)
                end if
            end do
            if (j < 6) then
                write(u, '(A)') '],'
            else
                write(u, '(A)') ']'
            end if
        end do
        write(u, '(A)') '  ],'

        ! components
        write(u, '(A)') '  "components": ['
        do i = 1, gmm%n_components
            write(u, '(A)') '    {'
            write(u, '(A,ES22.15,A)') '      "weight": ', gmm%components(i)%weight, ','

            write(u, '(A)', advance='no') '      "mean": ['
            do j = 1, 6
                if (j < 6) then
                    write(u, '(ES22.15,A)', advance='no') gmm%components(i)%mean(j), ', '
                else
                    write(u, '(ES22.15,A)') gmm%components(i)%mean(j), '],'
                end if
            end do

            ! cov as 6x6 2D row-major array
            write(u, '(A)') '      "cov": ['
            do j = 1, 6
                write(u, '(A)', advance='no') '        ['
                do k = 1, 6
                    if (k < 6) then
                        write(u, '(ES22.15,A)', advance='no') gmm%components(i)%cov(j,k), ', '
                    else
                        write(u, '(ES22.15)', advance='no') gmm%components(i)%cov(j,k)
                    end if
                end do
                if (j < 6) then
                    write(u, '(A)') '],'
                else
                    write(u, '(A)') ']'
                end if
            end do

            if (i < gmm%n_components) then
                write(u, '(A)') '      ]'
                write(u, '(A)') '    },'
            else
                write(u, '(A)') '      ]'
                write(u, '(A)') '    }'
            end if
        end do
        write(u, '(A)') '  ]'
        write(u, '(A)') '}'
        close(u)
        write(*,*) '    Saved: ', trim(filename)
    end subroutine write_gmm_json

    function int2str(n) result(s)
        integer, intent(in) :: n
        character(len=2) :: s
        write(s, '(I0)') n
    end function int2str

end program test_gmm_comparison
