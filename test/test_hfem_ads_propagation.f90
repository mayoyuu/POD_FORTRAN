program test_hfem_ads_propagation
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use pod_global, only: DP
    use pod_engine_module, only: pod_engine_init
    use pod_config, only: config
    use pod_data_format_module, only: load_initial_opm
    use pod_dace_classes, only: AlgebraicVector, da_var, dace_initialize, &
        operator(+), operator(*), operator(/), assignment(=)
    use pod_ads_split_module, only: patch_type, manifold_type, sh_count, &
        sh_center, sh_width, sh_contain, sh_map_point, patch_init, &
        patch_destroy, patch_get_trunc_err, patch_get_split_dir, patch_split, &
        mf_init, mf_destroy, mf_push, mf_pop_front
    use pod_da_integrator_module, only: da_adaptive_step_integrate, METHOD_RKF78
    use pod_da_force_model_module, only: set_propagation_epoch, &
        set_srp_scale_uncertainty, clear_srp_scale_uncertainty, &
        cleanup_gravity_network
    implicit none

    integer, parameter :: STATE_DIM = 6
    integer, parameter :: DA_ORDER = 4
    integer, parameter :: MAX_SPLIT_DEPTH = 8
    real(DP), parameter :: PROPAGATION_SECONDS = 604800.0_DP
    real(DP), parameter :: ETA_MEAN = 0.0_DP
    real(DP), parameter :: ETA_HALF_WIDTH = 0.1_DP
    character(len=*), parameter :: CONFIG_FILE = 'config/config.txt'
    character(len=*), parameter :: OPM_FILE = &
        'OPM/L1Halo-1/L1Halo-1_init.opm.json'

    type :: ads_run_stats
        integer :: iterations = 0
        integer :: depth_limited_patches = 0
        integer :: max_queue_size = 0
        integer :: split_counts(7) = 0
    end type ads_run_stats

    real(DP) :: epoch0, state0(STATE_DIM), covariance(STATE_DIM, STATE_DIM)
    real(DP) :: domain_half_width(STATE_DIM), error_tolerance(STATE_DIM)
    integer :: n_fail

    n_fail = 0
    domain_half_width(1:3) = 300.0_DP
    domain_half_width(4:6) = 9.0e-4_DP
    error_tolerance(1:3) = 1.0e-1_DP
    error_tolerance(4:6) = 1.0e-6_DP

    call pod_engine_init(CONFIG_FILE)
    call load_initial_opm(OPM_FILE, epoch0, state0, covariance)
    call set_propagation_epoch(epoch0)
    config%use_srp = .true.

    write(*,'(A)') '============================================================'
    write(*,'(A)') ' HFEM ADS propagation: L1Halo-1, seven days, RKF78'
    write(*,'(A)') '============================================================'

    call run_ads_case('orbit-6D', state0, domain_half_width, &
        error_tolerance, epoch0, 6, .false., n_fail)
    call run_ads_case('orbit+SRP-7D', state0, domain_half_width, &
        error_tolerance, epoch0, 7, .true., n_fail)

    call clear_srp_scale_uncertainty()
    call cleanup_gravity_network()

    if (n_fail /= 0) then
        write(*,'(A,I0)') 'FAIL: HFEM ADS checks failed: ', n_fail
        stop 1
    end if

    write(*,'(A)') 'PASS: HFEM ADS completed orbit-6D and orbit+SRP-7D cases'

contains

    subroutine run_ads_case(label, initial_state, half_width, err_toll, &
                            epoch, n_variables, use_srp_uncertainty, failures)
        character(len=*), intent(in) :: label
        real(DP), intent(in) :: initial_state(STATE_DIM)
        real(DP), intent(in) :: half_width(STATE_DIM), err_toll(STATE_DIM)
        real(DP), intent(in) :: epoch
        integer, intent(in) :: n_variables
        logical, intent(in) :: use_srp_uncertainty
        integer, intent(inout) :: failures

        type(manifold_type) :: result_domain
        type(ads_run_stats) :: stats
        real(DP), allocatable :: point(:)
        real(DP) :: center_state(STATE_DIM), low_eta_state(STATE_DIM)
        real(DP) :: high_eta_state(STATE_DIM), interior_state(STATE_DIM)
        real(DP) :: elapsed, srp_sensitivity
        integer :: count_start, count_end, count_rate, i
        logical :: found

        write(*,'(/,A,A)') '--- ', trim(label)
        write(*,'(A,I0)') 'DA independent variables : ', n_variables
        write(*,'(A,I0)') 'DA order                 : ', DA_ORDER
        write(*,'(A,3(ES12.4,1X))') 'position half-width km   : ', half_width(1:3)
        write(*,'(A,3(ES12.4,1X))') 'velocity half-width km/s : ', half_width(4:6)
        if (use_srp_uncertainty) then
            write(*,'(A,ES12.4)') 'eta_srp half-width       : ', ETA_HALF_WIDTH
        end if

        call system_clock(count_start, count_rate)
        call build_hfem_ads_domain(initial_state, half_width, err_toll, epoch, &
            n_variables, use_srp_uncertainty, result_domain, stats, failures)
        call system_clock(count_end)
        elapsed = real(count_end - count_start, DP) / real(count_rate, DP)

        call assert_true(result_domain%n_patches > 0, &
            trim(label)//' accepts at least one Patch', failures)
        call assert_true(stats%depth_limited_patches == 0, &
            trim(label)//' satisfies truncation tolerance before the depth cap', failures)
        write(*,'(A,I0)') 'accepted patches         : ', result_domain%n_patches
        write(*,'(A,I0)') 'BFS iterations           : ', stats%iterations
        write(*,'(A,I0)') 'maximum queue size       : ', stats%max_queue_size
        write(*,'(A,F12.3,A)') 'domain construction wall: ', elapsed, ' s'
        write(*,'(A,I0)') 'depth-limited patches    : ', stats%depth_limited_patches
        do i = 1, n_variables
            write(*,'(A,I0,A,I0)') 'split count variable ', i, ': ', stats%split_counts(i)
        end do

        allocate(point(n_variables))
        point = 0.0_DP
        call evaluate_ads_point(result_domain, point, center_state, found)
        call assert_true(found, trim(label)//' center is covered', failures)
        call assert_finite_state(center_state, trim(label)//' center state', failures)
        write(*,'(A,6(ES15.7,1X))') 'final center state       : ', center_state

        point = 0.0_DP
        point(1:min(STATE_DIM, n_variables)) = &
            [0.25_DP, -0.25_DP, 0.50_DP, -0.50_DP, 0.125_DP, -0.125_DP]
        call evaluate_ads_point(result_domain, point, interior_state, found)
        call assert_true(found, trim(label)//' interior point is covered', failures)
        call assert_finite_state(interior_state, trim(label)//' interior state', failures)

        if (use_srp_uncertainty) then
            point = 0.0_DP
            point(7) = -1.0_DP
            call evaluate_ads_point(result_domain, point, low_eta_state, found)
            call assert_true(found, trim(label)//' eta=-0.1 endpoint is covered', failures)
            call assert_finite_state(low_eta_state, trim(label)//' eta=-0.1 state', failures)

            point(7) = 1.0_DP
            call evaluate_ads_point(result_domain, point, high_eta_state, found)
            call assert_true(found, trim(label)//' eta=+0.1 endpoint is covered', failures)
            call assert_finite_state(high_eta_state, trim(label)//' eta=+0.1 state', failures)

            srp_sensitivity = maxval(abs(high_eta_state - low_eta_state))
            write(*,'(A,ES15.7)') 'max endpoint SRP response: ', srp_sensitivity
            call assert_true(srp_sensitivity > 1.0e-12_DP, &
                trim(label)//' has resolvable SRP sensitivity', failures)
        end if

        if (allocated(point)) deallocate(point)
        call mf_destroy(result_domain)
        call clear_srp_scale_uncertainty()
    end subroutine run_ads_case

    subroutine build_hfem_ads_domain(initial_state, half_width, err_toll, epoch, &
                                     n_variables, use_srp_uncertainty, &
                                     result_domain, stats, failures)
        real(DP), intent(in) :: initial_state(STATE_DIM)
        real(DP), intent(in) :: half_width(STATE_DIM), err_toll(STATE_DIM)
        real(DP), intent(in) :: epoch
        integer, intent(in) :: n_variables
        logical, intent(in) :: use_srp_uncertainty
        type(manifold_type), intent(out) :: result_domain
        type(ads_run_stats), intent(out) :: stats
        integer, intent(inout) :: failures

        type(manifold_type) :: queue
        type(patch_type) :: input_patch, propagated_patch, left, right
        type(AlgebraicVector) :: initial_da, final_da_nondim, final_da_physical
        real(DP), allocatable :: times(:), nominal_states(:,:)
        real(DP), allocatable :: patch_center(:), patch_width(:)
        real(DP), allocatable :: check_point(:), initial_eval(:)
        real(DP) :: errors(STATE_DIM), excess(STATE_DIM)
        real(DP) :: expected_initial(STATE_DIM), t_end_nondim
        integer :: i, n_steps, component(1), direction

        call dace_initialize(DA_ORDER, n_variables)
        call set_propagation_epoch(epoch)
        call mf_init(queue)
        call mf_init(result_domain)

        call initial_da%init(STATE_DIM)
        do i = 1, 3
            initial_da%elements(i) = &
                (initial_state(i) + half_width(i) * da_var(i)) / config%LU
            initial_da%elements(i+3) = &
                (initial_state(i+3) + half_width(i+3) * da_var(i+3)) / config%VU
        end do

        ! Verify the physical interval map before Patch ownership is transferred.
        allocate(check_point(n_variables))
        check_point = 0.0_DP
        check_point(1:6) = &
            [0.25_DP, -0.50_DP, 0.75_DP, -0.25_DP, 0.50_DP, -0.75_DP]
        if (n_variables >= 7) check_point(7) = 0.5_DP
        initial_eval = initial_da%eval(check_point)
        initial_eval(1:3) = initial_eval(1:3) * config%LU
        initial_eval(4:6) = initial_eval(4:6) * config%VU
        expected_initial = initial_state + half_width * check_point(1:6)
        call assert_true(maxval(abs(initial_eval - expected_initial)) < 1.0e-10_DP, &
            'normalized ADS coordinates recover the physical initial interval', failures)
        deallocate(check_point, initial_eval)

        call patch_init(input_patch, initial_da)
        call initial_da%destroy()
        call mf_push(queue, input_patch)

        stats%iterations = 0
        stats%max_queue_size = 0
        stats%split_counts = 0
        stats%depth_limited_patches = 0
        t_end_nondim = PROPAGATION_SECONDS / config%TU

        do while (queue%n_patches > 0)
            stats%iterations = stats%iterations + 1
            stats%max_queue_size = max(stats%max_queue_size, queue%n_patches)
            call mf_pop_front(queue, input_patch)

            if (use_srp_uncertainty) then
                patch_center = sh_center(input_patch%history)
                patch_width = sh_width(input_patch%history)
                call set_srp_scale_uncertainty(7, &
                    ETA_MEAN + ETA_HALF_WIDTH * patch_center(7), &
                    ETA_HALF_WIDTH * patch_width(7) / 2.0_DP)
            else
                call clear_srp_scale_uncertainty()
            end if

            call da_adaptive_step_integrate(input_patch%da_vec, 0.0_DP, &
                t_end_nondim, METHOD_RKF78, times, nominal_states, &
                final_da_nondim, n_steps)

            call final_da_physical%init(STATE_DIM)
            do i = 1, 3
                final_da_physical%elements(i) = &
                    final_da_nondim%elements(i) * config%LU
                final_da_physical%elements(i+3) = &
                    final_da_nondim%elements(i+3) * config%VU
            end do
            call final_da_nondim%destroy()
            if (allocated(times)) deallocate(times)
            if (allocated(nominal_states)) deallocate(nominal_states)

            call patch_init(propagated_patch, final_da_physical, input_patch%history)
            call final_da_physical%destroy()
            call patch_get_trunc_err(propagated_patch, DA_ORDER, errors)
            excess = max(0.0_DP, errors - err_toll)

            if (maxval(excess) <= 0.0_DP) then
                call mf_push(result_domain, propagated_patch)
                call patch_destroy(input_patch)
            else if (sh_count(input_patch%history, 0) >= MAX_SPLIT_DEPTH) then
                stats%depth_limited_patches = stats%depth_limited_patches + 1
                call mf_push(result_domain, propagated_patch)
                call patch_destroy(input_patch)
            else
                component = maxloc(excess)
                direction = patch_get_split_dir(propagated_patch, component(1), DA_ORDER)
                stats%split_counts(direction) = stats%split_counts(direction) + 1
                call patch_split(input_patch, direction, left, right)
                call mf_push(queue, left)
                call mf_push(queue, right)
                call patch_destroy(input_patch)
                call patch_destroy(propagated_patch)
            end if

            if (allocated(patch_center)) deallocate(patch_center)
            if (allocated(patch_width)) deallocate(patch_width)
        end do

        call clear_srp_scale_uncertainty()
        call mf_destroy(queue)
    end subroutine build_hfem_ads_domain

    subroutine evaluate_ads_point(domain, point_unit, state_out, found)
        type(manifold_type), intent(in) :: domain
        real(DP), intent(in) :: point_unit(:)
        real(DP), intent(out) :: state_out(STATE_DIM)
        logical, intent(out) :: found

        real(DP), allocatable :: local_point(:), values(:)
        integer :: i

        found = .false.
        state_out = 0.0_DP
        do i = 1, domain%n_patches
            if (sh_contain(domain%patches(i)%history, point_unit)) then
                allocate(local_point(size(point_unit)))
                local_point = point_unit
                call sh_map_point(domain%patches(i)%history, local_point)
                values = domain%patches(i)%da_vec%eval(local_point)
                state_out = values(1:STATE_DIM)
                found = .true.
                deallocate(local_point, values)
                return
            end if
        end do
    end subroutine evaluate_ads_point

    subroutine assert_finite_state(state, label, failures)
        real(DP), intent(in) :: state(STATE_DIM)
        character(len=*), intent(in) :: label
        integer, intent(inout) :: failures

        call assert_true(all(ieee_is_finite(state)), label//' is finite', failures)
    end subroutine assert_finite_state

    subroutine assert_true(condition, label, failures)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: label
        integer, intent(inout) :: failures

        if (.not. condition) then
            write(*,'(A,A)') 'FAIL: ', trim(label)
            failures = failures + 1
        end if
    end subroutine assert_true

end program test_hfem_ads_propagation
