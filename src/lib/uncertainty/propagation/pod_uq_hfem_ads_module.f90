!> HFEM uncertainty propagation through an automatically split DA manifold.
module pod_uq_hfem_ads_module
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use pod_global, only: DP
    use pod_config, only: config
    use pod_dace_classes, only: AlgebraicVector, da_var, dace_initialize, &
        operator(+), operator(*), operator(/), assignment(=)
    use pod_ads_split_module, only: patch_type, manifold_type, sh_count, &
        sh_center, sh_width, patch_init, patch_destroy, patch_get_trunc_err, &
        patch_get_split_dir, patch_split, mf_init, mf_destroy, mf_push, &
        mf_pop_front, mf_evaluate_points
    use pod_uq_ads_coordinates_module, only: ADS_COORD_COMPONENT, &
        ads_coordinate_map_type, ads_build_coordinate_map, &
        ads_physical_to_unit,ads_physical_parameters_type,ads_patch_mode_parameters, &
        ADS_PARAM_AUTO,ADS_PARAM_NONE,ADS_PARAM_SRP_ETA,ADS_PARAM_SRP_MODES
    use pod_da_integrator_module, only: da_adaptive_step_integrate, METHOD_RKF78
    use pod_da_force_model_module, only: set_propagation_epoch, &
        set_srp_scale_uncertainty, clear_srp_scale_uncertainty, &
        set_srp_mode_uncertainty,clear_srp_mode_uncertainty
    use pod_srp_mode_module, only: interpolate_srp_modes
    implicit none
    private

    integer, parameter :: STATE_DIM = 6

    type, public :: hfem_ads_options_type
        integer :: coordinate_mode = ADS_COORD_COMPONENT
        integer :: parameterization=ADS_PARAM_AUTO
        real(DP) :: domain_sigma = 3.0_DP
        real(DP) :: srp_sigma = 0.0_DP
        integer :: da_storage_variables = 0
        integer :: da_order = 4
        integer :: max_split_depth = 8
        real(DP) :: error_tolerance(6) = [1.0e-1_DP, 1.0e-1_DP, &
            1.0e-1_DP, 1.0e-6_DP, 1.0e-6_DP, 1.0e-6_DP]
        real(DP) :: rel_tol = 1.0e-12_DP
        real(DP) :: abs_tol = 1.0e-12_DP
        real(DP) :: dt_min = 1.0e-6_DP
        real(DP) :: dt_max = 3600.0_DP
        integer :: max_steps = 100000
    end type hfem_ads_options_type

    type, public :: hfem_ads_split_event_type
        integer :: bfs_iteration = 0
        integer :: parent_depth = -1
        integer :: component = 0
        integer :: direction = 0
        real(DP) :: errors(6) = 0.0_DP
        real(DP) :: tolerances(6) = 0.0_DP
        real(DP) :: excess(6) = 0.0_DP
        real(DP) :: direction_scores(8) = 0.0_DP
        real(DP) :: parent_center(8) = 0.0_DP
        real(DP) :: parent_width(8) = 0.0_DP
    end type hfem_ads_split_event_type

    type, public :: hfem_ads_stats_type
        integer :: requested_count = 0
        integer :: sampled_count = 0
        integer :: rejected_count = 0
        integer :: input_count = 0
        integer :: inside_count = 0
        integer :: outside_count = 0
        integer :: propagated_count = 0
        integer :: written_count = 0
        integer :: n_patches = 0
        integer :: bfs_iterations = 0
        integer :: max_queue_size = 0
        integer :: depth_limited_patches = 0
        integer :: coordinate_mode = ADS_COORD_COMPONENT
        integer :: n_variables = 6
        integer :: parameterization=ADS_PARAM_NONE
        real(DP) :: domain_sigma = 3.0_DP
        real(DP) :: srp_sigma = 0.0_DP
        real(DP) :: basis6(6,6) = 0.0_DP
        real(DP) :: elapsed_seconds = 0.0_DP
        logical :: split_occurred = .false.
        integer :: first_split_component = 0
        integer :: first_split_direction = 0
        integer :: first_split_depth = -1
        integer :: first_split_bfs_iteration = 0
        real(DP) :: first_split_errors(6) = 0.0_DP
        real(DP) :: first_split_tolerances(6) = 0.0_DP
        real(DP) :: first_split_excess(6) = 0.0_DP
        real(DP) :: first_split_direction_scores(8) = 0.0_DP
        real(DP) :: first_split_parent_center(8) = 0.0_DP
        real(DP) :: first_split_parent_width(8) = 0.0_DP
        integer, allocatable :: split_counts(:)
        type(hfem_ads_split_event_type), allocatable :: split_events(:)
    end type hfem_ads_stats_type

    ! Incremental ADS domain. Each Patch maps the original unit box to the
    ! physical Cartesian state at current_time; its splitting history is kept.
    type, public :: hfem_ads_history_type
        type(ads_coordinate_map_type) :: coordinate_map
        type(hfem_ads_options_type) :: options
        type(hfem_ads_stats_type) :: stats
        type(manifold_type) :: domain
        type(manifold_type) :: previous_domain
        real(DP) :: previous_time = 0.0_DP
        logical :: previous_available = .false.
        real(DP) :: epoch0 = 0.0_DP
        real(DP) :: current_time = 0.0_DP
        integer :: da_variables = 0
        logical :: initialized = .false.
    end type hfem_ads_history_type

    public :: hfem_ads_propagate
    public :: hfem_ads_history_init, hfem_ads_history_advance
    public :: hfem_ads_history_evaluate, hfem_ads_history_destroy

contains

    subroutine hfem_ads_propagate(nominal_state, covariance, epoch0, t_start, &
            t_end, input_samples, options, output_samples, stats, status, message)
        real(DP), intent(in) :: nominal_state(6), covariance(6,6)
        real(DP), intent(in) :: epoch0, t_start, t_end
        real(DP), intent(in) :: input_samples(:,:)
        type(hfem_ads_options_type), intent(in) :: options
        real(DP), allocatable, intent(out) :: output_samples(:,:)
        type(hfem_ads_stats_type), intent(out) :: stats
        integer, intent(out) :: status
        character(len=*), intent(out) :: message
        type(ads_coordinate_map_type) :: coordinate_map
        type(manifold_type) :: domain
        real(DP), allocatable :: unit_all(:,:), unit_inside(:,:), propagated(:,:)
        real(DP), allocatable :: unit_point(:)
        logical, allocatable :: is_inside(:), found(:)
        integer, allocatable :: inside_indices(:)
        integer :: i, j, n_inside, n_output, map_status, count_start, count_end, count_rate
        real(DP) :: unit_tolerance
        type(ads_physical_parameters_type) :: parameters

        status = 0
        message = ''
        call clear_srp_scale_uncertainty()
        call clear_srp_mode_uncertainty()
        stats = hfem_ads_stats_type()
        stats%requested_count = size(input_samples,2)
        stats%input_count = size(input_samples,2)
        if (.not. all(ieee_is_finite(input_samples))) then
            call fail(-1, 'HFEM ADS input samples must be finite')
            return
        end if
        if (.not. ieee_is_finite(epoch0) .or. .not. ieee_is_finite(t_start) .or. &
            .not. ieee_is_finite(t_end) .or. t_end < t_start) then
            call fail(-2, 'HFEM ADS end time must not precede start time')
            return
        end if

        call ads_build_coordinate_map(nominal_state, covariance, &
            options%coordinate_mode, options%domain_sigma, options%srp_sigma, &
            coordinate_map, map_status, message,options%parameterization)
        if (map_status /= 0) then
            status = map_status
            return
        end if
        stats%coordinate_mode = coordinate_map%mode
        stats%n_variables = coordinate_map%n_variables
        stats%parameterization=coordinate_map%parameterization
        stats%domain_sigma = coordinate_map%domain_sigma
        stats%srp_sigma = coordinate_map%srp_sigma
        stats%basis6 = coordinate_map%basis6
        call mode_preflight(coordinate_map,epoch0+t_start,epoch0+t_end,status,message)
        if(status/=0) return
        allocate(stats%split_counts(coordinate_map%n_variables), source=0)
        if (size(input_samples,1) /= coordinate_map%n_variables) then
            call fail(-3, 'HFEM ADS sample dimension does not match enabled uncertainties')
            return
        end if
        if (options%da_order < 1 .or. options%max_split_depth < 0 .or. &
            options%max_steps < 2 .or. any(options%error_tolerance <= 0.0_DP) .or. &
            .not. all(ieee_is_finite(options%error_tolerance)) .or. &
            .not. ieee_is_finite(options%rel_tol) .or. options%rel_tol <= 0.0_DP .or. &
            .not. ieee_is_finite(options%abs_tol) .or. options%abs_tol <= 0.0_DP .or. &
            .not. ieee_is_finite(options%dt_min) .or. options%dt_min <= 0.0_DP .or. &
            .not. ieee_is_finite(options%dt_max) .or. options%dt_max < options%dt_min) then
            call fail(-4, 'HFEM ADS options contain a non-positive limit')
            return
        end if

        allocate(unit_all(coordinate_map%n_variables, size(input_samples,2)))
        allocate(is_inside(size(input_samples,2)), source=.false.)
        unit_tolerance = 1000.0_DP*epsilon(1.0_DP)
        do j = 1, size(input_samples,2)
            parameters=ads_physical_parameters_type()
            if(coordinate_map%parameterization==ADS_PARAM_SRP_ETA) parameters%eta=input_samples(7,j)
            if(coordinate_map%parameterization==ADS_PARAM_SRP_MODES) then
                parameters%s=input_samples(7,j); parameters%u=input_samples(8,j)
            end if
            call ads_physical_to_unit(coordinate_map, &
                input_samples(1:6,j)-nominal_state, parameters, unit_point, map_status)
            if (map_status == 0) then
                unit_all(:,j) = unit_point
                is_inside(j) = maxval(abs(unit_point)) <= 1.0_DP + unit_tolerance
            end if
        end do
        n_inside = count(is_inside)
        stats%inside_count = n_inside
        stats%outside_count = stats%input_count - n_inside
        if (n_inside == 0) then
            allocate(output_samples(coordinate_map%n_variables,0))
            call fail(-5, 'HFEM ADS has no samples inside its finite root domain')
            return
        end if
        allocate(unit_inside(coordinate_map%n_variables,n_inside))
        allocate(inside_indices(n_inside))
        i = 0
        do j = 1, size(input_samples,2)
            if (.not. is_inside(j)) cycle
            i = i + 1
            unit_inside(:,i) = unit_all(:,j)
            inside_indices(i) = j
        end do

        call system_clock(count_start, count_rate)
        call build_domain(nominal_state, epoch0, t_start, t_end, coordinate_map, &
            options, domain, stats,status,message)
        if(status/=0) then
            call mf_destroy(domain)
            return
        end if
        allocate(propagated(6,n_inside), found(n_inside))
        call mf_evaluate_points(domain, unit_inside, propagated, found, map_status)
        call system_clock(count_end)
        stats%elapsed_seconds = real(count_end-count_start,DP)/real(count_rate,DP)
        if (map_status /= 0) then
            call mf_destroy(domain)
            call clear_srp_scale_uncertainty()
            call fail(-6, 'HFEM ADS manifold evaluation failed')
            return
        end if

        n_output = count(found)
        stats%propagated_count = n_output
        stats%written_count = n_output
        allocate(output_samples(coordinate_map%n_variables,n_output))
        i = 0
        do j = 1, n_inside
            if (.not. found(j)) cycle
            i = i + 1
            output_samples(1:6,i) = propagated(:,j)
            if(coordinate_map%n_variables>6) &
                output_samples(7:coordinate_map%n_variables,i)=input_samples(7:coordinate_map%n_variables,inside_indices(j))
        end do
        call mf_destroy(domain)
        call clear_srp_scale_uncertainty()
        if (n_output == 0) call fail(-7, 'HFEM ADS propagated no valid samples')

    contains
        subroutine fail(code, text)
            integer, intent(in) :: code
            character(len=*), intent(in) :: text
            status = code
            message = text
            call clear_srp_scale_uncertainty()
            call clear_srp_mode_uncertainty()
        end subroutine fail
    end subroutine hfem_ads_propagate

    subroutine build_domain(nominal_state, epoch0, t_start, t_end, coordinate_map, &
            options, result_domain, stats,status,message)
        real(DP), intent(in) :: nominal_state(6), epoch0, t_start, t_end
        type(ads_coordinate_map_type), intent(in) :: coordinate_map
        type(hfem_ads_options_type), intent(in) :: options
        type(manifold_type), intent(out) :: result_domain
        type(hfem_ads_stats_type), intent(inout) :: stats
        integer,intent(out) :: status
        character(len=*),intent(out) :: message
        type(manifold_type) :: queue
        type(patch_type) :: input_patch, propagated_patch, left, right
        type(AlgebraicVector) :: initial_da, final_da_nondim, final_da_physical
        real(DP), allocatable :: times(:), nominal_states(:,:)
        real(DP), allocatable :: patch_center(:), patch_width(:)
        real(DP) :: errors(6), excess(6), state_scale, eta_scale
        real(DP) :: direction_scores(8)
        integer :: i, j, n_steps, component(1), direction

        status=0; message=''
        call dace_initialize(options%da_order, coordinate_map%n_variables)
        call set_propagation_epoch(epoch0)
        call clear_srp_scale_uncertainty()
        call mf_init(queue)
        call mf_init(result_domain)
        call initial_da%init(6)
        do i = 1, 6
            state_scale = merge(config%LU, config%VU, i <= 3)
            initial_da%elements(i) = nominal_state(i)/state_scale
            do j = 1, 6
                initial_da%elements(i) = initial_da%elements(i) + &
                    (coordinate_map%basis6(i,j)/state_scale)*da_var(j)
            end do
        end do
        call patch_init(input_patch, initial_da)
        call mf_push(queue, input_patch)
        call patch_destroy(input_patch)

        do while (queue%n_patches > 0)
            stats%bfs_iterations = stats%bfs_iterations + 1
            stats%max_queue_size = max(stats%max_queue_size, queue%n_patches)
            call mf_pop_front(queue, input_patch)
            call prepare_patch_uncertainty(coordinate_map,input_patch,status,message)
            if(status/=0) exit

            if (t_end == t_start) then
                final_da_nondim = input_patch%da_vec
            else
                call da_adaptive_step_integrate(input_patch%da_vec, &
                    t_start/config%TU, t_end/config%TU, METHOD_RKF78, times, &
                    nominal_states, final_da_nondim, n_steps, options%max_steps, &
                    options%rel_tol, options%abs_tol, options%dt_min, options%dt_max)
                if(.not.integration_complete(times,nominal_states,n_steps,t_end/config%TU)) then
                    status=-10; message='HFEM ADS integration did not reach end time'; exit
                end if
            end if
            call final_da_physical%init(6)
            do i = 1, 3
                final_da_physical%elements(i) = final_da_nondim%elements(i)*config%LU
                final_da_physical%elements(i+3) = &
                    final_da_nondim%elements(i+3)*config%VU
            end do
            call final_da_nondim%destroy()
            if (allocated(times)) deallocate(times)
            if (allocated(nominal_states)) deallocate(nominal_states)

            call patch_init(propagated_patch, final_da_physical, input_patch%history)
            call patch_get_trunc_err(propagated_patch, options%da_order, errors)
            excess = max(0.0_DP, errors-options%error_tolerance)
            if (maxval(excess) <= 0.0_DP) then
                call mf_push(result_domain, propagated_patch)
                call patch_destroy(input_patch)
                call patch_destroy(propagated_patch)
            else if (sh_count(input_patch%history,0) >= options%max_split_depth) then
                stats%depth_limited_patches = stats%depth_limited_patches + 1
                call mf_push(result_domain, propagated_patch)
                call patch_destroy(input_patch)
                call patch_destroy(propagated_patch)
            else
                component = maxloc(excess)
                direction = patch_get_split_dir(propagated_patch, &
                    component(1), options%da_order, direction_scores(1:coordinate_map%n_variables))
                call record_split_event(stats, input_patch, errors, &
                    options%error_tolerance, excess, component(1), direction, &
                    direction_scores)
                stats%split_counts(direction) = stats%split_counts(direction) + 1
                call patch_split(input_patch, direction, left, right)
                call mf_push(queue, left)
                call mf_push(queue, right)
                call patch_destroy(left)
                call patch_destroy(right)
                call patch_destroy(input_patch)
                call patch_destroy(propagated_patch)
            end if
        end do
        stats%n_patches = result_domain%n_patches
        call initial_da%destroy()
        call final_da_nondim%destroy()
        call final_da_physical%destroy()
        call patch_destroy(input_patch)
        call clear_srp_scale_uncertainty()
        call clear_srp_mode_uncertainty()
        call mf_destroy(queue)
    end subroutine build_domain

    subroutine hfem_ads_history_init(history, nominal_state, covariance, epoch0, &
            options, status, message)
        type(hfem_ads_history_type), intent(inout) :: history
        real(DP), intent(in) :: nominal_state(6), covariance(6,6), epoch0
        type(hfem_ads_options_type), intent(in) :: options
        integer, intent(out) :: status
        character(len=*), intent(out) :: message
        type(AlgebraicVector) :: initial_da
        type(patch_type) :: initial_patch
        integer :: i, j, map_status

        status = 0
        message = ''
        if (history%initialized) call hfem_ads_history_destroy(history)
        call clear_srp_scale_uncertainty()
        call clear_srp_mode_uncertainty()
        if (.not. ieee_is_finite(epoch0)) then
            status = -1
            message = 'ADS history epoch must be finite'
            return
        end if
        if (options%da_order < 1 .or. options%max_split_depth < 0 .or. &
            options%max_steps < 2 .or. any(options%error_tolerance <= 0.0_DP) .or. &
            .not. all(ieee_is_finite(options%error_tolerance)) .or. &
            .not. ieee_is_finite(options%rel_tol) .or. options%rel_tol <= 0.0_DP .or. &
            .not. ieee_is_finite(options%abs_tol) .or. options%abs_tol <= 0.0_DP .or. &
            .not. ieee_is_finite(options%dt_min) .or. options%dt_min <= 0.0_DP .or. &
            .not. ieee_is_finite(options%dt_max) .or. options%dt_max < options%dt_min) then
            status = -2
            message = 'ADS history options contain an invalid limit'
            return
        end if
        call ads_build_coordinate_map(nominal_state, covariance, &
            options%coordinate_mode, options%domain_sigma, options%srp_sigma, &
            history%coordinate_map, map_status, message,options%parameterization)
        if (map_status /= 0) then
            status = map_status
            return
        end if

        history%da_variables = history%coordinate_map%n_variables
        call mode_preflight(history%coordinate_map,epoch0,epoch0,status,message)
        if(status/=0) return
        if (options%da_storage_variables > 0) then
            if (options%da_storage_variables < history%coordinate_map%n_variables .or. &
                options%da_storage_variables > 8) then
                status = -3
                message = 'ADS storage dimension must include every active variable'
                return
            end if
            history%da_variables = options%da_storage_variables
        end if

        history%options = options
        history%epoch0 = epoch0
        history%current_time = 0.0_DP
        history%stats = hfem_ads_stats_type()
        history%stats%coordinate_mode = history%coordinate_map%mode
        history%stats%n_variables = history%coordinate_map%n_variables
        history%stats%parameterization=history%coordinate_map%parameterization
        history%stats%domain_sigma = options%domain_sigma
        history%stats%srp_sigma = options%srp_sigma
        history%stats%basis6 = history%coordinate_map%basis6
        allocate(history%stats%split_counts(history%da_variables), source=0)
        call dace_initialize(options%da_order, history%da_variables)
        call set_propagation_epoch(epoch0)
        call clear_srp_scale_uncertainty()
        call initial_da%init(STATE_DIM)
        do i = 1, STATE_DIM
            initial_da%elements(i) = nominal_state(i)
            do j = 1, STATE_DIM
                initial_da%elements(i) = initial_da%elements(i) + &
                    history%coordinate_map%basis6(i,j)*da_var(j)
            end do
        end do
        call patch_init(initial_patch, initial_da)
        call mf_init(history%domain)
        call mf_push(history%domain, initial_patch)
        call patch_destroy(initial_patch)
        history%stats%n_patches = 1
        history%initialized = .true.
    end subroutine hfem_ads_history_init

    subroutine hfem_ads_history_advance(history, next_time, status, message)
        type(hfem_ads_history_type), intent(inout) :: history
        real(DP), intent(in) :: next_time
        integer, intent(out) :: status
        character(len=*), intent(out) :: message
        type(manifold_type) :: queue, accepted
        type(patch_type) :: input_patch, propagated_patch, left, right
        type(AlgebraicVector) :: starting_nd, ending_nd, ending_physical
        real(DP), allocatable :: times(:), nominal_states(:,:), patch_center(:), patch_width(:)
        real(DP) :: errors(6), excess(6), state_scale, eta_scale
        real(DP) :: direction_scores(8)
        integer :: i, component(1), direction, n_steps, count_start, count_end, count_rate

        status = 0
        message = ''
        call clear_srp_scale_uncertainty()
        call clear_srp_mode_uncertainty()
        if (.not. history%initialized) then
            status = -1
            message = 'ADS history is not initialized'
            return
        end if
        if (.not. ieee_is_finite(next_time) .or. next_time < history%current_time) then
            status = -2
            message = 'ADS history time must be finite and nondecreasing'
            return
        end if
        if (next_time == history%current_time) return
        call mode_preflight(history%coordinate_map,history%epoch0+history%current_time, &
            history%epoch0+next_time,status,message)
        if(status/=0) return

        call system_clock(count_start, count_rate)
        call mf_init(queue)
        call mf_init(accepted)
        do i = 1, history%domain%n_patches
            call mf_push(queue, history%domain%patches(i))
        end do
        history%stats%bfs_iterations = 0
        history%stats%max_queue_size = 0
        history%stats%depth_limited_patches = 0
        history%stats%split_counts = 0
        call reset_first_split(history%stats)
        do while (queue%n_patches > 0)
            history%stats%bfs_iterations = history%stats%bfs_iterations + 1
            history%stats%max_queue_size = max(history%stats%max_queue_size, queue%n_patches)
            call mf_pop_front(queue, input_patch)
            call prepare_patch_uncertainty(history%coordinate_map,input_patch,status,message)
            if(status/=0) exit

            call starting_nd%init(STATE_DIM)
            do i = 1, STATE_DIM
                state_scale = merge(config%LU, config%VU, i <= 3)
                starting_nd%elements(i) = input_patch%da_vec%elements(i)/state_scale
            end do
            call da_adaptive_step_integrate(starting_nd, &
                history%current_time/config%TU, next_time/config%TU, METHOD_RKF78, &
                times, nominal_states, ending_nd, n_steps, history%options%max_steps, &
                history%options%rel_tol, history%options%abs_tol, &
                history%options%dt_min, history%options%dt_max)
            if(.not.integration_complete(times,nominal_states,n_steps,next_time/config%TU)) then
                status=-10; message='HFEM ADS history integration did not reach end time'; exit
            end if
            call starting_nd%destroy()
            call ending_physical%init(STATE_DIM)
            do i = 1, STATE_DIM
                state_scale = merge(config%LU, config%VU, i <= 3)
                ending_physical%elements(i) = ending_nd%elements(i)*state_scale
            end do
            call ending_nd%destroy()
            if (allocated(times)) deallocate(times)
            if (allocated(nominal_states)) deallocate(nominal_states)

            call patch_init(propagated_patch, ending_physical, input_patch%history)
            call patch_get_trunc_err(propagated_patch, history%options%da_order, errors)
            excess = max(0.0_DP, errors-history%options%error_tolerance)
            if (maxval(excess) <= 0.0_DP) then
                call mf_push(accepted, propagated_patch)
            else if (sh_count(input_patch%history,0) >= history%options%max_split_depth) then
                history%stats%depth_limited_patches = history%stats%depth_limited_patches + 1
                call mf_push(accepted, propagated_patch)
            else
                component = maxloc(excess)
                direction = patch_get_split_dir(propagated_patch, &
                    component(1), history%options%da_order, direction_scores(1:history%da_variables))
                call record_split_event(history%stats, input_patch, errors, &
                    history%options%error_tolerance, excess, component(1), direction, &
                    direction_scores)
                history%stats%split_counts(direction) = &
                    history%stats%split_counts(direction) + 1
                call patch_split(input_patch, direction, left, right)
                call mf_push(queue, left)
                call mf_push(queue, right)
                call patch_destroy(left)
                call patch_destroy(right)
            end if
            call patch_destroy(input_patch)
            call patch_destroy(propagated_patch)
        end do
        call starting_nd%destroy()
        call ending_nd%destroy()
        call ending_physical%destroy()
        call patch_destroy(input_patch)
        call mf_destroy(queue)
        call clear_srp_scale_uncertainty()
        call clear_srp_mode_uncertainty()
        if(status/=0) then
            call mf_destroy(accepted)
            return
        end if
        if (history%previous_available) call mf_destroy(history%previous_domain)
        call move_alloc(history%domain%patches, history%previous_domain%patches)
        history%previous_domain%n_patches = history%domain%n_patches
        history%previous_time = history%current_time
        history%previous_available = .true.
        call move_alloc(accepted%patches, history%domain%patches)
        history%domain%n_patches = accepted%n_patches
        accepted%n_patches = 0
        history%stats%n_patches = history%domain%n_patches
        history%current_time = next_time
        call clear_srp_scale_uncertainty()
        call system_clock(count_end)
        history%stats%elapsed_seconds = real(count_end-count_start,DP)/real(count_rate,DP)
    end subroutine hfem_ads_history_advance

    subroutine reset_first_split(stats)
        type(hfem_ads_stats_type), intent(inout) :: stats

        stats%split_occurred = .false.
        stats%first_split_component = 0
        stats%first_split_direction = 0
        stats%first_split_depth = -1
        stats%first_split_bfs_iteration = 0
        stats%first_split_errors = 0.0_DP
        stats%first_split_tolerances = 0.0_DP
        stats%first_split_excess = 0.0_DP
        stats%first_split_direction_scores = 0.0_DP
        stats%first_split_parent_center = 0.0_DP
        stats%first_split_parent_width = 0.0_DP
        if (allocated(stats%split_events)) deallocate(stats%split_events)
        allocate(stats%split_events(0))
    end subroutine reset_first_split

    subroutine record_split_event(stats, parent, errors, tolerances, excess, &
            component, direction, direction_scores)
        type(hfem_ads_stats_type), intent(inout) :: stats
        type(patch_type), intent(in) :: parent
        real(DP), intent(in) :: errors(6), tolerances(6), excess(6)
        integer, intent(in) :: component, direction
        real(DP), intent(in) :: direction_scores(:)
        real(DP), allocatable :: center(:), width(:)
        type(hfem_ads_split_event_type), allocatable :: grown(:)
        type(hfem_ads_split_event_type) :: event
        integer :: nvars, n_events

        nvars = min(stats%n_variables, 8)
        center = sh_center(parent%history)
        width = sh_width(parent%history)
        event%bfs_iteration = stats%bfs_iterations
        event%parent_depth = sh_count(parent%history, 0)
        event%component = component
        event%direction = direction
        event%errors = errors
        event%tolerances = tolerances
        event%excess = excess
        event%direction_scores(1:nvars) = direction_scores(1:nvars)
        event%parent_center(1:nvars) = center(1:nvars)
        event%parent_width(1:nvars) = width(1:nvars)

        if (.not. allocated(stats%split_events)) allocate(stats%split_events(0))
        n_events = size(stats%split_events)
        allocate(grown(n_events+1))
        if (n_events > 0) grown(1:n_events) = stats%split_events
        grown(n_events+1) = event
        call move_alloc(grown, stats%split_events)

        if (stats%split_occurred) return
        stats%split_occurred = .true.
        stats%first_split_component = event%component
        stats%first_split_direction = event%direction
        stats%first_split_depth = event%parent_depth
        stats%first_split_bfs_iteration = event%bfs_iteration
        stats%first_split_errors = event%errors
        stats%first_split_tolerances = event%tolerances
        stats%first_split_excess = event%excess
        stats%first_split_direction_scores = event%direction_scores
        stats%first_split_parent_center = event%parent_center
        stats%first_split_parent_width = event%parent_width
    end subroutine record_split_event

    subroutine hfem_ads_history_evaluate(history, unit_points, values, found, status)
        type(hfem_ads_history_type), intent(in) :: history
        real(DP), intent(in) :: unit_points(:,:)
        real(DP), intent(out) :: values(:,:)
        logical, intent(out) :: found(:)
        integer, intent(out) :: status

        real(DP), allocatable :: embedded_points(:,:)

        status = -1
        found = .false.
        values = 0.0_DP
        if (.not. history%initialized) return
        if (size(unit_points,1) == history%da_variables) then
            call mf_evaluate_points(history%domain, unit_points, values, found, status)
        else if (size(unit_points,1) == history%coordinate_map%n_variables .and. &
                 history%da_variables > history%coordinate_map%n_variables) then
            allocate(embedded_points(history%da_variables,size(unit_points,2)),source=0.0_DP)
            embedded_points(1:size(unit_points,1),:) = unit_points
            call mf_evaluate_points(history%domain, embedded_points, values, found, status)
            deallocate(embedded_points)
        else
            status = -2
        end if
    end subroutine hfem_ads_history_evaluate

    subroutine hfem_ads_history_destroy(history)
        type(hfem_ads_history_type), intent(inout) :: history

        call clear_srp_scale_uncertainty()
        call clear_srp_mode_uncertainty()
        if (.not. history%initialized) return
        call mf_destroy(history%domain)
        if (history%previous_available) call mf_destroy(history%previous_domain)
        history%previous_available = .false.
        call clear_srp_scale_uncertainty()
        history%stats = hfem_ads_stats_type()
        history%da_variables = 0
        history%initialized = .false.
    end subroutine hfem_ads_history_destroy

    subroutine mode_preflight(map,start_epoch,end_epoch,status,message)
        type(ads_coordinate_map_type),intent(in) :: map
        real(DP),intent(in) :: start_epoch,end_epoch
        integer,intent(out) :: status
        character(len=*),intent(out) :: message
        real(DP) :: eE(3),eM(3)
        status=0; message=''
        if((map%parameterization==ADS_PARAM_SRP_MODES) .neqv. (trim(config%srp_model)=='time_varying_modes')) then
            status=-8; message='HFEM modes parameterization and force model must agree'; return
        end if
        if(map%parameterization/=ADS_PARAM_SRP_MODES) return
        call interpolate_srp_modes(start_epoch,eE,eM,status,message)
        if(status==0) call interpolate_srp_modes(end_epoch,eE,eM,status,message)
    end subroutine

    subroutine prepare_patch_uncertainty(map,patch,status,message)
        type(ads_coordinate_map_type),intent(in) :: map
        type(patch_type),intent(in) :: patch
        integer,intent(out) :: status
        character(len=*),intent(out) :: message
        real(DP),allocatable :: center(:),width(:)
        real(DP) :: scale,s0,ds,u0,du
        status=0; message=''
        call clear_srp_scale_uncertainty(); call clear_srp_mode_uncertainty()
        if(map%parameterization==ADS_PARAM_NONE) return
        center=sh_center(patch%history); width=sh_width(patch%history)
        select case(map%parameterization)
        case(ADS_PARAM_SRP_ETA)
            scale=map%domain_sigma*map%srp_sigma
            call set_srp_scale_uncertainty(7,scale*center(7),scale*width(7)/2)
        case(ADS_PARAM_SRP_MODES)
            call ads_patch_mode_parameters(map,center(1:8),width(1:8),s0,ds,u0,du,status)
            if(status/=0) then
                message='invalid global s/u patch mapping'; return
            end if
            call set_srp_mode_uncertainty(7,s0,ds,8,u0,du)
        end select
    end subroutine

    logical function integration_complete(times,states,n,end_time) result(ok)
        real(DP),allocatable,intent(in) :: times(:),states(:,:)
        integer,intent(in) :: n
        real(DP),intent(in) :: end_time
        ok=.false.
        if(.not.allocated(times) .or. .not.allocated(states) .or. n<1) return
        if(n>size(times) .or. size(states,1)/=6 .or. n>size(states,2)) return
        ok=all(ieee_is_finite(states(:,n))) .and. &
            abs(times(n)-end_time)<=1000*epsilon(1.0_DP)*max(1.0_DP,abs(end_time))
    end function
end module pod_uq_hfem_ads_module
