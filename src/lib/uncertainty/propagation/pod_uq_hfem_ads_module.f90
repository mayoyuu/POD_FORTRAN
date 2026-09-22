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
        ads_physical_to_unit
    use pod_da_integrator_module, only: da_adaptive_step_integrate, METHOD_RKF78
    use pod_da_force_model_module, only: set_propagation_epoch, &
        set_srp_scale_uncertainty, clear_srp_scale_uncertainty
    implicit none
    private

    integer, parameter :: STATE_DIM = 6

    type, public :: hfem_ads_options_type
        integer :: coordinate_mode = ADS_COORD_COMPONENT
        real(DP) :: domain_sigma = 3.0_DP
        real(DP) :: srp_sigma = 0.0_DP
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
        real(DP) :: domain_sigma = 3.0_DP
        real(DP) :: srp_sigma = 0.0_DP
        real(DP) :: basis6(6,6) = 0.0_DP
        real(DP) :: elapsed_seconds = 0.0_DP
        integer, allocatable :: split_counts(:)
    end type hfem_ads_stats_type

    public :: hfem_ads_propagate

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
        real(DP) :: eta, unit_tolerance

        status = 0
        message = ''
        call clear_srp_scale_uncertainty()
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
            coordinate_map, map_status, message)
        if (map_status /= 0) then
            status = map_status
            return
        end if
        stats%coordinate_mode = coordinate_map%mode
        stats%n_variables = coordinate_map%n_variables
        stats%domain_sigma = coordinate_map%domain_sigma
        stats%srp_sigma = coordinate_map%srp_sigma
        stats%basis6 = coordinate_map%basis6
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
            eta = 0.0_DP
            if (coordinate_map%n_variables == 7) eta = input_samples(7,j)
            call ads_physical_to_unit(coordinate_map, &
                input_samples(1:6,j)-nominal_state, eta, unit_point, map_status)
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
            options, domain, stats)
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
            if (coordinate_map%n_variables == 7) &
                output_samples(7,i) = input_samples(7,inside_indices(j))
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
        end subroutine fail
    end subroutine hfem_ads_propagate

    subroutine build_domain(nominal_state, epoch0, t_start, t_end, coordinate_map, &
            options, result_domain, stats)
        real(DP), intent(in) :: nominal_state(6), epoch0, t_start, t_end
        type(ads_coordinate_map_type), intent(in) :: coordinate_map
        type(hfem_ads_options_type), intent(in) :: options
        type(manifold_type), intent(out) :: result_domain
        type(hfem_ads_stats_type), intent(inout) :: stats
        type(manifold_type) :: queue
        type(patch_type) :: input_patch, propagated_patch, left, right
        type(AlgebraicVector) :: initial_da, final_da_nondim, final_da_physical
        real(DP), allocatable :: times(:), nominal_states(:,:)
        real(DP), allocatable :: patch_center(:), patch_width(:)
        real(DP) :: errors(6), excess(6), state_scale, eta_scale
        integer :: i, j, n_steps, component(1), direction

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
            if (coordinate_map%n_variables == 7) then
                patch_center = sh_center(input_patch%history)
                patch_width = sh_width(input_patch%history)
                eta_scale = coordinate_map%domain_sigma*coordinate_map%srp_sigma
                call set_srp_scale_uncertainty(7, eta_scale*patch_center(7), &
                    eta_scale*patch_width(7)/2.0_DP)
            else
                call clear_srp_scale_uncertainty()
            end if

            if (t_end == t_start) then
                final_da_nondim = input_patch%da_vec
            else
                call da_adaptive_step_integrate(input_patch%da_vec, &
                    t_start/config%TU, t_end/config%TU, METHOD_RKF78, times, &
                    nominal_states, final_da_nondim, n_steps, options%max_steps, &
                    options%rel_tol, options%abs_tol, options%dt_min, options%dt_max)
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
                    component(1), options%da_order)
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
        call clear_srp_scale_uncertainty()
        call mf_destroy(queue)
    end subroutine build_domain

end module pod_uq_hfem_ads_module
