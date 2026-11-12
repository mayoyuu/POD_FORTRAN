!> Time-advancing Automatic Domain Splitting (ADS) for the CRTBP + SRP study
!! force model. Splits the normalized parameter box (xi in [-1,1]) so that the
!! sixth-order state map stays below a truncation-error threshold as the DA
!! flow grows nonlinear in time. A patch carries [state(6); params(4)] as one
!! 10-component AlgebraicVector so patch_split applies the SAME affine xi
!! substitution to both, keeping the split time-continuous (no root re-init).
module pod_crtbp_srp_study_ads_module
    use pod_global, only: DP
    use pod_dace_classes, only: DA, AlgebraicVector, dace_initialize, &
        da_estim_norm, active_da_count, assignment(=)
    use pod_ads_split_module, only: splitting_history_type, patch_type, manifold_type, &
        sh_count, sh_center, sh_width, sh_map_point, &
        patch_init, patch_copy, patch_destroy, patch_split, &
        mf_init, mf_destroy, mf_push, mf_pop_front
    use pod_srp_study_parameters_module, only: SrpStudyDomain, &
        srp_study_init_da
    use pod_crtbp_srp_study_module, only: CrtbpSrpStudyConfig
    use pod_crtbp_srp_study_integrator_module, only: StudyIntegratorTempPool, &
        StudyIntegrationStatus, study_adaptive_integrate_da
    use pod_crtbp_srp_study_stats_module, only: study_endpoint_statistics
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private

    integer, parameter :: N_STATE = 6, N_PARAM = 4, N_TOTAL = 10

    !> Per-checkpoint ADS statistics: state covariance and split status at one
    !! of the fixed time checkpoints.
    type, public :: study_checkpoint_stats
        real(DP) :: t_nd = 0.0_DP
        real(DP) :: mean(6) = 0.0_DP
        real(DP) :: cov(6,6) = 0.0_DP
        integer  :: n_patches = 0
        integer  :: n_diverged = 0
    end type

    !> One ADS split event, recorded for split-time and direction analysis.
    type, public :: study_ads_split_event
        real(DP) :: time_nd = 0.0_DP
        real(DP) :: time_lower = 0.0_DP
        real(DP) :: time_upper = 0.0_DP
        integer  :: axis = 0
        integer  :: depth = 0
        integer  :: n_active = 0
        real(DP) :: center = 0.0_DP
        real(DP) :: halfwidth = 0.0_DP
        real(DP) :: criterion = 0.0_DP
        real(DP) :: center_all(4) = 0.0_DP
        real(DP) :: width_all(4) = 0.0_DP
    end type

    type, public :: study_ads_log
        integer  :: n_splits_total = 0
        integer  :: split_count_dir(4) = 0
        integer  :: n_max_active = 0
        integer  :: n_depth_limited = 0
        integer  :: n_numerical_failures = 0
        integer  :: n_diverged = 0
        logical  :: case_diverged = .false.
        real(DP) :: first_split_time = -1.0_DP
        real(DP) :: first_split_time_lower = -1.0_DP
        real(DP) :: first_split_time_upper = -1.0_DP
        real(DP) :: first_split_time_accurate = -1.0_DP
        real(DP) :: first_numerical_failure_time = -1.0_DP
        real(DP) :: first_divergence_time = -1.0_DP
        real(DP) :: t_end = 0.0_DP
        type(study_ads_split_event), allocatable :: split_events(:)
        type(study_checkpoint_stats), allocatable :: checkpoints(:)
    end type

    public :: study_ads_propagate
    public :: study_ads_patch_bounds
    public :: detect_first_split_time

contains

    !> Build the root patch vector [state(6); params(4)]. State is the constant
    !! initial condition; params are the affine xi map (srp_study_init_da).
    subroutine build_root_da_vec(x0, domain, index_map, da_vec)
        real(DP), intent(in) :: x0(6)
        type(SrpStudyDomain), intent(in) :: domain
        integer, intent(out) :: index_map(4)
        type(AlgebraicVector), intent(out) :: da_vec
        type(DA) :: params(4)
        integer :: i

        call srp_study_init_da(domain, params, index_map)
        call da_vec%init(N_TOTAL)
        do i = 1, N_STATE
            da_vec%elements(i) = x0(i)
        end do
        do i = 1, N_PARAM
            da_vec%elements(N_STATE+i) = params(i)
        end do
        do i = 1, N_PARAM
            call params(i)%destroy()
        end do
    end subroutine build_root_da_vec

    !> Advance the state block (components 1:6) of a patch vector from tau0 to
    !! tau1 using the current params block (components 7:10). The params block
    !! is unchanged over the arc; only the state block is written back.
    subroutine advance_patch_state(da_vec, tau0, tau1, cfg, pool, &
            rel_tol, abs_tol, dt_min, dt_max, max_steps, integration_status)
        type(AlgebraicVector), intent(inout) :: da_vec
        real(DP), intent(in) :: tau0, tau1, rel_tol, abs_tol, dt_min, dt_max
        integer, intent(in) :: max_steps
        type(CrtbpSrpStudyConfig), intent(in) :: cfg
        type(StudyIntegratorTempPool), intent(inout) :: pool
        type(StudyIntegrationStatus), intent(out) :: integration_status
        type(AlgebraicVector) :: state6, final6
        type(DA) :: params4(4)
        integer :: i

        call state6%init(N_STATE)
        do i = 1, N_STATE
            state6%elements(i) = da_vec%elements(i)
        end do
        do i = 1, N_PARAM
            call params4(i)%init()
            params4(i) = da_vec%elements(N_STATE+i)
        end do
        call final6%init(N_STATE)
        call study_adaptive_integrate_da(state6, tau0, tau1, rel_tol, abs_tol, &
            dt_min, dt_max, max_steps, cfg, params4, pool, final6, status=integration_status)
        do i = 1, N_STATE
            da_vec%elements(i) = final6%elements(i)
        end do
        call state6%destroy()
        do i = 1, N_PARAM
            call params4(i)%destroy()
        end do
        call final6%destroy()
    end subroutine advance_patch_state

    !> Move patches from one manifold to another without deep copies.
    subroutine move_manifold(src, dst)
        type(manifold_type), intent(inout) :: src, dst
        call mf_destroy(dst)
        call move_alloc(src%patches, dst%patches)
        dst%n_patches = src%n_patches
        src%n_patches = 0
    end subroutine move_manifold

    !> Map a DA split direction (1..nvar) onto the physical parameter index and
    !! bump that axis's split counter. Returns 0 for an unknown direction.
    integer function physical_axis_of_dir(dir, index_map) result(axis)
        integer, intent(in) :: dir, index_map(4)
        integer :: j
        axis = 0
        do j = 1, 4
            if (index_map(j) == dir) then
                axis = j
                return
            end if
        end do
    end function physical_axis_of_dir

    subroutine patch_error_ratio(p, da_order, err_toll, errors, ratio, finite)
        type(patch_type), intent(in) :: p
        integer, intent(in) :: da_order
        real(DP), intent(in) :: err_toll(N_STATE)
        real(DP), intent(out) :: errors(N_STATE), ratio
        logical, intent(out) :: finite
        integer :: i

        errors = 0.0_DP
        finite = .true.
        do i = 1, N_STATE
            call da_estim_norm(p%da_vec%elements(i)%handle, 0, da_order, errors(i))
            if (.not. ieee_is_finite(errors(i))) finite = .false.
        end do
        if (finite) then
            ratio = maxval(errors/err_toll)
        else
            ratio = huge(1.0_DP)
        end if
    end subroutine patch_error_ratio

    integer function choose_split_direction(p, nvar, da_order) result(dir)
        type(patch_type), intent(in) :: p
        integer, intent(in) :: nvar, da_order
        real(DP), allocatable :: scores(:), widths(:)
        real(DP) :: top
        integer :: i, j

        allocate(scores(nvar))
        scores = 0.0_DP
        do i = 1, N_STATE
            do j = 1, nvar
                call da_estim_norm(p%da_vec%elements(i)%handle, j, da_order, top)
                if (ieee_is_finite(top)) scores(j) = scores(j) + max(0.0_DP, top)
            end do
        end do
        if (maxval(scores) > 0.0_DP) then
            dir = maxloc(scores, dim=1)
        else
            widths = sh_width(p%history)
            dir = maxloc(widths(1:nvar), dim=1)
            deallocate(widths)
        end if
        deallocate(scores)
    end function choose_split_direction

    subroutine record_split(log, p, dir, axis, ratio, tau0, tau1, index_map, n_active)
        type(study_ads_log), intent(inout) :: log
        type(patch_type), intent(in) :: p
        integer, intent(in) :: dir, axis, index_map(N_PARAM), n_active
        real(DP), intent(in) :: ratio, tau0, tau1
        real(DP), allocatable :: gc(:), gw(:)
        integer :: j, event_index

        log%n_splits_total = log%n_splits_total + 1
        event_index = log%n_splits_total
        if (event_index > size(log%split_events)) &
            error stop 'study ADS: split-event capacity exceeded'
        if (axis >= 1 .and. axis <= N_PARAM) &
            log%split_count_dir(axis) = log%split_count_dir(axis) + 1
        if (log%first_split_time < 0.0_DP) then
            log%first_split_time = tau1
            log%first_split_time_lower = tau0
            log%first_split_time_upper = tau1
        end if

        gc = sh_center(p%history)
        gw = sh_width(p%history)
        log%split_events(event_index)%time_nd = tau1
        log%split_events(event_index)%time_lower = tau0
        log%split_events(event_index)%time_upper = tau1
        log%split_events(event_index)%axis = axis
        log%split_events(event_index)%depth = sh_count(p%history, 0)
        log%split_events(event_index)%n_active = n_active
        log%split_events(event_index)%center = gc(dir)
        log%split_events(event_index)%halfwidth = 0.5_DP*gw(dir)
        log%split_events(event_index)%criterion = ratio
        do j = 1, N_PARAM
            if (index_map(j) > 0) then
                log%split_events(event_index)%center_all(j) = gc(index_map(j))
                log%split_events(event_index)%width_all(j) = gw(index_map(j))
            end if
        end do
        deallocate(gc, gw)
    end subroutine record_split

    subroutine process_event_interval(active, tau0, tau1, cfg, index_map, da_order, &
            err_toll, tol_div, n_split_max, rel_tol, abs_tol, dt_min, dt_max, max_steps, &
            pool, result, log, verbose)
        type(manifold_type), intent(in) :: active
        real(DP), intent(in) :: tau0, tau1, err_toll(N_STATE), tol_div(N_STATE)
        integer, intent(in) :: index_map(N_PARAM), da_order, n_split_max, max_steps
        real(DP), intent(in) :: rel_tol, abs_tol, dt_min, dt_max
        type(CrtbpSrpStudyConfig), intent(in) :: cfg
        type(StudyIntegratorTempPool), intent(inout) :: pool
        type(manifold_type), intent(out) :: result
        type(study_ads_log), intent(inout) :: log
        logical, intent(in) :: verbose

        type(manifold_type) :: queue
        type(patch_type) :: p, q, left, right
        type(StudyIntegrationStatus) :: integration_status
        real(DP) :: errors(N_STATE), ratio
        integer :: i, nvar, depth, dir, axis
        logical :: finite, valid

        nvar = count(index_map > 0)
        call mf_init(queue)
        call mf_init(result)
        do i = 1, active%n_patches
            call mf_push(queue, active%patches(i))
        end do

        do while (queue%n_patches > 0)
            call mf_pop_front(queue, p)
            call patch_copy(p, q)
            call advance_patch_state(q%da_vec, tau0, tau1, cfg, pool, &
                rel_tol, abs_tol, dt_min, dt_max, max_steps, integration_status)
            call patch_error_ratio(q, da_order, err_toll, errors, ratio, finite)
            depth = sh_count(p%history, 0)
            valid = integration_status%completed .and. finite .and. ratio <= 1.0_DP

            if (valid) then
                call mf_push(result, q)
            else if (integration_status%completed .and. finite .and. depth >= n_split_max) then
                if (any(errors > tol_div)) then
                    q%diverged = .true.
                    if (log%first_divergence_time < 0.0_DP) log%first_divergence_time = tau1
                end if
                call mf_push(result, q)
            else if (depth < n_split_max) then
                if (integration_status%completed .and. finite) then
                    dir = choose_split_direction(q, nvar, da_order)
                else
                    dir = choose_split_direction(p, nvar, da_order)
                end if
                axis = physical_axis_of_dir(dir, index_map)
                call patch_split(p, dir, left, right)
                call mf_push(queue, left)
                call mf_push(queue, right)
                call record_split(log, p, dir, axis, ratio, tau0, tau1, index_map, &
                    result%n_patches + queue%n_patches)
                log%n_max_active = max(log%n_max_active, result%n_patches + queue%n_patches)
                if (verbose) write(*,'(A,ES12.4,A,I0,A,I0)') &
                    '[study ADS] split t=', tau1, ' dir=', dir, ' depth=', depth
                call patch_destroy(left)
                call patch_destroy(right)
            else
                log%n_numerical_failures = log%n_numerical_failures + 1
                if (log%first_numerical_failure_time < 0.0_DP) &
                    log%first_numerical_failure_time = tau1
            end if

            call patch_destroy(q)
            call patch_destroy(p)
        end do
        call mf_destroy(queue)
    end subroutine process_event_interval


    !> Time-advancing ADS propagation of the SRP parameter box over a fixed
    !! sequence of n_segments equal checkpoints.
    !!
    !! x0(6): constant initial state. domain: physical parameter box. The DA
    !! variables are the normalized xi coordinates of the active parameters.
    !! Each checkpoint advances the whole manifold from tau0 to tau1 with
    !! process_event_interval, which re-propagates every patch and splits any
    !! patch whose truncation error still exceeds err_toll (unless it is already
    !! at n_split_max depth, in which case it is kept and, when it exceeds
    !! tol_div, flagged diverged). Per-checkpoint mean/covariance and split
    !! status are recorded for covariance-time analysis.
    subroutine study_ads_propagate(x0, cfg, domain, index_map, da_order, err_toll, tol_div, &
            n_split_max, t_end, n_segments, n_samples, seed, &
            rel_tol, abs_tol, dt_min, dt_max, max_steps, &
            result_manifold, log, verbose)
        real(DP), intent(in) :: x0(N_STATE)
        type(CrtbpSrpStudyConfig), intent(in) :: cfg
        type(SrpStudyDomain), intent(in) :: domain
        integer, intent(out) :: index_map(N_PARAM)
        integer, intent(in) :: da_order, n_split_max, max_steps, n_segments, n_samples, seed
        real(DP), intent(in) :: err_toll(N_STATE), tol_div(N_STATE), t_end
        real(DP), intent(in) :: rel_tol, abs_tol, dt_min, dt_max
        type(manifold_type), intent(out) :: result_manifold
        type(study_ads_log), intent(out) :: log
        logical, intent(in) :: verbose

        type(manifold_type) :: active, next_active
        type(patch_type) :: p
        type(AlgebraicVector) :: root_da_vec
        type(StudyIntegratorTempPool) :: pool
        real(DP) :: tau0, tau1, dtau, mean6(N_STATE), cov66(N_STATE,N_STATE)
        real(DP) :: errors(N_STATE), ratio
        integer :: nvar, i, k, depth, split_capacity, stat, n_div
        logical :: finite

        if (cfg%attitude_coordinates /= domain%attitude_coordinates) &
            error stop 'study_ads: force and domain attitude coordinates must match'
        nvar = domain%n_variables()
        if (nvar < 1) error stop 'study_ads: domain has no active parameter'
        if (da_order < 1 .or. n_split_max < 0) error stop 'study_ads: invalid DA order or split depth'
        if (any(err_toll <= 0.0_DP) .or. .not. all(ieee_is_finite(err_toll))) &
            error stop 'study_ads: err_toll must be positive and finite'
        if (any(tol_div <= 0.0_DP) .or. .not. all(ieee_is_finite(tol_div))) &
            error stop 'study_ads: tol_div must be positive and finite'
        if (any(tol_div <= err_toll)) error stop 'study_ads: require tol_div > err_toll componentwise'
        if (n_segments < 1) error stop 'study_ads: n_segments must be >= 1'
        if (n_samples < 2) error stop 'study_ads: n_samples must be >= 2'

        call dace_initialize(da_order, nvar)
        call build_root_da_vec(x0, domain, index_map, root_da_vec)
        call mf_init(active)
        call patch_init(p, root_da_vec)
        call mf_push(active, p)
        call patch_destroy(p)
        call pool%init()

        log = study_ads_log()
        log%n_max_active = active%n_patches
        log%t_end = t_end
        split_capacity = max(1, 2**n_split_max - 1)
        allocate(log%split_events(split_capacity))
        allocate(log%checkpoints(n_segments))

        dtau = t_end / real(n_segments, DP)

        do k = 1, n_segments
            if (active%n_patches == 0) exit
            tau0 = real(k-1, DP) * dtau
            tau1 = real(k, DP) * dtau
            call process_event_interval(active, tau0, tau1, cfg, index_map, da_order, &
                err_toll, tol_div, n_split_max, rel_tol, abs_tol, dt_min, dt_max, max_steps, &
                pool, next_active, log, verbose)
            call move_manifold(next_active, active)
            log%n_max_active = max(log%n_max_active, active%n_patches)
            call study_endpoint_statistics(active, n_samples, seed, mean6, cov66, stat)
            log%checkpoints(k)%t_nd = tau1
            log%checkpoints(k)%mean = mean6
            log%checkpoints(k)%cov = cov66
            log%checkpoints(k)%n_patches = active%n_patches
            n_div = 0
            do i = 1, active%n_patches
                if (active%patches(i)%diverged) n_div = n_div + 1
            end do
            log%checkpoints(k)%n_diverged = n_div
        end do

        call move_manifold(active, result_manifold)
        log%n_depth_limited = 0
        log%n_diverged = 0
        do i = 1, result_manifold%n_patches
            depth = sh_count(result_manifold%patches(i)%history, 0)
            call patch_error_ratio(result_manifold%patches(i), da_order, err_toll, &
                errors, ratio, finite)
            if (depth >= n_split_max .and. finite .and. ratio > 1.0_DP) &
                log%n_depth_limited = log%n_depth_limited + 1
            if (result_manifold%patches(i)%diverged) log%n_diverged = log%n_diverged + 1
        end do
        log%case_diverged = log%n_diverged > 0
        call pool%destroy()

        if (verbose) write(*,'(A,I0,A,I0,A,I0,A,ES12.4)') '[study ADS] done: ', &
            result_manifold%n_patches, ' patches, splits=', log%n_splits_total, &
            ' diverged=', log%n_diverged, ' first split=', log%first_split_time
    end subroutine study_ads_propagate

    !> Accurate first-split time: bisect on the ROOT patch (no splitting) to find
    !! the earliest time its truncation error exceeds err_toll, to relative
    !! precision ~1e-3. Decoupled from the main propagation loop; returns -1 when
    !! the root patch stays within err_toll over the whole [0,t_end] arc.
    !! DACE must already be initialized (as it is after study_ads_propagate);
    !! it is NOT re-initialized here so the caller's manifold handles stay valid.
    subroutine detect_first_split_time(x0, cfg, domain, da_order, err_toll, &
            t_end, rel_tol, abs_tol, dt_min, dt_max, max_steps, &
            first_split_time, status)
        real(DP), intent(in) :: x0(N_STATE)
        type(CrtbpSrpStudyConfig), intent(in) :: cfg
        type(SrpStudyDomain), intent(in) :: domain
        integer, intent(in) :: da_order, max_steps
        real(DP), intent(in) :: err_toll(N_STATE), t_end, rel_tol, abs_tol, dt_min, dt_max
        real(DP), intent(out) :: first_split_time
        integer, intent(out) :: status
        type(AlgebraicVector) :: root_da_vec
        type(patch_type) :: p, q
        type(StudyIntegratorTempPool) :: pool
        type(StudyIntegrationStatus) :: integration_status
        real(DP) :: errors(N_STATE), ratio, t_lo, t_hi, t_mid
        integer :: index_map(N_PARAM), nvar, iter
        logical :: finite

        status = 0
        first_split_time = -1.0_DP
        nvar = domain%n_variables()
        if (nvar < 1) then
            status = -1
            return
        end if
        if (any(err_toll <= 0.0_DP) .or. .not. all(ieee_is_finite(err_toll))) then
            status = -2
            return
        end if
        if (t_end <= 0.0_DP) then
            status = -3
            return
        end if

        call build_root_da_vec(x0, domain, index_map, root_da_vec)
        call patch_init(p, root_da_vec)
        call pool%init()

        call patch_copy(p, q)
        call advance_patch_state(q%da_vec, 0.0_DP, t_end, cfg, pool, &
            rel_tol, abs_tol, dt_min, dt_max, max_steps, integration_status)
        call patch_error_ratio(q, da_order, err_toll, errors, ratio, finite)
        call patch_destroy(q)
        if (.not. (integration_status%completed .and. finite)) then
            status = -4
            call patch_destroy(p)
            call pool%destroy()
            return
        end if
        if (ratio <= 1.0_DP) then
            call patch_destroy(p)
            call pool%destroy()
            return
        end if

        t_lo = 0.0_DP
        t_hi = t_end
        do iter = 1, 12
            t_mid = 0.5_DP*(t_lo + t_hi)
            call patch_copy(p, q)
            call advance_patch_state(q%da_vec, 0.0_DP, t_mid, cfg, pool, &
                rel_tol, abs_tol, dt_min, dt_max, max_steps, integration_status)
            call patch_error_ratio(q, da_order, err_toll, errors, ratio, finite)
            call patch_destroy(q)
            if (integration_status%completed .and. finite .and. ratio > 1.0_DP) then
                t_hi = t_mid
            else
                t_lo = t_mid
            end if
        end do
        first_split_time = 0.5_DP*(t_lo + t_hi)

        call patch_destroy(p)
        call pool%destroy()
    end subroutine detect_first_split_time

    !> Physical [lower,upper] bounds of a patch from its split history and the
    !! root parameter domain. Used for split-direction and domain-bound logging.
    subroutine study_ads_patch_bounds(history, domain, index_map, lower, upper)
        type(splitting_history_type), intent(in) :: history
        type(SrpStudyDomain), intent(in) :: domain
        integer, intent(in) :: index_map(4)
        real(DP), intent(out) :: lower(4), upper(4)
        real(DP) :: c(4), h(4), gc(4), gw(4)
        integer :: j, k

        c = domain%center(); h = domain%halfwidth()
        gc = sh_center(history)   ! global xi center, length nvar
        gw = sh_width(history)    ! global xi width,  length nvar
        do j = 1, 4
            k = index_map(j)
            if (k > 0) then
                lower(j) = c(j) + h(j) * (gc(k) - 0.5_DP * gw(k))
                upper(j) = c(j) + h(j) * (gc(k) + 0.5_DP * gw(k))
            else
                lower(j) = c(j)
                upper(j) = c(j)
            end if
        end do
    end subroutine study_ads_patch_bounds

end module pod_crtbp_srp_study_ads_module
