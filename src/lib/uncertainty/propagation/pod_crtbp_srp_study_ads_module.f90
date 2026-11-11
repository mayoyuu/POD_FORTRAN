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
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private

    integer, parameter :: N_STATE = 6, N_PARAM = 4, N_TOTAL = 10

    !> Adaptive ADS macro-step controls. All values use the same time
    !! coordinate as t_end (nondimensional time in the CRTBP study).
    type, public :: study_ads_time_options
        real(DP) :: dt_initial = 0.0_DP
        real(DP) :: dt_min = 0.0_DP
        real(DP) :: dt_max = 0.0_DP
        real(DP) :: safety = 0.8_DP
        real(DP) :: growth_max = 2.0_DP
        real(DP) :: shrink_min = 0.2_DP
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
        integer  :: n_ads_steps_accepted = 0
        integer  :: n_ads_steps_rejected = 0
        integer  :: n_depth_limited = 0
        integer  :: n_numerical_failures = 0
        integer  :: n_diverged = 0
        real(DP) :: first_split_time = -1.0_DP
        real(DP) :: first_split_time_lower = -1.0_DP
        real(DP) :: first_split_time_upper = -1.0_DP
        real(DP) :: first_numerical_failure_time = -1.0_DP
        real(DP) :: first_divergence_time = -1.0_DP
        real(DP) :: min_accepted_dt = 0.0_DP
        real(DP) :: max_accepted_dt = 0.0_DP
        real(DP) :: t_end = 0.0_DP
        type(study_ads_split_event), allocatable :: split_events(:)
    end type

    public :: study_ads_propagate
    public :: study_ads_patch_bounds

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
            err_toll, n_split_max, rel_tol, abs_tol, dt_min, dt_max, max_steps, &
            pool, result, log, verbose)
        type(manifold_type), intent(in) :: active
        real(DP), intent(in) :: tau0, tau1, err_toll(N_STATE)
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
                log%n_depth_limited = log%n_depth_limited + 1
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


    !> Time-advancing ADS propagation of the SRP parameter box.
    !!
    !! x0(6): constant initial state. domain: physical parameter box. The DA
    !! variables are the normalized xi coordinates of the active parameters.
    !! A synchronized adaptive macro-step is accepted only while every live
    !! patch satisfies err_toll. Failed trials roll back and shrink in time;
    !! at dt_min the offending input patch is split and both children repropagate.
    subroutine study_ads_propagate(x0, cfg, domain, index_map, da_order, err_toll, &
            n_split_max, t_end, time_options, rel_tol, abs_tol, dt_min, dt_max, max_steps, &
            result_manifold, log, verbose)
        real(DP), intent(in) :: x0(N_STATE)
        type(CrtbpSrpStudyConfig), intent(in) :: cfg
        type(SrpStudyDomain), intent(in) :: domain
        integer, intent(out) :: index_map(N_PARAM)
        integer, intent(in) :: da_order, n_split_max, max_steps
        real(DP), intent(in) :: err_toll(N_STATE), t_end, rel_tol, abs_tol, dt_min, dt_max
        type(study_ads_time_options), intent(in) :: time_options
        type(manifold_type), intent(out) :: result_manifold
        type(study_ads_log), intent(out) :: log
        logical, intent(in) :: verbose

        type(manifold_type) :: active, trial, event_result
        type(patch_type) :: p, q
        type(AlgebraicVector) :: root_da_vec
        type(StudyIntegratorTempPool) :: pool
        type(StudyIntegrationStatus) :: integration_status
        real(DP) :: tau, tau1, dt_ads, dt_trial, errors(N_STATE)
        real(DP) :: ratio, max_ratio, fail_ratio, factor, exponent
        integer :: nvar, i, depth, split_capacity
        logical :: finite, failed

        if (cfg%attitude_coordinates /= domain%attitude_coordinates) &
            error stop 'study_ads: force and domain attitude coordinates must match'
        nvar = domain%n_variables()
        if (nvar < 1) error stop 'study_ads: domain has no active parameter'
        if (da_order < 1 .or. n_split_max < 0) error stop 'study_ads: invalid DA order or split depth'
        if (any(err_toll <= 0.0_DP) .or. .not. all(ieee_is_finite(err_toll))) &
            error stop 'study_ads: tolerances must be positive and finite'
        if (.not. ieee_is_finite(time_options%dt_min) .or. &
            .not. ieee_is_finite(time_options%dt_initial) .or. &
            .not. ieee_is_finite(time_options%dt_max) .or. &
            time_options%dt_min <= 0.0_DP .or. &
            time_options%dt_initial < time_options%dt_min .or. &
            time_options%dt_max < time_options%dt_initial) &
            error stop 'study_ads: require 0 < dt_min <= dt_initial <= dt_max'
        if (time_options%safety <= 0.0_DP .or. time_options%safety > 1.0_DP .or. &
            time_options%growth_max < 1.0_DP .or. time_options%shrink_min <= 0.0_DP .or. &
            time_options%shrink_min >= 1.0_DP) error stop 'study_ads: invalid step controller factors'

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

        tau = 0.0_DP
        dt_ads = min(time_options%dt_initial, time_options%dt_max, &
            max(t_end, time_options%dt_min))
        exponent = 1.0_DP/real(da_order+1, DP)

        do while (tau < t_end - 100.0_DP*epsilon(1.0_DP)*max(1.0_DP,abs(t_end)))
            if (active%n_patches == 0) exit
            dt_trial = min(dt_ads, t_end-tau)
            tau1 = tau + dt_trial
            call mf_init(trial)
            failed = .false.
            max_ratio = 0.0_DP
            fail_ratio = huge(1.0_DP)

            do i = 1, active%n_patches
                call patch_copy(active%patches(i), q)
                call advance_patch_state(q%da_vec, tau, tau1, cfg, pool, &
                    rel_tol, abs_tol, dt_min, dt_max, max_steps, integration_status)
                call patch_error_ratio(q, da_order, err_toll, errors, ratio, finite)
                depth = sh_count(active%patches(i)%history, 0)

                if (integration_status%completed .and. finite .and. ratio <= 1.0_DP) then
                    call mf_push(trial, q)
                    if (depth < n_split_max) max_ratio = max(max_ratio, ratio)
                else
                    failed = .true.
                    fail_ratio = ratio
                    call patch_destroy(q)
                    exit
                end if
                call patch_destroy(q)
            end do

            if (.not. failed) then
                call mf_destroy(active)
                call move_manifold(trial, active)
                tau = tau1
                log%n_ads_steps_accepted = log%n_ads_steps_accepted + 1
                if (log%min_accepted_dt == 0.0_DP) log%min_accepted_dt = dt_trial
                log%min_accepted_dt = min(log%min_accepted_dt, dt_trial)
                log%max_accepted_dt = max(log%max_accepted_dt, dt_trial)
                factor = time_options%safety*max(max_ratio, 1.0e-12_DP)**(-exponent)
                factor = max(time_options%shrink_min, min(time_options%growth_max, factor))
                dt_ads = min(time_options%dt_max, max(time_options%dt_min, dt_trial*factor))
            else
                call mf_destroy(trial)
                log%n_ads_steps_rejected = log%n_ads_steps_rejected + 1
                if (dt_trial > time_options%dt_min*(1.0_DP+10.0_DP*epsilon(1.0_DP))) then
                    if (ieee_is_finite(fail_ratio) .and. fail_ratio > 1.0_DP) then
                        factor = time_options%safety*fail_ratio**(-exponent)
                        factor = min(0.8_DP, max(time_options%shrink_min, factor))
                    else
                        factor = max(time_options%shrink_min, 0.5_DP)
                    end if
                    dt_ads = max(time_options%dt_min, dt_trial*factor)
                else
                    call process_event_interval(active, tau, tau1, cfg, index_map, da_order, &
                        err_toll, n_split_max, rel_tol, abs_tol, dt_min, dt_max, max_steps, &
                        pool, event_result, log, verbose)
                    call mf_destroy(active)
                    call move_manifold(event_result, active)
                    tau = tau1
                    log%n_ads_steps_accepted = log%n_ads_steps_accepted + 1
                    if (log%min_accepted_dt == 0.0_DP) log%min_accepted_dt = dt_trial
                    log%min_accepted_dt = min(log%min_accepted_dt, dt_trial)
                    log%max_accepted_dt = max(log%max_accepted_dt, dt_trial)
                    dt_ads = min(time_options%dt_max, max(time_options%dt_min, &
                        dt_trial*time_options%growth_max))
                end if
            end if
        end do

        call move_manifold(active, result_manifold)
        log%n_depth_limited = 0
        do i = 1, result_manifold%n_patches
            depth = sh_count(result_manifold%patches(i)%history, 0)
            call patch_error_ratio(result_manifold%patches(i), da_order, err_toll, &
                errors, ratio, finite)
            if (depth >= n_split_max .and. finite .and. ratio > 1.0_DP) &
                log%n_depth_limited = log%n_depth_limited + 1
        end do
        call pool%destroy()
        log%n_diverged = log%n_numerical_failures
        log%first_divergence_time = log%first_numerical_failure_time

        if (verbose) write(*,'(A,I0,A,I0,A,ES12.4)') '[study ADS] done: ', &
            result_manifold%n_patches, ' patches, splits=', log%n_splits_total, &
            ' first split=', log%first_split_time
    end subroutine study_ads_propagate

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
