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
        patch_init, patch_destroy, patch_split, &
        mf_init, mf_destroy, mf_push, mf_pop_front
    use pod_srp_study_parameters_module, only: SrpStudyDomain, SRP_PARAMETER_NAMES, &
        srp_study_init_da
    use pod_crtbp_srp_study_module, only: CrtbpSrpStudyConfig
    use pod_crtbp_srp_study_integrator_module, only: StudyIntegratorTempPool, &
        study_adaptive_integrate_da
    implicit none
    private

    integer, parameter :: N_STATE = 6, N_PARAM = 4, N_TOTAL = 10

    !> Split diagnostics collected during one ADS propagation.
    type, public :: study_ads_log
        integer  :: n_splits_total = 0
        integer  :: split_count_dir(4) = 0   ! per physical axis: beta, chi, u, v
        integer  :: n_max_active = 0         ! largest number of live patches
        real(DP) :: first_split_time = -1.0_DP
        real(DP) :: t_end = 0.0_DP
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
            rel_tol, abs_tol, dt_min, dt_max, max_steps)
        type(AlgebraicVector), intent(inout) :: da_vec
        real(DP), intent(in) :: tau0, tau1, rel_tol, abs_tol, dt_min, dt_max
        integer, intent(in) :: max_steps
        type(CrtbpSrpStudyConfig), intent(in) :: cfg
        type(StudyIntegratorTempPool), intent(inout) :: pool
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
            dt_min, dt_max, max_steps, cfg, params4, pool, final6)

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

    !> Time-advancing ADS propagation of the SRP parameter box.
    !!
    !! x0(6): constant initial state. domain: physical parameter box. The DA
    !! variables are the normalized xi coordinates of the active parameters.
    !! At each of n_segments time checkpoints the state map's truncation error
    !! is estimated; patches over err_toll are split along the dominant xi axis.
    subroutine study_ads_propagate(x0, cfg, domain, index_map, da_order, err_toll, &
            n_split_max, t_end, n_segments, rel_tol, abs_tol, dt_min, dt_max, max_steps, &
            result_manifold, log, verbose)
        real(DP), intent(in) :: x0(6)
        type(CrtbpSrpStudyConfig), intent(in) :: cfg
        type(SrpStudyDomain), intent(in) :: domain
        integer, intent(out) :: index_map(4)
        integer, intent(in) :: da_order, n_split_max, n_segments, max_steps
        real(DP), intent(in) :: err_toll(6), t_end, rel_tol, abs_tol, dt_min, dt_max
        type(manifold_type), intent(out) :: result_manifold
        type(study_ads_log), intent(out) :: log
        logical, intent(in) :: verbose

        type(manifold_type) :: active, next_active
        type(patch_type) :: p, left, right
        type(AlgebraicVector) :: root_da_vec
        type(StudyIntegratorTempPool) :: pool
        real(DP), allocatable :: scores(:)
        real(DP) :: tau0, tau1, dtau, errors(6), top, toll
        integer :: nvar, i, j, dir, axis, depth, k

        nvar = domain%n_variables()
        if (nvar < 1) error stop 'study_ads: domain has no active parameter'
        if (n_segments < 1) error stop 'study_ads: n_segments must be positive'

        call dace_initialize(da_order, nvar)

        call build_root_da_vec(x0, domain, index_map, root_da_vec)

        call mf_init(active)
        call patch_init(p, root_da_vec)
        call mf_push(active, p)
        call patch_destroy(p)

        call pool%init()
        allocate(scores(nvar))

        log%n_splits_total = 0
        log%split_count_dir = 0
        log%n_max_active = active%n_patches
        log%first_split_time = -1.0_DP
        log%t_end = t_end
        toll = maxval(err_toll)

        dtau = t_end / real(n_segments, DP)

        do k = 1, n_segments
            tau0 = real(k-1, DP) * dtau
            tau1 = real(k, DP) * dtau

            call mf_init(next_active)

            do while (active%n_patches > 0)
                call mf_pop_front(active, p)

                call advance_patch_state(p%da_vec, tau0, tau1, cfg, pool, &
                    rel_tol, abs_tol, dt_min, dt_max, max_steps)

                errors = 0.0_DP
                do i = 1, N_STATE
                    call da_estim_norm(p%da_vec%elements(i)%handle, 0, da_order, errors(i))
                end do

                depth = sh_count(p%history, 0)

                if (maxval(errors) <= toll .or. depth >= n_split_max) then
                    call mf_push(next_active, p)
                    if (verbose) then
                        write(*,'(A,I0,A,I0,A,ES10.3)') '[study ADS] t=', k, &
                            ' depth=', depth, ' err=', maxval(errors)
                    end if
                else
                    scores = 0.0_DP
                    do i = 1, N_STATE
                        do j = 1, nvar
                            call da_estim_norm(p%da_vec%elements(i)%handle, j, da_order, top)
                            scores(j) = scores(j) + top
                        end do
                    end do
                    dir = maxloc(scores, dim=1)
                    call patch_split(p, dir, left, right)
                    axis = physical_axis_of_dir(dir, index_map)
                    log%n_splits_total = log%n_splits_total + 1
                    if (axis >= 1 .and. axis <= 4) log%split_count_dir(axis) = log%split_count_dir(axis) + 1
                    if (log%first_split_time < 0.0_DP) log%first_split_time = tau1
                    call mf_push(next_active, left)
                    call mf_push(next_active, right)
                    if (verbose) then
                        write(*,'(A,I0,A,I0,A,A,A,I0,A,I0)') '[study ADS] t=', k, &
                            ' split dir=', dir, ' (', trim(SRP_PARAMETER_NAMES(axis)), &
                            ') depth=', depth, ' -> n_active=', next_active%n_patches
                    end if
                    call patch_destroy(left)
                    call patch_destroy(right)
                end if

                call patch_destroy(p)
            end do

            log%n_max_active = max(log%n_max_active, next_active%n_patches)
            call move_manifold(next_active, active)
        end do

        call move_manifold(active, result_manifold)

        deallocate(scores)
        call pool%destroy()

        if (verbose) then
            write(*,'(A,I0,A,I0,A,F10.5)') '[study ADS] done: ', result_manifold%n_patches, &
                ' patches, ', log%n_splits_total, ' splits, first split t=', log%first_split_time
        end if
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
