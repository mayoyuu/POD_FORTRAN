!> Endpoint statistics for the SRP study ADS propagation.
!! Mean/covariance come from a fixed-seed uniform sample of the normalized
!! parameter box evaluated on the split manifold. d(x)/d(beta) at beta=0 is a
!! high-precision one-sided finite difference (the beta=0 reference is the edge
!! of the [0,B] box, not the expansion centre, so it is NOT read off the DA).
module pod_crtbp_srp_study_stats_module
    use pod_global, only: DP
    use pod_dace_classes, only: dace_max_variables
    use pod_ads_split_module, only: manifold_type, mf_evaluate_points
    use pod_crtbp_srp_study_module, only: CrtbpSrpStudyConfig
    use pod_crtbp_srp_study_integrator_module, only: study_adaptive_integrate_real
    implicit none
    private

    public :: study_endpoint_statistics
    public :: study_beta_zero_sensitivity

contains

    !> Mean and covariance of the endpoint state over a fixed-seed uniform
    !! sample of xi in [-1,1]^nvar, evaluated on the split manifold. The affine
    !! xi<->physical map makes this identical to uniform physical sampling.
    subroutine study_endpoint_statistics(manifold, n_samples, seed, mean_state, cov_state, status)
        type(manifold_type), intent(in) :: manifold
        integer, intent(in) :: n_samples, seed
        real(DP), intent(out) :: mean_state(6), cov_state(6,6)
        integer, intent(out) :: status
        real(DP), allocatable :: xi(:,:), values(:,:)
        logical, allocatable :: found(:)
        real(DP) :: centered(6)
        integer :: nvar, n_out, i, j, k, n_ok

        status = 0
        mean_state = 0.0_DP
        cov_state = 0.0_DP
        if (n_samples < 2) then
            status = -1
            return
        end if
        nvar = dace_max_variables()
        if (manifold%n_patches < 1) then
            status = -2
            return
        end if
        n_out = manifold%patches(1)%da_vec%size
        if (n_out < 6) then
            status = -3
            return
        end if

        allocate(xi(nvar, n_samples), values(n_out, n_samples), found(n_samples))
        call fill_uniform_samples(nvar, n_samples, seed, xi)

        call mf_evaluate_points(manifold, xi, values, found, status)
        if (status /= 0) then
            deallocate(xi, values, found)
            return
        end if

        n_ok = 0
        do k = 1, n_samples
            if (.not. found(k)) then
                status = -4
                deallocate(xi, values, found)
                return
            end if
            n_ok = n_ok + 1
            do i = 1, 6
                mean_state(i) = mean_state(i) + values(i, k)
            end do
        end do
        mean_state = mean_state / real(n_ok, DP)

        do k = 1, n_samples
            do i = 1, 6
                centered(i) = values(i, k) - mean_state(i)
            end do
            do i = 1, 6
                do j = 1, 6
                    cov_state(i, j) = cov_state(i, j) + centered(i)*centered(j)
                end do
            end do
        end do
        cov_state = cov_state / real(n_ok - 1, DP)

        deallocate(xi, values, found)
    end subroutine study_endpoint_statistics

    !> d(x)/d(beta) at beta=0 (one-sided): (x(eps) - x(0)) / eps with u=v=0 and
    !! chi=chi0. At u=v=0 the flat-plate and cannonball accelerations coincide,
    !! so the result is model-independent for a given beta.
    subroutine study_beta_zero_sensitivity(x0, cfg, chi0, beta_eps, t_end, &
            rel_tol, abs_tol, dt_min, dt_max, max_steps, dsens)
        real(DP), intent(in) :: x0(6), chi0, beta_eps, t_end, rel_tol, abs_tol, dt_min, dt_max
        integer, intent(in) :: max_steps
        type(CrtbpSrpStudyConfig), intent(in) :: cfg
        real(DP), intent(out) :: dsens(6)
        real(DP) :: p0(4), pe(4), x_ref(6), x_eps(6)

        if (beta_eps <= 0.0_DP) error stop 'study_beta_zero_sensitivity: beta_eps must be positive'
        p0 = [0.0_DP, chi0, 0.0_DP, 0.0_DP]
        pe = [beta_eps, chi0, 0.0_DP, 0.0_DP]
        call study_adaptive_integrate_real(x0, 0.0_DP, t_end, rel_tol, abs_tol, &
            dt_min, dt_max, max_steps, cfg, p0, x_ref)
        call study_adaptive_integrate_real(x0, 0.0_DP, t_end, rel_tol, abs_tol, &
            dt_min, dt_max, max_steps, cfg, pe, x_eps)
        dsens = (x_eps - x_ref) / beta_eps
    end subroutine study_beta_zero_sensitivity

    !> Fill xi(nvar, n_samples) with deterministic uniform draws in [-1,1]
    !! from a fixed integer seed, using the intrinsic generator.
    subroutine fill_uniform_samples(nvar, n_samples, seed, xi)
        integer, intent(in) :: nvar, n_samples, seed
        real(DP), intent(out) :: xi(:,:)
        integer :: n_seed, i
        integer, allocatable :: s(:)
        call random_seed(size=n_seed)
        allocate(s(n_seed))
        do i = 1, n_seed
            s(i) = seed + i
        end do
        call random_seed(put=s)
        call random_number(xi)
        xi = 2.0_DP*xi - 1.0_DP
        deallocate(s)
    end subroutine fill_uniform_samples

end module pod_crtbp_srp_study_stats_module
