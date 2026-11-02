!> Validate the time-advancing ADS engine for the SRP study force model.
!! The root patch carries [state(6); params(4)] as one 10-component vector, so
!! patch_split applies the SAME affine xi substitution to state and params.
!! Check: the terminal manifold, evaluated at box corners + fixed interior
!! points, reproduces direct real propagation at the same physical parameters.
program test_crtbp_srp_study_ads
    use pod_global, only: DP
    use pod_dace_classes, only: dace_max_variables, active_da_count
    use pod_crtbp_module, only: set_crtbp_mu
    use pod_srp_study_parameters_module, only: SrpStudyDomain, srp_study_domain
    use pod_crtbp_srp_study_module, only: CrtbpSrpStudyConfig
    use pod_crtbp_srp_study_integrator_module, only: study_adaptive_integrate_real
    use pod_ads_split_module, only: manifold_type, mf_destroy, mf_evaluate_points
    use pod_crtbp_srp_study_ads_module, only: study_ads_log, study_ads_propagate
    implicit none
    real(DP), parameter :: orbit_period = 2.3017923284002024_DP
    real(DP), parameter :: rel_tol = 1.0e-10_DP, abs_tol = 1.0e-12_DP
    real(DP), parameter :: dt_min = 1.0e-13_DP
    integer,  parameter :: max_steps = 100000
    integer,  parameter :: da_order = 6
    real(DP) :: x0(6), tend, dt_max
    type(CrtbpSrpStudyConfig) :: cfg
    integer :: total_splits = 0, n_fail = 0

    call set_crtbp_mu(0.012153614091891635_DP)
    cfg%sun%alpha0_rad = 0.37_DP
    x0 = [1.1309265107780351_DP, 0.0_DP, 0.0_DP, 0.0_DP, -0.46540743845849059_DP, 0.0_DP]
    tend = 2.0_DP*orbit_period
    if (tend*cfg%sun%tu_s < 15.0_DP*86400.0_DP) tend = 3.0_DP*orbit_period
    dt_max = tend

    call run_ads_case('1D cannonball [0,0.2]', srp_study_domain(0.0_DP, 0.2_DP), 1.0e-10_DP, 8, 6)
    call run_ads_case('3D flatplate local', srp_study_domain(0.08_DP, 0.12_DP, tilt_max_deg=3.0_DP), 1.0e-9_DP, 6, 6)
    call run_ads_case('4D flatplate full', &
        srp_study_domain(0.08_DP, 0.12_DP, tilt_max_deg=3.0_DP, chi_lower=0.4_DP, chi_upper=0.8_DP), &
        1.0e-9_DP, 6, 6)
    ! Low order over a wide box forces the truncation-error estimate above toll,
    ! exercising the split path (same affine state+params substitution).
    call run_ads_case('1D cannonball [0,0.2] order-3', srp_study_domain(0.0_DP, 0.2_DP), 1.0e-8_DP, 10, 3)

    call check(active_da_count() == 0, 'all DA handles released after ADS')
    call check(total_splits > 0, 'ADS engine performs at least one split')
    if (n_fail /= 0) then
        print '(A,I0)', 'FAIL: CRTBP SRP study ADS checks: ', n_fail
        error stop 1
    end if
    print '(A,I0,A)', 'PASS: CRTBP SRP study ADS engine (', total_splits, ' total splits)'
contains

    subroutine check(ok, label)
        logical, intent(in) :: ok
        character(*), intent(in) :: label
        if (.not. ok) then
            n_fail = n_fail + 1
            print '(A,A)', 'FAIL: ', label
        end if
    end subroutine check

    subroutine run_ads_case(label, domain, err_toll, n_split_max, da_order_use)
        character(*), intent(in) :: label
        type(SrpStudyDomain), intent(in) :: domain
        real(DP), intent(in) :: err_toll
        integer, intent(in) :: n_split_max, da_order_use
        type(manifold_type) :: result_manifold
        type(study_ads_log) :: log
        integer :: index_map(4), nvar, n_segments, l, status, j
        real(DP) :: err_toll_vec(6), c(4), h(4), p(4), xout(6)
        real(DP) :: state_err, param_err, max_state_err, max_param_err
        real(DP), allocatable :: xi(:,:), values(:,:)
        logical, allocatable :: found(:)
        integer :: n_samples

        err_toll_vec = err_toll
        n_segments = 4

        call study_ads_propagate(x0, cfg, domain, index_map, da_order_use, err_toll_vec, &
            n_split_max, tend, n_segments, rel_tol, abs_tol, dt_min, dt_max, max_steps, &
            result_manifold, log, .false.)

        nvar = dace_max_variables()
        total_splits = total_splits + log%n_splits_total

        call generate_samples(nvar, xi, n_samples)
        allocate(values(10, n_samples), found(n_samples))
        call mf_evaluate_points(result_manifold, xi, values, found, status)
        call check(status == 0, label//': manifold evaluate status ok')

        c = domain%center(); h = domain%halfwidth()
        max_state_err = 0.0_DP
        max_param_err = 0.0_DP
        do l = 1, n_samples
            if (.not. found(l)) then
                n_fail = n_fail + 1
                print '(A,A,A,I0)', 'FAIL: ', label, ' sample outside manifold: ', l
                cycle
            end if
            p = c
            do j = 1, 4
                if (index_map(j) > 0) p(j) = c(j) + h(j)*xi(index_map(j), l)
            end do
            call study_adaptive_integrate_real(x0, 0.0_DP, tend, rel_tol, abs_tol, &
                dt_min, dt_max, max_steps, cfg, p, xout)
            state_err = maxval(abs(values(1:6, l) - xout))
            param_err = maxval(abs(values(7:10, l) - p))
            max_state_err = max(max_state_err, state_err)
            max_param_err = max(max_param_err, param_err)
        end do

        print '(A,1X,A,1X,A,ES10.3,A,ES10.3,A,I0,A,I0)', &
            '[study ADS test] ', label, 'state_err=', max_state_err, &
            ' param_err=', max_param_err, ' splits=', log%n_splits_total, &
            ' patches=', result_manifold%n_patches

        call check(max_state_err < 1.0e-6_DP, label//': manifold state matches direct integration')
        call check(max_param_err < 1.0e-10_DP, label//': manifold params reproduce physical box')

        call mf_destroy(result_manifold)
        deallocate(xi, values, found)
    end subroutine run_ads_case

    subroutine generate_samples(nvar, xi, n_samples)
        integer, intent(in) :: nvar
        real(DP), allocatable, intent(out) :: xi(:,:)
        integer, intent(out) :: n_samples
        integer :: n_corner, i, j, k
        real(DP) :: pattern(8,4)
        n_corner = 2**nvar
        n_samples = n_corner + 8
        allocate(xi(nvar, n_samples))
        xi = 0.0_DP
        do i = 1, n_corner
            do j = 1, nvar
                if (btest(i-1, j-1)) then
                    xi(j, i) = 1.0_DP
                else
                    xi(j, i) = -1.0_DP
                end if
            end do
        end do
        pattern = 0.0_DP
        pattern(1,:) = [ 0.37_DP, -0.21_DP,  0.55_DP, -0.63_DP]
        pattern(2,:) = [-0.19_DP,  0.44_DP, -0.27_DP,  0.71_DP]
        pattern(3,:) = [ 0.61_DP,  0.13_DP, -0.48_DP,  0.05_DP]
        pattern(4,:) = [-0.33_DP, -0.57_DP,  0.23_DP, -0.19_DP]
        pattern(5,:) = [ 0.09_DP,  0.81_DP, -0.15_DP,  0.39_DP]
        pattern(6,:) = [-0.77_DP,  0.29_DP,  0.41_DP, -0.51_DP]
        pattern(7,:) = [ 0.28_DP, -0.66_DP, -0.73_DP,  0.12_DP]
        pattern(8,:) = [-0.47_DP,  0.03_DP,  0.65_DP,  0.88_DP]
        do k = 1, 8
            do j = 1, nvar
                xi(j, n_corner + k) = pattern(k, j)
            end do
        end do
    end subroutine generate_samples

end program test_crtbp_srp_study_ads
