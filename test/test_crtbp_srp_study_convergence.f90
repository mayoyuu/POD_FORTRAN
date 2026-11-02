!> Validation for the SRP study ADS engine at the fixed design
!! (DA order 8, err_tol 1e-8, n_split_max 12). For each representative orbit
!! (one per family) we report the max |ADS - direct integration| endpoint
!! state error, the patch/split/divergence counts, and wall time.
program test_crtbp_srp_study_convergence
    use pod_global, only: DP
    use pod_dace_classes, only: dace_max_variables, active_da_count
    use pod_crtbp_module, only: set_crtbp_mu
    use pod_srp_study_parameters_module, only: SrpStudyDomain, srp_study_domain
    use pod_crtbp_srp_study_module, only: CrtbpSrpStudyConfig, &
        STUDY_CANNONBALL, STUDY_FLAT_PLATE
    use pod_crtbp_srp_study_integrator_module, only: study_adaptive_integrate_real
    use pod_ads_split_module, only: manifold_type, mf_destroy, mf_evaluate_points
    use pod_crtbp_srp_study_ads_module, only: study_ads_log, study_ads_propagate
    implicit none

    real(DP), parameter :: MU = 0.0121536140918916_DP
    real(DP), parameter :: TU_S = 375237.080249162_DP
    real(DP), parameter :: REL_TOL = 1.0e-10_DP, ABS_TOL = 1.0e-12_DP
    real(DP), parameter :: DT_MIN = 1.0e-13_DP
    integer,  parameter :: MAX_STEPS = 100000
    integer,  parameter :: N_SEGMENTS = 48
    integer,  parameter :: N_SPLIT_MAX = 12
    integer,  parameter :: N_ORBITS = 5

    type :: orbit_t
        character(len=14) :: label
        real(DP) :: x0(6)
        real(DP) :: period_nd
    end type

    type(orbit_t) :: orbits(N_ORBITS)
    integer :: i, n_fail

    call set_crtbp_mu(MU)

    orbits(1) = make_orbit('DRO',     1.0919046585407821_DP, 0.0_DP, 0.0_DP, 0.0_DP, &
        -0.46197686112063113_DP, 0.0_DP, 1.551142049282021_DP)
    orbits(2) = make_orbit('Halo_L1', 0.82336585706409793_DP, 0.0_DP, 0.013024366884938094_DP, 0.0_DP, &
        0.12911335596440643_DP, 0.0_DP, 2.744124888564603_DP)
    orbits(3) = make_orbit('Halo_L2', 1.180747161876686_DP, 0.0_DP, 0.013009837283522473_DP, 0.0_DP, &
        -0.15684501482174368_DP, 0.0_DP, 3.414175928347213_DP)
    orbits(4) = make_orbit('Lyap_L1', 0.8421130711339766_DP, 0.0_DP, 0.0_DP, 0.0_DP, &
        -0.042066572390419173_DP, 0.0_DP, 2.696693232624513_DP)
    orbits(5) = make_orbit('Lyap_L2', 1.1504585909779124_DP, 0.0_DP, 0.0_DP, 0.0_DP, &
        0.027899769469798077_DP, 0.0_DP, 3.374369373601701_DP)

    n_fail = 0
    write(*,'(A)') '# SRP study ADS convergence'
    write(*,'(A12,1X,A4,1X,A2,1X,A5,1X,A8,1X,A9,1X,A9,1X,A8,1X,A7,1X,A9,1X,A7,1X,A9,1X,A7,1X,A7,1X,A7)') &
        'orbit', 'dim', 'mdl', 'order', 'toll', 'rel_pos', 'rel_vel', 'patches', 'splits', &
        't_split_d', 'nmax', 't_div_d', 'ndiv', 'wall_s', 'nout'

    ! Q1 1D cannonball at fixed design (order 8, err_tol 1e-8, all orbits).
    do i = 1, N_ORBITS
        call run_case(orbits(i), '1D', 'CB', 8, 1.0e-8_DP, &
            srp_study_domain(0.0_DP, 0.1_DP))
    end do

    ! Q2 3D flatplate at fixed design (DRO, Halo L1).
    do i = 1, N_ORBITS
        if (i /= 1 .and. i /= 2) cycle
        call run_case(orbits(i), '3D', 'FP', 8, 1.0e-8_DP, &
            srp_study_domain(0.08_DP, 0.12_DP, tilt_max_deg=10.0_DP))
    end do

    ! Q4 4D flatplate at fixed design (Halo L1).
    call run_case(orbits(2), '4D', 'FP', 8, 1.0e-8_DP, &
        srp_study_domain(0.08_DP, 0.12_DP, tilt_max_deg=10.0_DP, &
        chi_lower=0.4_DP, chi_upper=0.8_DP))

    call check(active_da_count() == 0, 'all DA handles released')
    if (n_fail /= 0) then
        print '(A,I0)', 'FAIL: SRP study convergence checks: ', n_fail
        error stop 1
    end if
    print '(A)', 'PASS: SRP study ADS convergence sweep (see table above)'
contains

    function make_orbit(label, x, y, z, vx, vy, vz, period) result(o)
        character(*), intent(in) :: label
        real(DP), intent(in) :: x, y, z, vx, vy, vz, period
        type(orbit_t) :: o
        o%label = label
        o%x0 = [x, y, z, vx, vy, vz]
        o%period_nd = period
    end function make_orbit

    subroutine check(ok, label)
        logical, intent(in) :: ok
        character(*), intent(in) :: label
        if (.not. ok) then
            n_fail = n_fail + 1
            print '(A,A)', 'FAIL: ', label
        end if
    end subroutine check

    subroutine run_case(orbit, dim, mdl, da_order, toll, domain)
        type(orbit_t), intent(in) :: orbit
        character(*), intent(in) :: dim, mdl
        integer, intent(in) :: da_order
        real(DP), intent(in) :: toll
        type(SrpStudyDomain), intent(in) :: domain
        type(CrtbpSrpStudyConfig) :: cfg
        type(manifold_type) :: result_manifold
        type(study_ads_log) :: log
        integer :: index_map(4), nvar, n_samples, status, l, j, n_outside
        real(DP) :: err_toll_vec(6), tend, dt_max, c(4), h(4), p(4), xout(6)
        real(DP) :: max_rel_pos, max_rel_vel, dr_norm, dv_norm, r_norm, v_norm
        real(DP), allocatable :: xi(:,:), values(:,:)
        logical, allocatable :: found(:)
        real(DP) :: wall_t0, wall_t1

        err_toll_vec = toll
        tend = 2.0_DP*orbit%period_nd
        if (tend*TU_S < 15.0_DP*86400.0_DP) tend = 3.0_DP*orbit%period_nd
        dt_max = tend

        cfg%sun%alpha0_rad = 0.0_DP
        cfg%sun%tu_s = TU_S
        if (mdl == 'CB') then
            cfg%model = STUDY_CANNONBALL
        else
            cfg%model = STUDY_FLAT_PLATE
        end if

        call cpu_time(wall_t0)
        call study_ads_propagate(orbit%x0, cfg, domain, index_map, da_order, err_toll_vec, &
            N_SPLIT_MAX, tend, N_SEGMENTS, REL_TOL, ABS_TOL, DT_MIN, dt_max, MAX_STEPS, &
            result_manifold, log, .false.)
        call cpu_time(wall_t1)

        nvar = dace_max_variables()
        call generate_samples(nvar, xi, n_samples)
        allocate(values(10, n_samples), found(n_samples))
        call mf_evaluate_points(result_manifold, xi, values, found, status)
        call check(status == 0, dim//' manifold evaluate ok')

        c = domain%center(); h = domain%halfwidth()
        max_rel_pos = 0.0_DP
        max_rel_vel = 0.0_DP
        n_outside = 0
        do l = 1, n_samples
            if (.not. found(l)) then
                n_outside = n_outside + 1
                cycle
            end if
            p = c
            do j = 1, 4
                if (index_map(j) > 0) p(j) = c(j) + h(j)*xi(index_map(j), l)
            end do
            call study_adaptive_integrate_real(orbit%x0, 0.0_DP, tend, REL_TOL, ABS_TOL, &
                DT_MIN, dt_max, MAX_STEPS, cfg, p, xout)
            dr_norm = norm2(values(1:3, l) - xout(1:3))
            dv_norm = norm2(values(4:6, l) - xout(4:6))
            r_norm = norm2(xout(1:3))
            v_norm = norm2(xout(4:6))
            max_rel_pos = max(max_rel_pos, dr_norm / r_norm)
            max_rel_vel = max(max_rel_vel, dv_norm / v_norm)
        end do

        write(*,'(A12,1X,A4,1X,A2,1X,I5,1X,ES8.1,1X,ES9.2,1X,ES9.2,1X,I8,1X,I7,1X,F9.2,1X,I7,1X,F9.2,1X,I7,1X,F7.2,1X,I7)') &
            trim(orbit%label), dim, mdl, da_order, toll, max_rel_pos, max_rel_vel, &
            result_manifold%n_patches, log%n_splits_total, &
            log%first_split_time*TU_S/86400.0_DP, log%n_max_active, &
            log%first_divergence_time*TU_S/86400.0_DP, log%n_diverged, wall_t1 - wall_t0, n_outside

        call mf_destroy(result_manifold)
        deallocate(xi, values, found)
    end subroutine run_case

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

end program test_crtbp_srp_study_convergence
