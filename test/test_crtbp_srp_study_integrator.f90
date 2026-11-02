!> Validate the study adaptive RKF78 integrator (real + DA).
!! Real path: RK4 step-halving reference vs adaptive RKF78.
!! DA path (1D/3D/4D): DA center/constant == real, DA polynomial == real at
!! samples, and DA first derivative == central finite difference. Narrow boxes
!! keep a sixth-order DA accurate; the full [0,.1]/[0,.2] study boxes are the
!! job of the ADS engine, not this plumbing test.
program test_crtbp_srp_study_integrator
    use pod_global, only: DP
    use pod_dace_classes
    use pod_crtbp_module, only: set_crtbp_mu
    use pod_srp_study_parameters_module
    use pod_crtbp_srp_study_module
    use pod_crtbp_srp_study_integrator_module
    implicit none
    real(DP), parameter :: orbit_period = 2.3017923284002024_DP
    real(DP), parameter :: rel_tol = 1.0e-10_DP, abs_tol = 1.0e-12_DP
    real(DP), parameter :: dt_min = 1.0e-13_DP
    integer,  parameter :: max_steps = 100000
    real(DP) :: x0(6), tend, dt_max
    type(CrtbpSrpStudyConfig) :: cfg

    call set_crtbp_mu(0.012153614091891635_DP)
    cfg%sun%alpha0_rad = 0.37_DP
    x0 = [1.1309265107780351_DP, 0.0_DP, 0.0_DP, 0.0_DP, -0.46540743845849059_DP, 0.0_DP]
    tend = 2.0_DP*orbit_period
    if (tend*cfg%sun%tu_s < 15.0_DP*86400.0_DP) tend = 3.0_DP*orbit_period
    dt_max = tend

    call run_real_convergence()
    call run_da_case(1)
    call run_da_case(3)
    call run_da_case(4)

    print '(A,F10.5)', 'Propagation duration (days): ', tend*cfg%sun%tu_s/86400.0_DP
    print *, 'PASS: CRTBP SRP study integrator (real + DA RKF78)'
contains

    subroutine require(ok, label)
        logical, intent(in) :: ok
        character(*), intent(in) :: label
        if (.not. ok) then
            print *, 'FAIL: '//label
            error stop 1
        end if
    end subroutine require

    subroutine step_real(x, t, h, p, config)
        real(DP), intent(inout) :: x(6)
        real(DP), intent(in) :: t, h, p(4)
        type(CrtbpSrpStudyConfig), intent(in) :: config
        real(DP) :: k(6,4)
        call crtbp_srp_study_derivatives_real(x, k(:,1), t, config, p)
        call crtbp_srp_study_derivatives_real(x+h*k(:,1)/2.0_DP, k(:,2), t+h/2.0_DP, config, p)
        call crtbp_srp_study_derivatives_real(x+h*k(:,2)/2.0_DP, k(:,3), t+h/2.0_DP, config, p)
        call crtbp_srp_study_derivatives_real(x+h*k(:,3), k(:,4), t+h, config, p)
        x = x + h/6.0_DP*(k(:,1)+2.0_DP*k(:,2)+2.0_DP*k(:,3)+k(:,4))
    end subroutine step_real

    subroutine propagate_rk4(xstart, tend2, nsteps, p, config, xout)
        real(DP), intent(in) :: xstart(6), tend2, p(4)
        integer, intent(in) :: nsteps
        type(CrtbpSrpStudyConfig), intent(in) :: config
        real(DP), intent(out) :: xout(6)
        real(DP) :: h, t
        integer :: k
        xout = xstart
        h = tend2/nsteps
        do k = 1, nsteps
            t = (k-1)*h
            call step_real(xout, t, h, p, config)
        end do
    end subroutine propagate_rk4

    subroutine propagate_adaptive(p, config, xout)
        real(DP), intent(in) :: p(4)
        type(CrtbpSrpStudyConfig), intent(in) :: config
        real(DP), intent(out) :: xout(6)
        call study_adaptive_integrate_real(x0, 0.0_DP, tend, rel_tol, abs_tol, &
            dt_min, dt_max, max_steps, config, p, xout)
    end subroutine propagate_adaptive

    subroutine run_real_convergence()
        real(DP) :: p(4), ref_coarse(6), ref_fine(6), adapt(6), refinement
        p = [0.1_DP, 0.6_DP, 0.0_DP, 0.0_DP]
        call propagate_rk4(x0, tend, 2048, p, cfg, ref_coarse)
        call propagate_rk4(x0, tend, 4096, p, cfg, ref_fine)
        refinement = maxval(abs(ref_fine-ref_coarse))
        call propagate_adaptive(p, cfg, adapt)
        print '(A,ES12.4)', 'RK4 2048 vs 4096 step-halving difference (ND): ', refinement
        print '(A,ES12.4)', 'Adaptive RKF78 vs 4096-step RK4 (ND):            ', maxval(abs(adapt-ref_fine))
        call require(refinement < 1.0e-8_DP, 'RK4 step-halving refinement small')
        call require(maxval(abs(adapt-ref_fine)) < 1.0e-8_DP, 'adaptive RKF78 matches RK4 reference')
    end subroutine run_real_convergence

    subroutine run_da_case(nvar)
        integer, intent(in) :: nvar
        type(SrpStudyDomain) :: domain
        type(StudyIntegratorTempPool) :: ipool
        type(DA) :: params(4)
        type(AlgebraicVector) :: dstate, dfinal
        integer :: index_map(4)
        real(DP) :: c(4), w(4), physical(4,5), direct(6,5), center_real(6)
        real(DP) :: mapped(6), da_center(6), err, fd(6), da_deriv
        real(DP) :: p_plus(4), p_minus(4), x_plus(6), x_minus(6)
        real(DP), allocatable :: samples(:,:)
        real(DP), parameter :: delta_xi = 1.0e-4_DP
        integer :: i, j, k, l

        if (nvar == 1) then
            domain = srp_study_domain(0.08_DP, 0.12_DP)
        else if (nvar == 3) then
            domain = srp_study_domain(0.08_DP, 0.12_DP, tilt_max_deg=3.0_DP)
        else
            domain = srp_study_domain(0.08_DP, 0.12_DP, tilt_max_deg=3.0_DP, &
                chi_lower=0.4_DP, chi_upper=0.8_DP)
        end if

        call dace_initialize(6, nvar)
        call ipool%init()
        call srp_study_init_da(domain, params, index_map)
        c = domain%center(); w = domain%halfwidth()

        allocate(samples(nvar, 5))
        samples = 0.0_DP
        do i = 1, nvar
            samples(i,2) = 0.4_DP
            samples(i,3) = -0.4_DP
        end do
        samples(:,4) = 1.0_DP
        samples(:,5) = -1.0_DP
        do l = 1, 5
            physical(:,l) = c
            do j = 1, 4
                if (index_map(j) > 0) physical(j,l) = c(j) + w(j)*samples(index_map(j),l)
            end do
        end do

        call dstate%init(6)
        do i = 1, 6
            dstate%elements(i) = x0(i)
        end do
        call dfinal%init(6)
        call study_adaptive_integrate_da(dstate, 0.0_DP, tend, rel_tol, abs_tol, &
            dt_min, dt_max, max_steps, cfg, params, ipool, dfinal)

        ! DA constant term == real at center
        call propagate_adaptive(c, cfg, center_real)
        do i = 1, 6
            da_center(i) = dfinal%elements(i)%cons()
        end do
        err = maxval(abs(da_center-center_real))
        print '(A,I0,A,ES12.4)', 'DA ', nvar, 'D center vs real (ND): ', err
        call require(err < 1.0e-8_DP, 'DA constant term equals real center propagation')

        ! DA polynomial == real at samples
        err = 0.0_DP
        do l = 1, 5
            call propagate_adaptive(physical(:,l), cfg, direct(:,l))
            do i = 1, 6
                mapped(i) = dfinal%elements(i)%eval(samples(:,l))
            end do
            err = max(err, maxval(abs(mapped-direct(:,l))))
        end do
        print '(A,I0,A,ES12.4)', 'DA ', nvar, 'D sample vs real (ND): ', err
        call require(err < 1.0e-7_DP, 'DA polynomial matches real propagation at samples')

        ! DA first derivative == central finite difference (d/dxi)
        err = 0.0_DP
        do j = 1, 4
            if (index_map(j) == 0) cycle
            k = index_map(j)
            p_plus = c;  p_plus(j)  = c(j) + w(j)*delta_xi
            p_minus = c; p_minus(j) = c(j) - w(j)*delta_xi
            call propagate_adaptive(p_plus, cfg, x_plus)
            call propagate_adaptive(p_minus, cfg, x_minus)
            do i = 1, 6
                fd(i) = (x_plus(i)-x_minus(i))/(2.0_DP*delta_xi)
                da_deriv = dfinal%elements(i)%get_deriv_value(k)
                err = max(err, abs(fd(i)-da_deriv)/max(1.0_DP, abs(da_deriv)))
            end do
        end do
        print '(A,I0,A,ES12.4)', 'DA ', nvar, 'D derivative vs FD (rel): ', err
        call require(err < 1.0e-5_DP, 'DA first derivative matches central finite difference')

        deallocate(samples)
        call dstate%destroy()
        call dfinal%destroy()
        do i = 1, 4
            call params(i)%destroy()
        end do
        call ipool%destroy()
        call require(active_da_count() == 0, 'all DA handles released')
    end subroutine run_da_case
end program test_crtbp_srp_study_integrator
