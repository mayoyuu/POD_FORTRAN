! Test-only fixed-step RK4: production integrators and ADS are not changed.
! Small parameter boxes test DA plumbing, not convergence over every study box.
program test_crtbp_flatplate_propagation
    use pod_global, only: DP
    use pod_dace_classes
    use pod_crtbp_module, only: set_crtbp_mu
    use pod_srp_study_parameters_module
    use pod_crtbp_srp_study_module
    implicit none
    integer, parameter :: nsteps=4096, nsamples=5
    real(DP), parameter :: orbit_period=2.3017923284002024_DP
    real(DP) :: x0(6),tend,dt
    type(CrtbpSrpStudyConfig) :: cfg
    call set_crtbp_mu(0.012153614091891635_DP)
    cfg%sun%alpha0_rad=0.37_DP
    ! Existing DRO seed; period is used only to set the regression duration.
    ! Adding SRP does not make this a newly targeted periodic orbit.
    x0=[1.1309265107780351_DP,0.0_DP,0.0_DP,0.0_DP,-0.46540743845849059_DP,0.0_DP]
    tend=2*orbit_period
    if (tend*cfg%sun%tu_s<15*86400.0_DP) tend=3*orbit_period
    dt=tend/nsteps
    call run_real_regression()
    call run_da_case(3)
    call run_da_case(4)
    print '(A,F10.5)', 'Propagation duration (days): ',tend*cfg%sun%tu_s/86400.0_DP
    print *, 'PASS: CRTBP flat-plate trajectory and parameter-only DA regressions'
contains
    subroutine require(ok,label)
        logical,intent(in)::ok
        character(*),intent(in)::label
        if (.not.ok) then
            print *, 'FAIL: '//label
            error stop 1
        end if
    end subroutine

    subroutine step_real(x,t,h,p,config)
        real(DP),intent(inout)::x(6)
        real(DP),intent(in)::t,h,p(4)
        type(CrtbpSrpStudyConfig),intent(in)::config
        real(DP)::k(6,4)
        call crtbp_srp_study_derivatives_real(x,k(:,1),t,config,p)
        call crtbp_srp_study_derivatives_real(x+h*k(:,1)/2,k(:,2),t+h/2,config,p)
        call crtbp_srp_study_derivatives_real(x+h*k(:,2)/2,k(:,3),t+h/2,config,p)
        call crtbp_srp_study_derivatives_real(x+h*k(:,3),k(:,4),t+h,config,p)
        x=x+h/6*(k(:,1)+2*k(:,2)+2*k(:,3)+k(:,4))
    end subroutine

    subroutine run_real_regression()
        type(CrtbpSrpStudyConfig)::cb,off,rotated
        real(DP)::baseline(6),zero(6),ball(6),plate(6,3),tilted(6),rotated_state(6)
        real(DP)::p(4),p_rotated(4),faceon_error,zero_error,rotation_error,angle,chi(3),t
        integer::k,j
        cb=cfg; cb%model=STUDY_CANNONBALL
        off=cfg; off%enabled=.false.
        angle=0.63_DP; rotated=cfg; rotated%tangent_rotation_rad=angle
        p=[0.1_DP,0.6_DP,0.1_DP,-0.2_DP]; p_rotated=p
        p_rotated(3)=cos(angle)*p(3)+sin(angle)*p(4)
        p_rotated(4)=-sin(angle)*p(3)+cos(angle)*p(4)
        chi=[0.0_DP,0.6_DP,1.0_DP]
        baseline=x0; zero=x0; ball=x0; tilted=x0; rotated_state=x0
        do j=1,3
            plate(:,j)=x0
        end do
        faceon_error=0.0_DP; zero_error=0.0_DP; rotation_error=0.0_DP
        do k=1,nsteps
            t=(k-1)*dt
            call step_real(baseline,t,dt,p,off)
            call step_real(zero,t,dt,[0.0_DP,0.6_DP,0.1_DP,-0.2_DP],cfg)
            call step_real(ball,t,dt,[0.1_DP,0.6_DP,0.0_DP,0.0_DP],cb)
            do j=1,3
                call step_real(plate(:,j),t,dt,[0.1_DP,chi(j),0.0_DP,0.0_DP],cfg)
                faceon_error=max(faceon_error,maxval(abs(plate(:,j)-ball)))
            end do
            call step_real(tilted,t,dt,p,cfg)
            call step_real(rotated_state,t,dt,p_rotated,rotated)
            zero_error=max(zero_error,maxval(abs(zero-baseline)))
            rotation_error=max(rotation_error,maxval(abs(tilted-rotated_state)))
        end do
        call require(zero_error==0.0_DP,'beta=0 reproduces disabled CRTBP throughout arc')
        call require(faceon_error<2e-12_DP,'face-on CB/FP and chi invariance throughout arc')
        call require(rotation_error<2e-12_DP,'same physical plate in rotated tangent basis throughout arc')
        print '(A,ES12.4)', 'Full-arc CB/FP maximum absolute difference (ND): ',faceon_error
        print '(A,ES12.4)', 'Full-arc basis-rotation maximum difference (ND): ',rotation_error
    end subroutine

    subroutine run_da_case(nvar)
        integer,intent(in)::nvar
        type(SrpStudyDomain)::domain
        type(SrpStudyTempPool)::pool
        type(DA)::state(6),stages(6,4),tmp(6),params(4)
        real(DP)::samples(nvar,nsamples),physical(4,nsamples),direct(6,nsamples),refined(6)
        real(DP)::mapped(6),c(4),w(4),error,relative_error,scaled_error,refinement_error,t
        integer::i,j,k,indices(4),live
        ! Narrow tilt permits a sixth-order polynomial to be checked at all corners.
        ! Large [0,.1]/[0,.2] study domains will require their own ADS validation.
        if (nvar==3) then
            domain=srp_study_domain(0.08_DP,0.12_DP,tilt_max_deg=3.0_DP)
        else
            domain=srp_study_domain(0.08_DP,0.12_DP,tilt_max_deg=3.0_DP,chi_lower=0.4_DP,chi_upper=0.8_DP)
        end if
        call dace_initialize(6,nvar)
        call pool%init()
        call srp_study_init_da(domain,params,indices)
        c=domain%center(); w=domain%halfwidth()
        samples(:,1)=0.0_DP
        do i=1,nvar
            samples(i,2)=(-1.0_DP)**i*0.2_DP
        end do
        samples(:,3)=1.0_DP; samples(:,4)=-1.0_DP
        samples(:,5)=1.0_DP; samples(1,5)=-1.0_DP
        do j=1,nsamples
            physical(:,j)=c
            do i=1,4
                if (indices(i)>0) physical(i,j)=c(i)+w(i)*samples(indices(i),j)
            end do
            direct(:,j)=x0
        end do
        do i=1,6
            state(i)=x0(i) ! No DA initial-position/velocity uncertainty.
        end do
        do k=1,nsteps
            t=(k-1)*dt
            do j=1,nsamples
                call step_real(direct(:,j),t,dt,physical(:,j),cfg)
            end do
            call crtbp_srp_study_derivatives_da(state,stages(:,1),t,cfg,params,pool)
            do i=1,6
                tmp(i)=state(i)+(dt/2)*stages(i,1)
            end do
            call crtbp_srp_study_derivatives_da(tmp,stages(:,2),t+dt/2,cfg,params,pool)
            do i=1,6
                tmp(i)=state(i)+(dt/2)*stages(i,2)
            end do
            call crtbp_srp_study_derivatives_da(tmp,stages(:,3),t+dt/2,cfg,params,pool)
            do i=1,6
                tmp(i)=state(i)+dt*stages(i,3)
            end do
            call crtbp_srp_study_derivatives_da(tmp,stages(:,4),t+dt,cfg,params,pool)
            do i=1,6
                state(i)=state(i)+(dt/6)*(stages(i,1)+2.0_DP*stages(i,2)+2.0_DP*stages(i,3)+stages(i,4))
            end do
            if (k==1) live=active_da_count()
            call require(active_da_count()==live,'propagation DA handles stable')
        end do
        error=0.0_DP; relative_error=0.0_DP; scaled_error=0.0_DP
        do j=1,nsamples
            do i=1,6
                mapped(i)=state(i)%eval(samples(:,j))
            end do
            if (j==1) call require(maxval(abs(mapped-direct(:,j)))<2e-11_DP,'DA center equals real trajectory')
            error=max(error,maxval(abs(mapped-direct(:,j))))
            relative_error=max(relative_error,norm2(mapped-direct(:,j))/norm2(direct(:,j)))
            ! Position and velocity blocks are scaled separately, with a floor near zero.
            scaled_error=max(scaled_error,norm2(mapped(1:3)-direct(1:3,j))/max(norm2(direct(1:3,j)),0.1_DP), &
                norm2(mapped(4:6)-direct(4:6,j))/max(norm2(direct(4:6,j)),0.1_DP))
        end do
        refined=x0
        do k=1,2*nsteps
            call step_real(refined,(k-1)*dt/2,dt/2,physical(:,3),cfg)
        end do
        refinement_error=maxval(abs(refined-direct(:,3)))
        print '(A,I0,A,ES12.4)', 'DA ',nvar,'D endpoint absolute error (ND): ',error
        print '(A,ES12.4)', 'Endpoint vector relative error:              ',relative_error
        print '(A,ES12.4)', 'Position/velocity block-scaled error:        ',scaled_error
        print '(A,ES12.4)', 'Real corner step-halving difference (ND):    ',refinement_error
        call require(error<5e-9_DP,'sixth-order DA center/interior/corner propagation')
        call require(refinement_error<2e-9_DP,'independent RK4 time-step refinement')
        do i=1,6
            call state(i)%destroy(); call tmp(i)%destroy()
            do j=1,4
                call stages(i,j)%destroy()
            end do
        end do
        do i=1,4
            call params(i)%destroy()
        end do
        call pool%destroy()
        call require(active_da_count()==0,'all integration handles released before reinitializing DACE')
    end subroutine
end program
