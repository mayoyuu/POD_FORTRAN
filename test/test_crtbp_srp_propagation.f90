! Test-only integration harness. Does not change production integrators or ADS.
program test_crtbp_srp_propagation
    use pod_global, only: DP
    use pod_dace_classes
    use pod_crtbp_module, only: crtbp_mu
    use pod_crtbp_srp_module
    use pod_bicircular_sun_module, only: bicircular_sun_position, bicircular_sun_phase
    implicit none
    integer, parameter :: nsteps=2048
    type(CrtbpSrpConfig) :: cfg,off_cfg
    type(CrtbpSrpForceTempPool) :: pool
    type(DA) :: state(6),kd(6,4),tmp(6),params(3)
    real(DP) :: x0(6),nominal(6),off(6),direct(6,4),refined(6),mapped(6),sensitivity(6)
    real(DP) :: endpoint_error,refinement_error,fd_error,h,t,tend,eta_values(4),sun(3),a(3),scale
    integer :: i,j,k,unit,ios
    character(512) :: csv_path
    ! Existing seed from test_uq_crtbp_comparison.f90. No claim of periodicity.
    x0=[1.1309265107780351_DP,0.0_DP,0.0_DP,0.0_DP,-0.46540743845849059_DP,0.0_DP]
    tend=15.0_DP*86400.0_DP/cfg%sun%tu_s
    h=tend/nsteps
    scale=0.30_DP
    eta_values=[-0.05_DP,0.05_DP,-0.001_DP,0.001_DP]
    off_cfg=cfg; off_cfg%enabled=.false.
    nominal=x0; off=x0
    do j=1,4
        direct(:,j)=x0
    end do
    call dace_initialize(4,1)
    call pool%init()
    call srp_init_da_parameters([.true.,.false.,.false.],[scale,0.0_DP,0.0_DP],params)
    do i=1,6
        state(i)=x0(i)
        call require(state(i)%get_coeff([1])==0.0_DP,'fixed initial state')
    end do
    call get_command_argument(1,csv_path)
    if (len_trim(csv_path)>0) then
        open(newunit=unit,file=trim(csv_path),status='replace',action='write',iostat=ios)
        call require(ios==0,'open diagnostic CSV')
        write(unit,'(A)') 'tau,t_seconds,alpha_s,sun_x,sun_y,sun_z,srp_ax_nd,srp_ay_nd,srp_az_nd,'// &
            'nom_x,nom_y,nom_z,nom_vx,nom_vy,nom_vz,off_x,off_y,off_z,off_vx,off_vy,off_vz,'// &
            'minus_x,minus_y,minus_z,minus_vx,minus_vy,minus_vz,plus_x,plus_y,plus_z,plus_vx,plus_vy,plus_vz'
        call write_row(0.0_DP)
    end if
    do k=1,nsteps
        t=(k-1)*h
        call step_real(nominal,t,h,0.0_DP,cfg)
        call step_real(off,t,h,0.0_DP,off_cfg)
        do j=1,4
            call step_real(direct(:,j),t,h,eta_values(j),cfg)
        end do
        call step_da(t,h)
        if (len_trim(csv_path)>0 .and. (mod(k,32)==0 .or. k==nsteps)) call write_row(t+h)
    end do
    if (len_trim(csv_path)>0) close(unit)
    do i=1,6
        mapped(i)=state(i)%cons()
        sensitivity(i)=state(i)%get_coeff([1])/scale
    end do
    call require(maxval(abs(mapped-nominal))<2e-11_DP,'DA nominal equals real nominal after 15 days')
    endpoint_error=0.0_DP
    do j=1,2
        do i=1,6
            mapped(i)=state(i)%eval([eta_values(j)/scale])
        end do
        endpoint_error=max(endpoint_error,maxval(abs(mapped-direct(:,j))))
    end do
    call require(endpoint_error<1e-10_DP,'fourth-order DA endpoints at eta=+/-0.05')
    fd_error=maxval(abs(sensitivity-(direct(:,4)-direct(:,3))/0.002_DP))
    call require(fd_error<1e-8_DP,'physical eta sensitivity vs central finite difference')
    refined=x0
    do k=1,2*nsteps
        call step_real(refined,(k-1)*h/2,h/2,0.0_DP,cfg)
    end do
    refinement_error=maxval(abs(refined-nominal))
    call require(refinement_error<2e-9_DP,'time step refinement')
    call require(abs(jacobi(off)-jacobi(x0))<1e-9_DP,'SRP off Jacobi conservation')
    call require(maxval(abs(nominal-off))>1e-7_DP,'nominal SRP produces an orbit change')
    print '(A,ES12.4)', '15-day DA/real endpoint max error (ND): ',endpoint_error
    print '(A,ES12.4)', 'Real step-halving max difference (ND): ',refinement_error
    print '(A,ES12.4)', 'dX/deta finite-difference max error:   ',fd_error
    print '(A,ES12.4)', 'SRP off Jacobi drift:                 ',abs(jacobi(off)-jacobi(x0))
    print '(A,ES12.4)', 'Nominal SRP position effect (m):      ',norm2(nominal(1:3)-off(1:3))*cfg%sun%du_m
    print '(A,I0,A,I0,A,F6.3)', 'DA order=',4,', variables=',dace_max_variables(),', eta half-width=',scale
    do i=1,6
        call state(i)%destroy(); call tmp(i)%destroy()
        do j=1,4
            call kd(i,j)%destroy()
        end do
    end do
    do i=1,3
        call params(i)%destroy()
    end do
    call pool%destroy()
    call require(active_da_count()==0,'all integration handles released')
    print *, 'PASS: 15-day CRTBP SRP parameter-only DA propagation'
contains
    subroutine require(ok,label)
        logical, intent(in) :: ok
        character(*), intent(in) :: label
        if (.not.ok) then
            print *, 'FAIL: '//label
            error stop 1
        end if
    end subroutine
    subroutine step_real(x,t,dt,eta_r,config)
        real(DP), intent(inout) :: x(6)
        real(DP), intent(in) :: t,dt,eta_r
        type(CrtbpSrpConfig), intent(in) :: config
        real(DP) :: stages(6,4),eta(3)
        eta=[eta_r,0.0_DP,0.0_DP]
        call crtbp_srp_derivatives_real(x,stages(:,1),t,config,eta)
        call crtbp_srp_derivatives_real(x+dt*stages(:,1)/2,stages(:,2),t+dt/2,config,eta)
        call crtbp_srp_derivatives_real(x+dt*stages(:,2)/2,stages(:,3),t+dt/2,config,eta)
        call crtbp_srp_derivatives_real(x+dt*stages(:,3),stages(:,4),t+dt,config,eta)
        x=x+dt/6*(stages(:,1)+2*stages(:,2)+2*stages(:,3)+stages(:,4))
    end subroutine
    subroutine step_da(t,dt)
        real(DP), intent(in) :: t,dt
        integer :: m
        call crtbp_srp_derivatives_da(state,kd(:,1),t,cfg,params,pool)
        do m=1,6
            tmp(m)=state(m)+(dt/2)*kd(m,1)
        end do
        call crtbp_srp_derivatives_da(tmp,kd(:,2),t+dt/2,cfg,params,pool)
        do m=1,6
            tmp(m)=state(m)+(dt/2)*kd(m,2)
        end do
        call crtbp_srp_derivatives_da(tmp,kd(:,3),t+dt/2,cfg,params,pool)
        do m=1,6
            tmp(m)=state(m)+dt*kd(m,3)
        end do
        call crtbp_srp_derivatives_da(tmp,kd(:,4),t+dt,cfg,params,pool)
        do m=1,6
            state(m)=state(m)+(dt/6)*(kd(m,1)+2.0_DP*kd(m,2)+2.0_DP*kd(m,3)+kd(m,4))
        end do
    end subroutine
    function jacobi(x) result(c)
        real(DP), intent(in) :: x(6)
        real(DP) :: c,r1,r2
        r1=norm2(x(1:3)+[crtbp_mu,0.0_DP,0.0_DP])
        r2=norm2(x(1:3)+[crtbp_mu-1.0_DP,0.0_DP,0.0_DP])
        c=x(1)**2+x(2)**2+2*(1-crtbp_mu)/r1+2*crtbp_mu/r2-sum(x(4:6)**2)
    end function
    subroutine write_row(tau)
        real(DP), intent(in) :: tau
        call bicircular_sun_position(tau,cfg%sun,sun)
        call srp_acceleration_real(tau,nominal(1:3),[0.0_DP,0.0_DP,0.0_DP],cfg,a)
        write(unit,'(*(ES24.16,:,","))') tau,tau*cfg%sun%tu_s,bicircular_sun_phase(tau,cfg%sun), &
            sun,a,nominal,off,direct(:,1),direct(:,2)
    end subroutine
end program test_crtbp_srp_propagation
