program test_crtbp_flatplate
    use pod_global, only: DP
    use pod_dace_classes
    use pod_crtbp_module, only: set_crtbp_mu, crtbp_derivatives_real
    use pod_bicircular_sun_module, only: bicircular_sun_position
    use pod_crtbp_srp_module, only: CrtbpSrpConfig, srp_acceleration_real
    use pod_srp_study_parameters_module
    use pod_crtbp_srp_study_module
    use pod_spacecraft_geometry, only: srp_optical_properties_type, compute_panel_force_real
    implicit none
    type(SrpStudyDomain) :: domain
    type(CrtbpSrpStudyConfig) :: cfg, cb
    type(CrtbpSrpConfig) :: legacy
    type(SrpStudyTempPool) :: pool
    type(DA) :: pd(4), rd(3), ad(3), si_da(3), xd(6), fd(6)
    real(DP) :: p(4), r(3), a(3), b(3), si(3), sun(3), pressure_nd, alpha, rs, rf, absorb
    real(DP) :: original(4), angle, x(6), f(6), g(6), h(4), eps, ap(3), am(3), expected(3)
    integer :: i,j,k,indices(4),powers(4),warm
    character(32) :: scenario
    call get_command_argument(1,scenario)
    if (scenario=='bad_domain') then
        domain=srp_study_domain(0.0_DP,0.1_DP,chi_lower=-0.1_DP)
        stop 0
    else if (scenario=='bad_tilt') then
        domain=srp_study_domain(0.0_DP,0.1_DP,tilt_max_deg=90.0_DP)
        stop 0
    else if (scenario=='polar_basis') then
        call bicircular_sun_position(0.0_DP,cfg%sun,sun)
        call srp_study_acceleration_real(0.0_DP,sun+[0.0_DP,0.0_DP,1.0_DP], &
            [0.1_DP,0.6_DP,0.1_DP,0.0_DP],cfg,a)
        stop 0
    else if (scenario=='wrong_dimension') then
        call dace_initialize(4,4)
        domain=srp_study_domain(0.0_DP,0.1_DP)
        call srp_study_init_da(domain,pd)
        stop 0
    end if
    call set_crtbp_mu(0.012153614091891635_DP)
    call require(cfg%sun%tu_s==375237.080249162_DP,'catalog TU default')
    call require(abs(srp_face_on_beta(80.0_DP,1200.0_DP)-0.1_DP)<1e-16_DP,'80 m2 reference')
    call require(abs(srp_face_on_beta(160.0_DP,1200.0_DP)-0.2_DP)<1e-16_DP,'160 m2 reference')
    call check_existing_panel()
    do k=0,10
        call srp_study_optics(real(k,DP)/10,rs,rf,absorb)
        call require(min(rs,rf,absorb)>=0.0_DP,'physical reflection fractions')
        call require(abs(rs+rf+absorb-1.0_DP)<1e-15_DP,'energy fractions sum to one')
        call require(abs(1+rs+2*rf/3-1.5_DP)<1e-15_DP,'constant face-on coefficient')
    end do
    domain=srp_study_domain(0.0_DP,0.1_DP)
    call require(domain%n_variables()==1,'Q1 one variable')
    call require(all(domain%center()==[0.05_DP,0.6_DP,0.0_DP,0.0_DP]),'midpoint differs from zero reference')
    domain=srp_study_domain(0.08_DP,0.12_DP,tilt_max_deg=10.0_DP)
    call require(domain%n_variables()==3,'Q2 three variables')
    h=domain%halfwidth()
    call require(abs(atan(sqrt(h(3)**2+h(4)**2))*180/acos(-1.0_DP)-10.0_DP)<1e-13_DP,'box tilt bound')
    cb=cfg; cb%model=STUDY_CANNONBALL
    legacy%sun=cfg%sun; legacy%cr0=1.0_DP; legacy%area_mass=0.1_DP
    r=[0.81_DP,0.12_DP,0.03_DP]
    p=[0.1_DP,0.6_DP,0.0_DP,0.0_DP]
    call srp_study_acceleration(0.7_DP,r,p,cb,a)
    call srp_acceleration_real(0.7_DP,r,[0.0_DP,0.0_DP,0.0_DP],legacy,b)
    call close(a,b,1e-18_DP,'absolute-beta CB matches existing complete-distance force')
    do k=0,10
        p(2)=real(k,DP)/10
        call srp_study_acceleration_real(0.7_DP,r,p,cfg,b)
        call close(a,b,1e-18_DP,'face-on FP equals CB for every chi')
    end do
    p=[0.0_DP,0.6_DP,0.1_DP,-0.2_DP]
    call srp_study_acceleration_real(0.7_DP,r,p,cfg,a)
    call require(all(a==0.0_DP),'zero beta exact zero force')
    p(1)=0.1_DP
    original=p
    call srp_study_acceleration_real(0.7_DP,r,p,cfg,a,si)
    call close(si,a*cfg%sun%du_m/cfg%sun%tu_s**2,1e-20_DP,'SI diagnostics')
    angle=0.37_DP; cfg%tangent_rotation_rad=angle
    p(3)=cos(angle)*original(3)+sin(angle)*original(4)
    p(4)=-sin(angle)*original(3)+cos(angle)*original(4)
    call srp_study_acceleration_real(0.7_DP,r,p,cfg,b)
    call close(a,b,1e-18_DP,'basis rotation with inverse parameter transform')
    cfg%tangent_rotation_rad=0.0_DP
    p=original
    x=[r,0.0_DP,0.1_DP,0.0_DP]
    call crtbp_derivatives_real(x,g,0.7_DP)
    call crtbp_srp_study_derivatives_real(x,f,0.7_DP,cfg,p)
    call close(f(1:3),g(1:3),0.0_DP,'kinematic RHS unchanged')
    call close(f(4:6)-g(4:6),a,1e-16_DP,'force added exactly once')
    cfg%enabled=.false.
    call crtbp_srp_study_derivatives_real(x,f,0.7_DP,cfg,p)
    call close(f,g,0.0_DP,'off recovers original CRTBP')
    cfg%enabled=.true.
    print *, 'PASS: real plate, optics, parameter domains, basis and CRTBP composition'

    ! Q4: four independent physical parameters, deterministic initial state.
    call dace_initialize(6,4)
    call pool%init()
    domain=srp_study_domain(0.08_DP,0.12_DP,tilt_max_deg=10.0_DP,chi_lower=0.4_DP,chi_upper=0.8_DP)
    call srp_study_init_da(domain,pd,indices)
    call require(all(indices==[1,2,3,4]),'Q4 physical-to-DA index map')
    h=domain%halfwidth()
    r=0.0_DP
    do i=1,3
        rd(i)=r(i)
    end do
    call srp_study_acceleration_da(0.0_DP,rd,pd,cfg,ad,pool,si_da)
    pressure_nd=cfg%pressure_1au*cfg%sun%tu_s**2/cfg%sun%du_m
    ! At tau=0: s=+x, e1=+y, e2=+z. Analytic small-tilt tests.
    call require(abs(ad(1)%cons()+pressure_nd*0.1_DP)<1e-18_DP,'DA face-on magnitude')
    call require(abs(ad(2)%get_coeff([0,0,1,0])/h(3)+pressure_nd*0.1_DP*8/15)<1e-18_DP,'analytic tilt derivative')
    call require(abs(ad(2)%get_coeff([0,1,1,0])/(h(2)*h(3))+pressure_nd*0.1_DP/3)<1e-18_DP, &
        'analytic chi-u mixed derivative')
    call require(abs(ad(3)%get_coeff([1,0,0,1])/(h(1)*h(4))+pressure_nd*8/15)<1e-18_DP, &
        'analytic beta-v mixed derivative')
    do i=1,3
        call require(abs(ad(i)%get_coeff([0,1,0,0]))<1e-18_DP,'chi alone inactive face-on')
    end do
    ! Independent central differences at a tilted state, not the same DA formula.
    p=[0.1_DP,0.6_DP,0.02_DP,-0.03_DP]
    do j=1,4
        pd(j)=p(j)+0.01_DP*da_var(j)
    end do
    call srp_study_acceleration_da(0.0_DP,rd,pd,cfg,ad,pool)
    eps=1e-5_DP
    do j=1,4
        original=p; original(j)=p(j)+eps
        call srp_study_acceleration_real(0.0_DP,r,original,cfg,ap)
        original(j)=p(j)-eps
        call srp_study_acceleration_real(0.0_DP,r,original,cfg,am)
        powers=0; powers(j)=1
        do i=1,3
            expected(i)=ad(i)%get_coeff(powers)/0.01_DP
        end do
        call close(expected,(ap-am)/(2*eps),2e-13_DP,'DA physical derivative vs finite difference')
    end do
    ! Geometry carries a propagated parameter dependence, including the tangent basis.
    rd(3)=0.03_DP+0.2_DP*da_var(1)
    call srp_study_acceleration_da(0.0_DP,rd,pd,cfg,ad,pool)
    original=p; original(1)=p(1)+0.01_DP*eps
    call srp_study_acceleration_real(0.0_DP,[0.0_DP,0.0_DP,0.03_DP+0.2_DP*eps],original,cfg,ap)
    original(1)=p(1)-0.01_DP*eps
    call srp_study_acceleration_real(0.0_DP,[0.0_DP,0.0_DP,0.03_DP-0.2_DP*eps],original,cfg,am)
    do i=1,3
        expected(i)=ad(i)%get_coeff([1,0,0,0])
    end do
    call close(expected,(ap-am)/(2*eps),2e-13_DP,'DA state-dependent geometry preserved')
    do i=1,6
        xd(i)=x(i)
    end do
    call crtbp_srp_study_derivatives_da(xd,fd,0.7_DP,cfg,pd,pool)
    warm=active_da_count()
    do k=1,100
        call crtbp_srp_study_derivatives_da(xd,fd,0.7_DP,cfg,pd,pool)
    end do
    call require(active_da_count()==warm,'RHS pool handle count stable')
    call cleanup()
    call require(active_da_count()==0,'four-parameter handles released')

    call dace_initialize(6,3)
    domain=srp_study_domain(0.0_DP,0.2_DP,tilt_max_deg=10.0_DP)
    call srp_study_init_da(domain,pd,indices)
    call require(all(indices==[1,0,2,3]),'Q2 split labels beta/u/v')
    call require(abs(pd(2)%cons()-0.6_DP)<1e-16_DP,'inactive chi has physical value 0.6')
    call cleanup()
    call dace_initialize(6,1)
    domain=srp_study_domain(0.0_DP,0.1_DP)
    call srp_study_init_da(domain,pd,indices)
    call require(all(indices==[1,0,0,0]),'Q1 split label beta')
    call require(abs(pd(1)%eval([-1.0_DP]))<1e-16_DP,'zero-beta endpoint')
    call pool%init()
    do i=1,3
        rd(i)=0.0_DP
    end do
    ! Zero constant beta must NOT skip its nonzero DA derivative.
    pd(1)=0.01_DP*da_var(1)
    call srp_study_acceleration_da(0.0_DP,rd,pd,cfg,ad,pool)
    call require(abs(ad(1)%get_coeff([1])+pressure_nd*0.01_DP)<1e-18_DP,'zero-center beta derivative survives')
    call cleanup()
    call require(active_da_count()==0,'all test handles released')
    print *, 'PASS: DA one/three/four parameters, mixed derivatives and resource lifecycle'
contains
    ! Cross-check the new effective-beta model against the independent panel API
    ! in SI force units. Its general optical inputs also cover pure-material limits
    ! that intentionally lie outside the constant-Ceq chi family used in this study.
    subroutine check_existing_panel()
        type(srp_optical_properties_type)::optical
        real(DP)::direction(3),normal(3),force(3),actual(3),unused(3),pars(4),tilt
        integer::m
        direction=[1.0_DP,0.0_DP,0.0_DP]
        tilt=0.4_DP; normal=[cos(tilt),sin(tilt),0.0_DP]
        optical%absorptivity=1.0_DP
        optical%specular_reflectivity=0.0_DP; optical%diffuse_reflectivity=0.0_DP
        call compute_panel_force_real(direction,normal,80.0_DP,cfg%pressure_1au,optical,force)
        call close(force,-cfg%pressure_1au*80*cos(tilt)*direction,1e-18_DP,'pure absorber limit')
        optical%absorptivity=0.0_DP; optical%specular_reflectivity=1.0_DP
        call compute_panel_force_real(direction,normal,80.0_DP,cfg%pressure_1au,optical,force)
        call close(force,-2*cfg%pressure_1au*80*cos(tilt)**2*normal,1e-18_DP,'pure specular limit')
        optical%specular_reflectivity=0.0_DP; optical%diffuse_reflectivity=1.0_DP
        call compute_panel_force_real(direction,direction,80.0_DP,cfg%pressure_1au,optical,force)
        call close(force,-(5.0_DP/3)*cfg%pressure_1au*80*direction,1e-18_DP,'pure Lambert face-on limit')
        call compute_panel_force_real(direction,-direction,80.0_DP,cfg%pressure_1au,optical,force)
        call require(all(force==0.0_DP),'general panel back face is unlit')
        do m=0,10
            pars=[0.1_DP,real(m,DP)/10,tan(tilt),0.0_DP]
            call srp_study_optics(pars(2),optical%specular_reflectivity,optical%diffuse_reflectivity, &
                optical%absorptivity)
            call compute_panel_force_real(direction,normal,80.0_DP,cfg%pressure_1au,optical,force)
            call srp_study_acceleration_real(0.0_DP,[0.0_DP,0.0_DP,0.0_DP],pars,cfg,unused,actual)
            call close(actual,force/1200.0_DP,1e-20_DP,'tilted study plate agrees with SI panel force/mass')
        end do
    end subroutine

    subroutine require(ok,label)
        logical,intent(in)::ok
        character(*),intent(in)::label
        if (.not.ok) then
            print *, 'FAIL: '//label
            error stop 1
        end if
    end subroutine
    subroutine close(actual,expected,tol,label)
        real(DP),intent(in)::actual(:),expected(:),tol
        character(*),intent(in)::label
        call require(maxval(abs(actual-expected))<=tol,label)
    end subroutine
    subroutine cleanup()
        integer :: m
        do m=1,4
            call pd(m)%destroy()
        end do
        do m=1,3
            call rd(m)%destroy(); call ad(m)%destroy(); call si_da(m)%destroy()
        end do
        do m=1,6
            call xd(m)%destroy(); call fd(m)%destroy()
        end do
        call pool%destroy()
    end subroutine
end program
