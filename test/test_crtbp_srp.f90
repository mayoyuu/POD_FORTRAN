program test_crtbp_srp
    use pod_global, only: DP
    use pod_dace_classes
    use pod_crtbp_module, only: crtbp_derivatives_real, crtbp_derivatives_da
    use pod_bicircular_sun_module
    use pod_crtbp_srp_module
    implicit none
    type(CrtbpSrpConfig) :: cfg
    type(CrtbpSrpForceTempPool) :: pool
    type(DA) :: xd(6), fd(6), gd(6), eta_da(3), ad(3), asi_da(3)
    real(DP) :: sun(3), sun2(3), a(3), b(3), asi(3), r(3), eta(3), f(6), g(6), x(6)
    real(DP) :: pi, omega, gamma, period, eps, ap(3), am(3), scales(3), expected(3)
    integer :: i, j, mode, before, warm
    logical :: fallback
    character(32) :: arg

    call get_command_argument(1,arg)
    pi = acos(-1.0_DP)
    omega = 2*pi*cfg%sun%tu_s/cfg%sun%period_s - 1.0_DP
    period = 2*pi/abs(omega)
    gamma = cfg%pressure_1au * cfg%cr0 * cfg%area_mass * cfg%sun%tu_s**2/cfg%sun%du_m
    eta = 0.0_DP
    r = 0.0_DP
    if (trim(arg) == 'invalid_units') then
        cfg%sun%du_m = 0.0_DP
        call bicircular_sun_position(0.0_DP, cfg%sun, sun)
        stop 0
    else if (trim(arg) == 'sun_collision') then
        call bicircular_sun_position(0.0_DP, cfg%sun, sun)
        call srp_acceleration_real(0.0_DP, sun, eta, cfg, a)
        stop 0
    else if (trim(arg) == 'invalid_parameters') then
        call dace_initialize(4,3)
        call srp_init_da_parameters([.true.,.false.,.false.],[0.3_DP,0.0_DP,0.0_DP],eta_da)
        stop 0
    else if (trim(arg) == 'da_sun_collision') then
        call dace_initialize(4,1)
        call pool%init()
        call bicircular_sun_position(0.0_DP,cfg%sun,sun)
        call srp_init_da_parameters([.true.,.false.,.false.],[0.3_DP,0.0_DP,0.0_DP],eta_da)
        do i=1,3
            xd(i)=sun(i)
        end do
        call srp_acceleration_da(0.0_DP,xd(1:3),eta_da,cfg,ad,pool)
        stop 0
    end if
    call require(omega < 0.0_DP, 'Sun rotates backwards in synodic frame')
    call close_vec([cfg%sun%du_m,cfg%sun%tu_s,0.0_DP], &
        [384400000.0_DP,375190.26189464843_DP,0.0_DP],0.0_DP,'existing CRTBP units')
    call require(abs(sun_elapsed_tau(110.0_DP,10.0_DP,cfg%sun)-100.0_DP/cfg%sun%tu_s)<1e-18_DP, &
        'epoch offset applied once')
    call bicircular_sun_position(0.0_DP,cfg%sun,sun)
    call bicircular_sun_position(period/2,cfg%sun,sun2)
    call close_vec(sun2,-sun,1e-12_DP,'half-period Sun reversal')
    do i=-3,3
        call bicircular_sun_position(real(i,DP)*0.7_DP,cfg%sun,sun2)
        call require(abs(norm2(sun2)-cfg%sun%radius_m/cfg%sun%du_m)<1e-12_DP,'Sun radius')
    end do
    call srp_acceleration_real(0.0_DP,r,eta,cfg,a,asi)
    call close_vec(a,[-gamma,0.0_DP,0.0_DP],1e-18_DP,'Sun +x gives force -x')
    call require(abs(norm2(asi)-cfg%pressure_1au*cfg%cr0*cfg%area_mass)<1e-20_DP,'SI amplitude')
    cfg%sun%alpha0_rad = pi/2
    call srp_acceleration_real(0.0_DP,r,eta,cfg,a)
    call close_vec(a,[0.0_DP,-gamma,0.0_DP],1e-18_DP,'Sun +y gives force -y')
    cfg%sun%alpha0_rad = 0.0_DP
    cfg%sun%radius_m = 2*cfg%sun%au_m
    call srp_acceleration_real(0.0_DP,r,eta,cfg,a)
    call close_vec(a,[-gamma/4,0.0_DP,0.0_DP],1e-18_DP,'inverse square at 2 AU')
    cfg%sun%radius_m = cfg%sun%au_m
    r = [0.8_DP,0.1_DP,0.2_DP]
    call srp_acceleration_real(0.0_DP,r,eta,cfg,a)
    call require(a(3)>0.0_DP,'3D parallax gives positive z acceleration')
    cfg%model = SRP_CANNONBALL_FARFIELD
    call srp_acceleration_real(0.0_DP,r,eta,cfg,b)
    call require(norm2(a-b)/norm2(b)<0.006_DP,'full/far-field proximity')
    call require(b(3)==0.0_DP,'far-field has zero direct z force')
    call srp_acceleration_real(period/2,r,eta,cfg,a)
    call close_vec(a,-b,1e-18_DP,'far-field half-period reversal')
    cfg%model = SRP_CANNONBALL_FULL
    eta = [0.2_DP,0.0_DP,0.0_DP]
    call srp_acceleration_real(0.3_DP,r,eta,cfg,a)
    cfg%model = SRP_EMPIRICAL_THREE_COMPONENT
    call srp_acceleration_real(0.3_DP,r,eta,cfg,b)
    call close_vec(a,b,1e-18_DP,'M2 radial limit equals M1')
    cfg%sun%alpha0_rad = pi
    eta = [0.1_DP,0.2_DP,0.3_DP]
    call srp_acceleration_real(0.0_DP,[0.0_DP,0.0_DP,0.0_DP],eta,cfg,a)
    call close_vec(a,gamma*[1.1_DP,0.2_DP,0.3_DP],1e-18_DP,'M2 signs and no renormalization')
    cfg%sun%alpha0_rad = 0.0_DP
    call bicircular_sun_position(0.0_DP,cfg%sun,sun)
    call srp_acceleration_real(0.0_DP,sun+[0.0_DP,0.0_DP,1.0_DP],eta,cfg,a, &
        used_fallback=fallback)
    call require(fallback,'polar geometry uses fallback axis')
    x = [0.8_DP,0.1_DP,0.2_DP,0.01_DP,0.12_DP,-0.02_DP]
    call crtbp_derivatives_real(x,g,0.3_DP)
    do mode=1,5
        cfg = CrtbpSrpConfig()
        select case(mode)
        case(1)
            cfg%enabled = .false.
        case(2)
            cfg%pressure_1au = 0.0_DP
        case(3)
            cfg%area_mass = 0.0_DP
        case(4)
            cfg%cr0 = 0.0_DP
        case(5)
            cfg%model = SRP_OFF
        end select
        call crtbp_srp_derivatives_real(x,f,0.3_DP,cfg,eta)
        call close_vec(f,g,0.0_DP,'zero SRP exactly recovers CRTBP')
    end do
    cfg = CrtbpSrpConfig()
    eta = 0.0_DP
    call srp_acceleration(0.0_DP,[0.0_DP,0.0_DP,0.0_DP],eta,cfg,a, &
        prescribed_sun_nd=[0.0_DP,cfg%sun%radius_m/cfg%sun%du_m,0.0_DP])
    call close_vec(a,[0.0_DP,-gamma,0.0_DP],1e-18_DP,'replaceable Sun geometry boundary')
    call srp_acceleration_real(0.3_DP,x(1:3),eta,cfg,a)
    call crtbp_srp_derivatives_real(x,f,0.3_DP,cfg,eta)
    call close_vec(f(1:3),g(1:3),0.0_DP,'kinematic RHS unchanged')
    call close_vec(f(4:6)-g(4:6),a,1e-16_DP,'SRP added once to rotating RHS')
    print *, 'PASS: real geometry, signs, units, models and RHS'

    call dace_initialize(4,1)
    before = active_da_count()
    call pool%init()
    do i=1,6
        xd(i) = x(i)
    end do
    scales = [0.3_DP,0.0_DP,0.0_DP]
    call srp_init_da_parameters([.true.,.false.,.false.],scales,eta_da)
    call require(dace_max_variables()==1,'only one DA variable')
    do i=1,6
        call require(xd(i)%get_coeff([1])==0.0_DP,'initial state constant')
    end do
    do mode=SRP_CANNONBALL_FULL,SRP_EMPIRICAL_THREE_COMPONENT
        cfg%model = mode
        call srp_acceleration_da(0.3_DP,xd(1:3),eta_da,cfg,ad,pool)
        call srp_acceleration_real(0.3_DP,x(1:3),eta,cfg,a)
        do i=1,3
            call require(abs(ad(i)%cons()-a(i))<1e-18_DP,'DA constant equals real')
            call require(abs(ad(i)%get_coeff([1])/scales(1)-a(i))<1e-18_DP,'radial sensitivity scaling')
        end do
    end do
    ! Synthetic propagated state dependence: never drop the DA part of r-r_S.
    cfg%model = SRP_CANNONBALL_FULL
    xd(3) = x(3) + 0.02_DP*da_var(1)
    call srp_acceleration_da(0.3_DP,xd(1:3),eta_da,cfg,ad,pool)
    eps = 1e-4_DP
    r = x(1:3); r(3) = r(3)+0.02_DP*eps
    call srp_acceleration_real(0.3_DP,r,[scales(1)*eps,0.0_DP,0.0_DP],cfg,ap)
    r(3) = x(3)-0.02_DP*eps
    call srp_acceleration_real(0.3_DP,r,[-scales(1)*eps,0.0_DP,0.0_DP],cfg,am)
    do i=1,3
        call require(abs(ad(i)%get_coeff([1])-(ap(i)-am(i))/(2*eps))<1e-15_DP,'DA geometry chain rule')
    end do
    call crtbp_srp_derivatives_da(xd,fd,0.3_DP,cfg,eta_da,pool)
    warm = active_da_count()
    do j=1,100
        call crtbp_srp_derivatives_da(xd,fd,0.3_DP,cfg,eta_da,pool)
    end do
    call require(active_da_count()==warm,'RHS repeated calls do not leak DA handles')
    cfg%enabled = .false.
    call crtbp_srp_derivatives_da(xd,fd,0.3_DP,cfg,eta_da,pool)
    call crtbp_derivatives_da(xd,gd,0.3_DP,pool%gravity)
    do i=1,6
        call require(abs(fd(i)%eval([0.2_DP])-gd(i)%eval([0.2_DP]))<1e-16_DP,'DA off limit')
    end do
    call cleanup()
    call require(active_da_count()==before,'pool and caller DA handles released')
    print *, 'PASS: one-parameter DA and handle lifecycle'
    call test_three_parameters()
    print *, 'PASS: all CRTBP SRP force tests'
contains
    subroutine require(ok,label)
        logical, intent(in) :: ok
        character(*), intent(in) :: label
        if (.not.ok) then
            print *, 'FAIL: '//label
            error stop 1
        end if
    end subroutine
    subroutine close_vec(actual,wanted,tol,label)
        real(DP), intent(in) :: actual(:),wanted(:),tol
        character(*), intent(in) :: label
        call require(maxval(abs(actual-wanted))<=tol,label)
    end subroutine
    subroutine cleanup()
        do i=1,6
            call xd(i)%destroy(); call fd(i)%destroy(); call gd(i)%destroy()
        end do
        do i=1,3
            call eta_da(i)%destroy(); call ad(i)%destroy(); call asi_da(i)%destroy()
        end do
        call pool%destroy()
    end subroutine
    subroutine test_three_parameters()
        call dace_initialize(4,3)
        call pool%init()
        cfg = CrtbpSrpConfig()
        cfg%model = SRP_EMPIRICAL_THREE_COMPONENT
        scales = [0.3_DP,0.1_DP,0.1_DP]
        call srp_init_da_parameters([.true.,.true.,.true.],scales,eta_da)
        do i=1,6
            xd(i) = x(i)
        end do
        call srp_acceleration_da(0.3_DP,xd(1:3),eta_da,cfg,ad,pool)
        eta = [0.02_DP,-0.03_DP,0.04_DP]
        call srp_acceleration_real(0.3_DP,x(1:3),eta,cfg,a)
        do i=1,3
            expected(i) = ad(i)%eval(eta/scales)
        end do
        call close_vec(expected,a,1e-18_DP,'three DA parameters evaluate to real M2')
        ! Verify each physical parameter sensitivity independently.
        do j=1,3
            eta=0.0_DP; eta(j)=1e-4_DP
            call srp_acceleration_real(0.3_DP,x(1:3),eta,cfg,ap)
            eta(j)=-1e-4_DP
            call srp_acceleration_real(0.3_DP,x(1:3),eta,cfg,am)
            do i=1,3
                block
                    integer :: powers(3)
                    powers=0; powers(j)=1
                    expected(i)=ad(i)%get_coeff(powers)/scales(j)
                end block
            end do
            call close_vec(expected,(ap-am)/2e-4_DP,1e-15_DP,'M2 parameter finite difference')
        end do
        call bicircular_sun_position(0.0_DP,cfg%sun,sun)
        r=sun+[0.0_DP,0.0_DP,1.0_DP]
        do i=1,3
            xd(i)=r(i)
        end do
        call srp_acceleration(0.0_DP,xd(1:3),eta_da,cfg,ad,pool,asi_da,fallback)
        call require(fallback,'DA polar fallback')
        eta=[0.02_DP,-0.03_DP,0.04_DP]
        call srp_acceleration_real(0.0_DP,r,eta,cfg,a,asi)
        do i=1,3
            expected(i)=ad(i)%eval(eta/scales)
            call require(abs(asi_da(i)%eval(eta/scales)-asi(i))<1e-15_DP,'DA SI acceleration')
        end do
        call close_vec(expected,a,1e-12_DP,'DA fallback agrees with real')
        call cleanup()
        call require(active_da_count()==0,'three-parameter handles released')
        call dace_initialize(4,1)
        call srp_init_da_parameters([.false.,.false.,.true.],[0.0_DP,0.0_DP,0.1_DP],eta_da)
        call require(eta_da(1)%eval([1.0_DP])==0.0_DP .and. eta_da(2)%eval([1.0_DP])==0.0_DP, &
            'inactive parameters remain zero')
        call require(abs(eta_da(3)%get_coeff([1])-0.1_DP)<1e-16_DP,'active subset packed into variable 1')
        do i=1,3
            call eta_da(i)%destroy()
        end do
    end subroutine
end program test_crtbp_srp
