! CRTBP + prescribed bicircular Sun geometry + SRP, all in synodic ND units.
! No solar gravity, shadow, state uncertainty initialization or integration here.
module pod_crtbp_srp_module
    use pod_global, only: DP
    use pod_dace_classes
    use pod_crtbp_module, only: CrtbpForceTempPool, crtbp_derivatives_real, crtbp_derivatives_da
    use pod_bicircular_sun_module
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private
    integer, parameter, public :: SRP_OFF=0, SRP_CANNONBALL_FULL=1
    integer, parameter, public :: SRP_CANNONBALL_FARFIELD=2, SRP_EMPIRICAL_THREE_COMPONENT=3
    real(DP), parameter :: min_distance_nd=1e-12_DP, min_transverse_squared=1e-20_DP

    type, public :: CrtbpSrpConfig
        type(BicircularSunConfig) :: sun
        logical :: enabled = .true.
        integer :: model = SRP_CANNONBALL_FULL
        real(DP) :: pressure_1au = 4.56e-6_DP ! N/m^2
        real(DP) :: cr0 = 1.2_DP
        real(DP) :: area_mass = 0.02_DP       ! m^2/kg
    end type

    type, public :: CrtbpSrpForceTempPool
        type(CrtbpForceTempPool) :: gravity
        type(DA) :: d(3), er(3), et(3), en(3)
        type(DA) :: distance2, distance, q, transverse2, transverse
    contains
        procedure :: init => srp_pool_init
        procedure :: destroy => srp_pool_destroy
    end type

    public :: srp_acceleration, srp_acceleration_real, srp_acceleration_da
    public :: crtbp_srp_derivatives, crtbp_srp_derivatives_real, crtbp_srp_derivatives_da
    public :: srp_init_da_parameters, validate_crtbp_srp
    interface srp_acceleration
        module procedure srp_acceleration_real, srp_acceleration_da
    end interface
    interface crtbp_srp_derivatives
        module procedure crtbp_srp_derivatives_real, crtbp_srp_derivatives_da
    end interface
contains
    subroutine srp_pool_init(this)
        class(CrtbpSrpForceTempPool), intent(inout) :: this
        integer :: i
        call this%gravity%init()
        do i=1,3
            call this%d(i)%init(); call this%er(i)%init()
            call this%et(i)%init(); call this%en(i)%init()
        end do
        call this%distance2%init(); call this%distance%init(); call this%q%init()
        call this%transverse2%init(); call this%transverse%init()
    end subroutine

    subroutine srp_pool_destroy(this)
        class(CrtbpSrpForceTempPool), intent(inout) :: this
        integer :: i
        call this%gravity%destroy()
        do i=1,3
            call this%d(i)%destroy(); call this%er(i)%destroy()
            call this%et(i)%destroy(); call this%en(i)%destroy()
        end do
        call this%distance2%destroy(); call this%distance%destroy(); call this%q%destroy()
        call this%transverse2%destroy(); call this%transverse%destroy()
    end subroutine

    subroutine validate_crtbp_srp(cfg)
        type(CrtbpSrpConfig), intent(in) :: cfg
        if (cfg%model<SRP_OFF .or. cfg%model>SRP_EMPIRICAL_THREE_COMPONENT) &
            error stop 'CRTBP SRP: unknown model'
        if (.not.all(ieee_is_finite([cfg%pressure_1au,cfg%cr0,cfg%area_mass])) .or. &
            min(cfg%pressure_1au,cfg%cr0,cfg%area_mass)<0.0_DP) &
            error stop 'CRTBP SRP: pressure, Cr and A/m must be finite and nonnegative'
        call validate_bicircular_sun(cfg%sun)
    end subroutine

    logical function srp_is_off(cfg) result(off)
        type(CrtbpSrpConfig), intent(in) :: cfg
        off = .not.cfg%enabled .or. cfg%model==SRP_OFF .or. &
            cfg%pressure_1au==0.0_DP .or. cfg%cr0==0.0_DP .or. cfg%area_mass==0.0_DP
    end function

    ! Caller initializes DACE(order,count(active)) BEFORE creating any DA handles.
    ! Inactive coefficients are exactly zero. Active variables are packed r,t,n.
    ! Scales are half-widths for xi in [-1,1], not Gaussian standard deviations.
    subroutine srp_init_da_parameters(active,scales,eta)
        logical, intent(in) :: active(3)
        real(DP), intent(in) :: scales(3)
        type(DA), intent(inout) :: eta(3)
        integer :: i,j
        if (count(active)<1 .or. dace_max_variables()/=count(active)) &
            error stop 'CRTBP SRP: initialize DACE with exactly the number of active SRP parameters'
        if (.not.all(ieee_is_finite(scales)) .or. any(scales<0.0_DP) .or. &
            any(active .and. scales<=0.0_DP)) error stop 'CRTBP SRP: active parameter scales must be positive'
        j=0
        do i=1,3
            eta(i)=0.0_DP
            if (.not.active(i)) cycle
            j=j+1
            eta(i)=scales(i)*da_var(j)
        end do
    end subroutine

    ! Optional prescribed Sun is a provider boundary: barycenter -> Sun, rotating ND.
    ! Far-field uses its radius/direction too, ignoring spacecraft parallax only.
    subroutine get_sun(tau,cfg,sun,prescribed_sun_nd)
        real(DP), intent(in) :: tau
        type(CrtbpSrpConfig), intent(in) :: cfg
        real(DP), intent(out) :: sun(3)
        real(DP), optional, intent(in) :: prescribed_sun_nd(3)
        call validate_crtbp_srp(cfg)
        if (present(prescribed_sun_nd)) then
            sun=prescribed_sun_nd
            if (.not.all(ieee_is_finite(sun))) error stop 'CRTBP SRP: non-finite prescribed Sun'
        else
            call bicircular_sun_position(tau,cfg%sun,sun)
        end if
    end subroutine

    subroutine check_distance(distance2)
        real(DP), intent(in) :: distance2
        if (.not.ieee_is_finite(distance2) .or. distance2<=min_distance_nd**2) &
            error stop 'CRTBP SRP: Sun-spacecraft distance is singular or non-finite'
    end subroutine

    subroutine srp_acceleration_real(tau,r,eta,cfg,a,a_si,used_fallback,prescribed_sun_nd)
        real(DP), intent(in) :: tau,r(3),eta(3)
        type(CrtbpSrpConfig), intent(in) :: cfg
        real(DP), intent(out) :: a(3)
        real(DP), optional, intent(out) :: a_si(3)
        logical, optional, intent(out) :: used_fallback
        real(DP), optional, intent(in) :: prescribed_sun_nd(3)
        real(DP) :: sun(3),d(3),er(3),et(3),en(3),d2,q,gamma
        logical :: fallback
        a=0.0_DP
        if (present(a_si)) a_si=0.0_DP
        if (present(used_fallback)) used_fallback=.false.
        if (srp_is_off(cfg)) return
        call get_sun(tau,cfg,sun,prescribed_sun_nd)
        if (.not.all(ieee_is_finite(eta))) error stop 'CRTBP SRP: non-finite parameter'
        d=r-sun
        if (cfg%model==SRP_CANNONBALL_FARFIELD) d=-sun
        d2=sum(d*d)
        call check_distance(d2)
        er=d/sqrt(d2)
        q=(cfg%sun%au_m/cfg%sun%du_m)**2/d2
        gamma=cfg%pressure_1au*cfg%cr0*cfg%area_mass*cfg%sun%tu_s**2/cfg%sun%du_m
        a=(1.0_DP+eta(1))*er
        if (cfg%model==SRP_EMPIRICAL_THREE_COMPONENT) then
            et=[-er(2),er(1),0.0_DP] ! k cross er
            fallback=sum(et*et)<min_transverse_squared
            if (fallback) et=[0.0_DP,-er(3),er(2)] ! x-axis cross er
            et=et/sqrt(sum(et*et))
            en=[er(2)*et(3)-er(3)*et(2),er(3)*et(1)-er(1)*et(3),er(1)*et(2)-er(2)*et(1)]
            a=a+eta(2)*et+eta(3)*en
            if (present(used_fallback)) used_fallback=fallback
        end if
        a=gamma*q*a
        if (present(a_si)) a_si=a*(cfg%sun%du_m/cfg%sun%tu_s**2)
    end subroutine

    subroutine srp_acceleration_da(tau,r,eta,cfg,a,pool,a_si,used_fallback,prescribed_sun_nd)
        real(DP), intent(in) :: tau
        type(DA), intent(in) :: r(3),eta(3)
        type(CrtbpSrpConfig), intent(in) :: cfg
        type(DA), intent(inout) :: a(3)
        type(CrtbpSrpForceTempPool), intent(inout) :: pool
        type(DA), optional, intent(inout) :: a_si(3)
        logical, optional, intent(out) :: used_fallback
        real(DP), optional, intent(in) :: prescribed_sun_nd(3)
        real(DP) :: sun(3),gamma
        integer :: i
        logical :: fallback
        do i=1,3
            a(i)=0.0_DP
            if (present(a_si)) a_si(i)=0.0_DP
        end do
        if (present(used_fallback)) used_fallback=.false.
        if (srp_is_off(cfg)) return
        call get_sun(tau,cfg,sun,prescribed_sun_nd)
        do i=1,3
            if (.not.ieee_is_finite(eta(i)%cons())) error stop 'CRTBP SRP: non-finite DA parameter constant'
            if (cfg%model==SRP_CANNONBALL_FARFIELD) then
                pool%d(i)=-sun(i)
            else
                pool%d(i)=r(i)-sun(i)
            end if
        end do
        pool%distance2=pool%d(1)*pool%d(1)+pool%d(2)*pool%d(2)+pool%d(3)*pool%d(3)
        ! Constant part only selects validity/branch, NEVER replaces DA geometry.
        call check_distance(pool%distance2%cons())
        pool%distance=sqrt(pool%distance2)
        pool%q=(cfg%sun%au_m/cfg%sun%du_m)**2/pool%distance2
        gamma=cfg%pressure_1au*cfg%cr0*cfg%area_mass*cfg%sun%tu_s**2/cfg%sun%du_m
        do i=1,3
            pool%er(i)=pool%d(i)/pool%distance
            a(i)=(1.0_DP+eta(1))*pool%er(i)
        end do
        if (cfg%model==SRP_EMPIRICAL_THREE_COMPONENT) then
            pool%et(1)=-pool%er(2); pool%et(2)=pool%er(1); pool%et(3)=0.0_DP
            pool%transverse2=pool%et(1)*pool%et(1)+pool%et(2)*pool%et(2)
            fallback=pool%transverse2%cons()<min_transverse_squared
            if (fallback) then
                pool%et(1)=0.0_DP; pool%et(2)=-pool%er(3); pool%et(3)=pool%er(2)
                pool%transverse2=pool%et(2)*pool%et(2)+pool%et(3)*pool%et(3)
            end if
            pool%transverse=sqrt(pool%transverse2)
            do i=1,3
                pool%et(i)=pool%et(i)/pool%transverse
            end do
            pool%en(1)=pool%er(2)*pool%et(3)-pool%er(3)*pool%et(2)
            pool%en(2)=pool%er(3)*pool%et(1)-pool%er(1)*pool%et(3)
            pool%en(3)=pool%er(1)*pool%et(2)-pool%er(2)*pool%et(1)
            do i=1,3
                a(i)=a(i)+eta(2)*pool%et(i)+eta(3)*pool%en(i)
            end do
            if (present(used_fallback)) used_fallback=fallback
        end if
        do i=1,3
            a(i)=gamma*pool%q*a(i)
            if (present(a_si)) a_si(i)=a(i)*(cfg%sun%du_m/cfg%sun%tu_s**2)
        end do
    end subroutine

    subroutine crtbp_srp_derivatives_real(x,dxdt,tau,cfg,eta)
        real(DP), intent(in) :: x(6),tau,eta(3)
        real(DP), intent(out) :: dxdt(6)
        type(CrtbpSrpConfig), intent(in) :: cfg
        real(DP) :: a(3)
        call crtbp_derivatives_real(x,dxdt,tau)
        call srp_acceleration_real(tau,x(1:3),eta,cfg,a)
        ! This is an external force vector already expressed in rotating axes.
        ! The base RHS already contains ALL centrifugal and Coriolis terms.
        dxdt(4:6)=dxdt(4:6)+a
    end subroutine

    subroutine crtbp_srp_derivatives_da(x,dxdt,tau,cfg,eta,pool)
        type(DA), intent(in) :: x(6),eta(3)
        type(DA), intent(inout) :: dxdt(6)
        real(DP), intent(in) :: tau
        type(CrtbpSrpConfig), intent(in) :: cfg
        type(CrtbpSrpForceTempPool), intent(inout) :: pool
        ! Separate output from pool to avoid aliasing an intent(inout) component
        ! with its owning object. These three handles are explicitly destroyed.
        type(DA) :: a(3)
        integer :: i
        call crtbp_derivatives_da(x,dxdt,tau,pool%gravity)
        if (srp_is_off(cfg)) return
        call srp_acceleration_da(tau,x(1:3),eta,cfg,a,pool)
        do i=1,3
            dxdt(i+3)=dxdt(i+3)+a(i)
            call a(i)%destroy()
        end do
    end subroutine
end module pod_crtbp_srp_module
