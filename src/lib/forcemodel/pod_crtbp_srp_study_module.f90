!> CRTBP 会合系中的绝对 beta 炮弹球/平板模型（Q1-Q4）。
!! All positions, time and returned accelerations are nondimensional.
!! Parameters [beta,chi,u,v] or [beta,chi,theta,phi] are constant over an arc.
!! beta is in m2/kg; theta and phi are in radians relative to the solar tangent basis.
!! The Sun supplies geometry only. Full sunlight, no solar gravity or torque.
!! Existing relative-eta models remain in pod_crtbp_srp_module unchanged.
module pod_crtbp_srp_study_module
    use pod_global, only: DP
    use pod_dace_classes
    use pod_crtbp_module, only: CrtbpForceTempPool, crtbp_derivatives_real, crtbp_derivatives_da
    use pod_bicircular_sun_module
    use pod_srp_study_parameters_module, only: SRP_BETA,SRP_CHI,SRP_U,SRP_V, &
        SRP_FACE_ON_COEFFICIENT,srp_study_optics,SRP_THETA,SRP_PHI,SRP_ATTITUDE_UV,SRP_ATTITUDE_THETA_PHI
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private
    integer, parameter, public :: STUDY_CANNONBALL=1, STUDY_FLAT_PLATE=2
    real(DP), parameter :: min_distance2=1e-24_DP, min_tangent2=1e-20_DP

    type, public :: CrtbpSrpStudyConfig
        ! This study uses the orbit catalog's TU. Old force defaults are preserved.
        type(BicircularSunConfig) :: sun=BicircularSunConfig(tu_s=375237.080249162_DP)
        integer :: model=STUDY_FLAT_PLATE
        logical :: enabled=.true.
        real(DP) :: pressure_1au=4.56e-6_DP ! N/m2, always scaled by (AU/distance)^2
        ! Fixed axis in SYNODIC coordinates, chosen once for the entire arc/domain.
        ! e1=normalize(axis cross s), e2=s cross e1. No implicit axis switching.
        real(DP) :: tangent_reference_axis(3)=[0.0_DP,0.0_DP,1.0_DP]
        real(DP) :: tangent_rotation_rad=0.0_DP ! constant rotation of e1,e2 for coordinate tests
        integer :: attitude_coordinates=SRP_ATTITUDE_UV ! preserve existing p(3:4)=u,v callers
    end type

    type, public :: SrpStudyTempPool
        type(CrtbpForceTempPool) :: gravity
        type(DA) :: d(3),s(3),b1(3),b2(3),e1(3),e2(3),normal(3),acceleration(3)
        type(DA) :: distance2,distance,q,tangent2,tangent,cos_incidence,rho_s,rho_d,reflection,factor
        type(DA) :: sin_theta,cos_phi,sin_phi ! persistent temporaries for angle-mode RHS
    contains
        procedure :: init => study_pool_init
        procedure :: destroy => study_pool_destroy
    end type
    public :: srp_study_acceleration,srp_study_acceleration_real,srp_study_acceleration_da
    public :: crtbp_srp_study_derivatives,crtbp_srp_study_derivatives_real,crtbp_srp_study_derivatives_da
    interface srp_study_acceleration
        module procedure srp_study_acceleration_real,srp_study_acceleration_da
    end interface
    interface crtbp_srp_study_derivatives
        module procedure crtbp_srp_study_derivatives_real,crtbp_srp_study_derivatives_da
    end interface
contains
    subroutine study_pool_init(this)
        class(SrpStudyTempPool),intent(inout)::this
        integer :: i
        call this%gravity%init()
        do i=1,3
            call this%d(i)%init(); call this%s(i)%init()
            call this%b1(i)%init(); call this%b2(i)%init()
            call this%e1(i)%init(); call this%e2(i)%init()
            call this%normal(i)%init(); call this%acceleration(i)%init()
        end do
        call this%distance2%init(); call this%distance%init(); call this%q%init()
        call this%tangent2%init(); call this%tangent%init(); call this%cos_incidence%init()
        call this%sin_theta%init(); call this%cos_phi%init(); call this%sin_phi%init()
        call this%rho_s%init(); call this%rho_d%init(); call this%reflection%init(); call this%factor%init()
    end subroutine

    subroutine study_pool_destroy(this)
        class(SrpStudyTempPool),intent(inout)::this
        integer :: i
        call this%gravity%destroy()
        do i=1,3
            call this%d(i)%destroy(); call this%s(i)%destroy()
            call this%b1(i)%destroy(); call this%b2(i)%destroy()
            call this%e1(i)%destroy(); call this%e2(i)%destroy()
            call this%normal(i)%destroy(); call this%acceleration(i)%destroy()
        end do
        call this%distance2%destroy(); call this%distance%destroy(); call this%q%destroy()
        call this%tangent2%destroy(); call this%tangent%destroy(); call this%cos_incidence%destroy()
        call this%sin_theta%destroy(); call this%cos_phi%destroy(); call this%sin_phi%destroy()
        call this%rho_s%destroy(); call this%rho_d%destroy(); call this%reflection%destroy(); call this%factor%destroy()
    end subroutine

    subroutine validate_study_config(cfg)
        type(CrtbpSrpStudyConfig),intent(in)::cfg
        call validate_bicircular_sun(cfg%sun)
        if (cfg%attitude_coordinates/=SRP_ATTITUDE_UV .and. &
            cfg%attitude_coordinates/=SRP_ATTITUDE_THETA_PHI) error stop 'SRP study: unknown attitude coordinates'
        if (cfg%model/=STUDY_CANNONBALL .and. cfg%model/=STUDY_FLAT_PLATE) error stop 'SRP study: unknown model'
        if (.not.ieee_is_finite(cfg%pressure_1au) .or. cfg%pressure_1au<0.0_DP) &
            error stop 'SRP study: pressure must be finite and nonnegative'
        if (.not.all(ieee_is_finite(cfg%tangent_reference_axis)) .or. &
            .not.ieee_is_finite(cfg%tangent_rotation_rad)) error stop 'SRP study: non-finite tangent frame'
        if (norm2(cfg%tangent_reference_axis)<=1e-12_DP) error stop 'SRP study: zero tangent reference axis'
    end subroutine

    subroutine validate_parameter_center(p,cfg)
        real(DP),intent(in)::p(4)
        type(CrtbpSrpStudyConfig),intent(in)::cfg
        if (.not.all(ieee_is_finite(p))) error stop 'SRP study: non-finite parameters'
        if (p(SRP_BETA)<0.0_DP) error stop 'SRP study: beta must be nonnegative'
        if (p(SRP_CHI)<0.0_DP .or. p(SRP_CHI)>1.0_DP) error stop 'SRP study: chi must be in [0,1]'
        if (cfg%model==STUDY_FLAT_PLATE .and. cfg%attitude_coordinates==SRP_ATTITUDE_THETA_PHI) then
            if (p(SRP_THETA)<0.0_DP .or. p(SRP_THETA)>=acos(-1.0_DP)/2) &
                error stop 'SRP study: theta must lie in [0,pi/2) radians'
        end if
    end subroutine

    subroutine check_geometry(distance2)
        real(DP),intent(in)::distance2
        if (.not.ieee_is_finite(distance2) .or. distance2<=min_distance2) &
            error stop 'SRP study: singular Sun-spacecraft distance'
    end subroutine

    subroutine check_tangent(tangent2)
        real(DP),intent(in)::tangent2
        if (.not.ieee_is_finite(tangent2) .or. tangent2<=min_tangent2) &
            error stop 'SRP study: tangent axis parallel to Sun direction; choose a fixed regular axis'
    end subroutine

    pure function cross(a,b) result(c)
        real(DP),intent(in)::a(3),b(3)
        real(DP)::c(3)
        c=[a(2)*b(3)-a(3)*b(2),a(3)*b(1)-a(1)*b(3),a(1)*b(2)-a(2)*b(1)]
    end function

    !> s points SPACECRAFT -> SUN, so acceleration carries a negative sign.
    !! n=(s+u*e1+v*e2)/sqrt(1+u^2+v^2); finite slopes keep cos_incidence>0.
    !! In angle mode n=cos(theta)*s+sin(theta)*(cos(phi)*e1+sin(phi)*e2).
    !! Attitude parameters are fixed relative to Sun tracking, not inertial normals.
    subroutine srp_study_acceleration_real(tau,r,p,cfg,a,a_si)
        real(DP),intent(in)::tau,r(3),p(4)
        type(CrtbpSrpStudyConfig),intent(in)::cfg
        real(DP),intent(out)::a(3)
        real(DP),optional,intent(out)::a_si(3)
        real(DP)::sun(3),d(3),s(3),axis(3),b1(3),b2(3),e1(3),e2(3),normal(3)
        real(DP)::d2,q,pressure_nd,cinc,rho_s,rho_d,absorptivity,angle
        call validate_study_config(cfg)
        a=0.0_DP
        if (present(a_si)) a_si=0.0_DP
        if (.not.cfg%enabled .or. cfg%pressure_1au==0.0_DP) return
        call validate_parameter_center(p,cfg)
        call bicircular_sun_position(tau,cfg%sun,sun)
        d=sun-r; d2=sum(d*d)
        call check_geometry(d2)
        s=d/sqrt(d2)
        q=(cfg%sun%au_m/cfg%sun%du_m)**2/d2
        pressure_nd=cfg%pressure_1au*cfg%sun%tu_s**2/cfg%sun%du_m
        if (cfg%model==STUDY_CANNONBALL) then
            a=-pressure_nd*q*p(SRP_BETA)*s
        else
            axis=cfg%tangent_reference_axis/norm2(cfg%tangent_reference_axis)
            b1=cross(axis,s)
            call check_tangent(sum(b1*b1))
            b1=b1/norm2(b1); b2=cross(s,b1)
            angle=cfg%tangent_rotation_rad
            e1=cos(angle)*b1+sin(angle)*b2; e2=-sin(angle)*b1+cos(angle)*b2
            if (cfg%attitude_coordinates==SRP_ATTITUDE_THETA_PHI) then
                cinc=cos(p(SRP_THETA))
                normal=cinc*s+sin(p(SRP_THETA))*(cos(p(SRP_PHI))*e1+sin(p(SRP_PHI))*e2)
            else
                cinc=1.0_DP/sqrt(1.0_DP+p(SRP_U)**2+p(SRP_V)**2)
                normal=cinc*(s+p(SRP_U)*e1+p(SRP_V)*e2)
            end if
            call srp_study_optics(p(SRP_CHI),rho_s,rho_d,absorptivity)
            ! beta/C_eq is the effective A/m. Do NOT multiply by another Cr.
            a=-pressure_nd*q*(p(SRP_BETA)/SRP_FACE_ON_COEFFICIENT)*cinc * &
                ((1.0_DP-rho_s)*s+2.0_DP*(rho_s*cinc+rho_d/3.0_DP)*normal)
        end if
        if (present(a_si)) a_si=a*(cfg%sun%du_m/cfg%sun%tu_s**2)
    end subroutine

    ! Private kernel writes only into its own pool: no output/pool aliasing.
    ! Constant terms are inspected for validity only; ALL geometry stays DA.
    subroutine study_da_kernel(tau,r,p,cfg,pool)
        real(DP),intent(in)::tau
        type(DA),intent(in)::r(3),p(4)
        type(CrtbpSrpStudyConfig),intent(in)::cfg
        type(SrpStudyTempPool),intent(inout)::pool
        real(DP)::sun(3),axis(3),values(4),pressure_nd,c,s
        integer::i
        call validate_study_config(cfg)
        do i=1,3
            pool%acceleration(i)=0.0_DP
        end do
        if (.not.cfg%enabled .or. cfg%pressure_1au==0.0_DP) return
        do i=1,4
            values(i)=p(i)%cons()
        end do
        call validate_parameter_center(values,cfg)
        call bicircular_sun_position(tau,cfg%sun,sun)
        do i=1,3
            pool%d(i)=sun(i)-r(i)
        end do
        pool%distance2=pool%d(1)*pool%d(1)+pool%d(2)*pool%d(2)+pool%d(3)*pool%d(3)
        call check_geometry(pool%distance2%cons())
        pool%distance=sqrt(pool%distance2)
        pool%q=(cfg%sun%au_m/cfg%sun%du_m)**2/pool%distance2
        pressure_nd=cfg%pressure_1au*cfg%sun%tu_s**2/cfg%sun%du_m
        do i=1,3
            pool%s(i)=pool%d(i)/pool%distance
        end do
        ! DO NOT short-circuit on beta%cons()==0: its DA derivative can be nonzero.
        if (cfg%model==STUDY_CANNONBALL) then
            pool%factor=-pressure_nd*pool%q*p(SRP_BETA)
            do i=1,3
                pool%acceleration(i)=pool%factor*pool%s(i)
            end do
            return
        end if
        axis=cfg%tangent_reference_axis/norm2(cfg%tangent_reference_axis)
        pool%b1(1)=axis(2)*pool%s(3)-axis(3)*pool%s(2)
        pool%b1(2)=axis(3)*pool%s(1)-axis(1)*pool%s(3)
        pool%b1(3)=axis(1)*pool%s(2)-axis(2)*pool%s(1)
        pool%tangent2=pool%b1(1)*pool%b1(1)+pool%b1(2)*pool%b1(2)+pool%b1(3)*pool%b1(3)
        call check_tangent(pool%tangent2%cons())
        pool%tangent=sqrt(pool%tangent2)
        do i=1,3
            pool%b1(i)=pool%b1(i)/pool%tangent
        end do
        pool%b2(1)=pool%s(2)*pool%b1(3)-pool%s(3)*pool%b1(2)
        pool%b2(2)=pool%s(3)*pool%b1(1)-pool%s(1)*pool%b1(3)
        pool%b2(3)=pool%s(1)*pool%b1(2)-pool%s(2)*pool%b1(1)
        c=cos(cfg%tangent_rotation_rad); s=sin(cfg%tangent_rotation_rad)
        do i=1,3
            pool%e1(i)=c*pool%b1(i)+s*pool%b2(i)
            pool%e2(i)=-s*pool%b1(i)+c*pool%b2(i)
        end do
        if (cfg%attitude_coordinates==SRP_ATTITUDE_THETA_PHI) then
            ! Forward trig map stays smooth at theta=0, including mixed terms.
            ! Do not wrap phi or replace either angle with its DA constant term.
            pool%cos_incidence=cos(p(SRP_THETA))
            pool%sin_theta=sin(p(SRP_THETA))
            pool%cos_phi=cos(p(SRP_PHI)); pool%sin_phi=sin(p(SRP_PHI))
            do i=1,3
                pool%normal(i)=pool%cos_incidence*pool%s(i)+pool%sin_theta* &
                    (pool%cos_phi*pool%e1(i)+pool%sin_phi*pool%e2(i))
            end do
        else
            pool%cos_incidence=1.0_DP/sqrt(1.0_DP+p(SRP_U)*p(SRP_U)+p(SRP_V)*p(SRP_V))
            do i=1,3
                pool%normal(i)=pool%cos_incidence*(pool%s(i)+p(SRP_U)*pool%e1(i)+p(SRP_V)*pool%e2(i))
            end do
        end if
        ! Same chi family as srp_study_optics, now retaining mixed DA coefficients.
        pool%rho_s=0.5_DP*p(SRP_CHI)
        pool%rho_d=0.75_DP*(1.0_DP-p(SRP_CHI))
        pool%reflection=2.0_DP*(pool%rho_s*pool%cos_incidence+pool%rho_d/3.0_DP)
        pool%factor=(-pressure_nd/SRP_FACE_ON_COEFFICIENT)*pool%q*p(SRP_BETA)*pool%cos_incidence
        do i=1,3
            pool%acceleration(i)=pool%factor*((1.0_DP-pool%rho_s)*pool%s(i)+pool%reflection*pool%normal(i))
        end do
    end subroutine

    !> Caller owns/init/destroys the pool. Inputs, outputs and pool must not alias.
    !! a_si is optional and also DA; no truncation to nominal acceleration.
    subroutine srp_study_acceleration_da(tau,r,p,cfg,a,pool,a_si)
        real(DP),intent(in)::tau
        type(DA),intent(in)::r(3),p(4)
        type(CrtbpSrpStudyConfig),intent(in)::cfg
        type(DA),intent(inout)::a(3)
        type(SrpStudyTempPool),intent(inout)::pool
        type(DA),optional,intent(inout)::a_si(3)
        integer::i
        call study_da_kernel(tau,r,p,cfg,pool)
        do i=1,3
            a(i)=pool%acceleration(i)
            if (present(a_si)) a_si(i)=a(i)*(cfg%sun%du_m/cfg%sun%tu_s**2)
        end do
    end subroutine

    !> Reuse current crtbp_mu (set_crtbp_mu) and unchanged gravity/Coriolis terms.
    !! SRP is already in synodic components and receives no extra inertial terms.
    subroutine crtbp_srp_study_derivatives_real(x,dxdt,tau,cfg,p)
        real(DP),intent(in)::x(6),tau,p(4)
        real(DP),intent(out)::dxdt(6)
        type(CrtbpSrpStudyConfig),intent(in)::cfg
        real(DP)::a(3)
        call crtbp_derivatives_real(x,dxdt,tau)
        call srp_study_acceleration_real(tau,x(1:3),p,cfg,a)
        dxdt(4:6)=dxdt(4:6)+a
    end subroutine

    subroutine crtbp_srp_study_derivatives_da(x,dxdt,tau,cfg,p,pool)
        type(DA),intent(in)::x(6),p(4)
        type(DA),intent(inout)::dxdt(6)
        real(DP),intent(in)::tau
        type(CrtbpSrpStudyConfig),intent(in)::cfg
        type(SrpStudyTempPool),intent(inout)::pool
        integer::i
        call crtbp_derivatives_da(x,dxdt,tau,pool%gravity)
        call study_da_kernel(tau,x(1:3),p,cfg,pool)
        do i=1,3
            dxdt(i+3)=dxdt(i+3)+pool%acceleration(i)
        end do
    end subroutine
end module pod_crtbp_srp_study_module
