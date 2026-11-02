!> Q1-Q4 光压实验的物理参数和 DA 参数域。
!! Physical order is ALWAYS [beta, chi, u, v]; active DA indices are packed.
!! beta: face-on equivalent Cr*A/m [m2/kg], not a solar-gravity ratio.
!! chi: reflection composition; u,v: dimensionless solar-tangent slopes.
!! Bounds describe a deterministic domain, NOT a probability distribution.
module pod_srp_study_parameters_module
    use pod_global, only: DP
    use pod_dace_classes, only: DA, da_var, dace_max_variables, assignment(=), operator(+), operator(*)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private
    integer, parameter, public :: SRP_BETA=1, SRP_CHI=2, SRP_U=3, SRP_V=4
    character(len=4), parameter, public :: SRP_PARAMETER_NAMES(4)=[character(len=4)::'beta','chi','u','v']
    real(DP), parameter, public :: SRP_FACE_ON_COEFFICIENT=1.5_DP
    type, public :: SrpStudyDomain
        real(DP) :: lower(4)=[0.0_DP,0.6_DP,0.0_DP,0.0_DP]
        real(DP) :: upper(4)=[0.1_DP,0.6_DP,0.0_DP,0.0_DP]
    contains
        procedure :: center => domain_center
        procedure :: halfwidth => domain_halfwidth
        procedure :: n_variables => domain_n_variables
    end type
    public :: srp_study_domain, srp_study_init_da, validate_srp_study_domain
    public :: srp_study_optics, srp_face_on_beta
contains
    !> Reflection family with constant C_eq=1.5; chi=.6 gives rho_s=rho_d=.3.
    !! No extra Cr multiplier belongs in the flat-plate force.
    subroutine srp_study_optics(chi,rho_s,rho_d,absorptivity)
        real(DP), intent(in) :: chi
        real(DP), intent(out) :: rho_s,rho_d,absorptivity
        if (.not.ieee_is_finite(chi) .or. chi<0.0_DP .or. chi>1.0_DP) &
            error stop 'SRP study: chi must be in [0,1]'
        rho_s=0.5_DP*chi
        rho_d=0.75_DP*(1.0_DP-chi)
        absorptivity=0.25_DP+0.25_DP*chi
    end subroutine

    !> Reference geometry only: 80/160 m2 and 1200 kg give beta=.1/.2.
    !! In the study beta is the independent effective-strength coordinate.
    !! A varying beta is NOT a material uncertainty of an otherwise fixed plate.
    function srp_face_on_beta(area_m2,mass_kg) result(beta)
        real(DP), intent(in) :: area_m2,mass_kg
        real(DP) :: beta
        if (.not.all(ieee_is_finite([area_m2,mass_kg])) .or. area_m2<0.0_DP .or. mass_kg<=0.0_DP) &
            error stop 'SRP study: area must be nonnegative and mass positive'
        beta=SRP_FACE_ON_COEFFICIENT*area_m2/mass_kg
    end function

    !> Construct Q1 (beta), Q2/Q3 (beta,u,v), or Q4 (beta,chi,u,v) domains.
    !! theta_max bounds the CORNERS of the uv box, not a uniformly sampled cone.
    !! Fixed parameters have equal lower/upper bounds, e.g. chi=.6 in Q2/Q3.
    function srp_study_domain(beta_lower,beta_upper,tilt_max_deg,chi_lower,chi_upper) result(domain)
        real(DP), intent(in) :: beta_lower,beta_upper
        real(DP), optional, intent(in) :: tilt_max_deg,chi_lower,chi_upper
        type(SrpStudyDomain) :: domain
        real(DP) :: q
        domain%lower(SRP_BETA)=beta_lower; domain%upper(SRP_BETA)=beta_upper
        if (present(chi_lower)) domain%lower(SRP_CHI)=chi_lower
        if (present(chi_upper)) domain%upper(SRP_CHI)=chi_upper
        if (present(tilt_max_deg)) then
            if (.not.ieee_is_finite(tilt_max_deg) .or. tilt_max_deg<0.0_DP .or. tilt_max_deg>=90.0_DP) &
                error stop 'SRP study: maximum tilt must be in [0,90) degrees'
            q=tan(tilt_max_deg*acos(-1.0_DP)/180.0_DP)/sqrt(2.0_DP)
            domain%lower(SRP_U:SRP_V)=-q
            domain%upper(SRP_U:SRP_V)=q
        end if
        call validate_srp_study_domain(domain)
    end function

    subroutine validate_srp_study_domain(domain)
        class(SrpStudyDomain), intent(in) :: domain
        if (.not.all(ieee_is_finite(domain%lower)) .or. .not.all(ieee_is_finite(domain%upper))) &
            error stop 'SRP study: parameter bounds must be finite'
        if (any(domain%upper<domain%lower)) error stop 'SRP study: reversed parameter bounds'
        if (domain%lower(SRP_BETA)<0.0_DP) error stop 'SRP study: beta lower bound must be nonnegative'
        if (domain%lower(SRP_CHI)<0.0_DP .or. domain%upper(SRP_CHI)>1.0_DP) &
            error stop 'SRP study: chi bounds must lie in [0,1]'
    end subroutine

    function domain_center(this) result(center)
        class(SrpStudyDomain), intent(in) :: this
        real(DP) :: center(4)
        call validate_srp_study_domain(this)
        center=0.5_DP*this%lower+0.5_DP*this%upper
    end function

    function domain_halfwidth(this) result(width)
        class(SrpStudyDomain), intent(in) :: this
        real(DP) :: width(4)
        call validate_srp_study_domain(this)
        width=0.5_DP*this%upper-0.5_DP*this%lower
    end function

    integer function domain_n_variables(this) result(n)
        class(SrpStudyDomain), intent(in) :: this
        call validate_srp_study_domain(this)
        n=count(this%upper>this%lower)
    end function

    !> Caller first initializes DACE(order,domain%n_variables()), with no live handles.
    !! 六维初始状态由调用者设为 DA 常数；此函数只为光压参数分配自变量。
    !! index_map(physical_index)=DA_index, or 0 for a fixed parameter. This lets
    !! future ADS logs name a split correctly: in 3D the order is beta,u,v.
    !! After splitting, affine-substitute BOTH the state and parameters using
    !! the same local variables; do not reinitialize the original full domain.
    subroutine srp_study_init_da(domain,parameters,index_map)
        type(SrpStudyDomain), intent(in) :: domain
        type(DA), intent(inout) :: parameters(4)
        integer, optional, intent(out) :: index_map(4)
        real(DP) :: c(4),h(4)
        integer :: i,j,n
        n=domain%n_variables()
        if (n<1 .or. dace_max_variables()/=n) &
            error stop 'SRP study: DACE variable count must match active parameters'
        c=domain%center(); h=domain%halfwidth()
        if (present(index_map)) index_map=0
        j=0
        do i=1,4
            parameters(i)=c(i)
            if (h(i)==0.0_DP) cycle
            j=j+1
            parameters(i)=c(i)+h(i)*da_var(j)
            if (present(index_map)) index_map(i)=j
        end do
    end subroutine
end module pod_srp_study_parameters_module
