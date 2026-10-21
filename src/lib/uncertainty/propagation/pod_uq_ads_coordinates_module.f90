!> Coordinate maps between physical uncertainty and the ADS unit box.
module pod_uq_ads_coordinates_module
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use pod_global, only: DP
    use pod_basicmath_module, only: dpotrf
    implicit none
    private

    integer, parameter, public :: ADS_COORD_COMPONENT = 1
    integer, parameter, public :: ADS_COORD_WHITENED = 2
    integer, parameter, public :: ADS_PARAM_AUTO=0, ADS_PARAM_NONE=1, ADS_PARAM_SRP_ETA=2, ADS_PARAM_SRP_MODES=3

    type, public :: ads_physical_parameters_type
        real(DP) :: eta=0.0_DP, s=0.0_DP, u=0.0_DP
    end type

    type, public :: ads_coordinate_map_type
        integer :: mode = ADS_COORD_COMPONENT
        integer :: n_variables = 6
        integer :: parameterization=ADS_PARAM_NONE
        real(DP) :: domain_sigma = 3.0_DP
        real(DP) :: srp_sigma = 0.0_DP
        real(DP) :: mean6(6) = 0.0_DP
        real(DP) :: basis6(6,6) = 0.0_DP
    end type ads_coordinate_map_type

    public :: ads_build_coordinate_map
    public :: ads_physical_to_unit
    public :: ads_unit_to_physical
    public :: ads_patch_mode_parameters

    interface ads_physical_to_unit
        module procedure physical_to_unit_eta, physical_to_unit_parameters
    end interface
    interface ads_unit_to_physical
        module procedure unit_to_physical_eta, unit_to_physical_parameters
    end interface

contains

    subroutine ads_build_coordinate_map(mean6,covariance,mode,domain_sigma, &
                                        srp_sigma,map,status,message,parameterization)
        real(DP), intent(in) :: mean6(6),covariance(6,6)
        integer, intent(in) :: mode
        real(DP), intent(in) :: domain_sigma,srp_sigma
        type(ads_coordinate_map_type), intent(out) :: map
        integer, intent(out) :: status
        character(len=*), intent(out) :: message
        integer,intent(in),optional :: parameterization
        real(DP) :: symmetric(6,6),scale,symmetry_tolerance
        integer :: i,j,info,resolved

        map=ads_coordinate_map_type()
        status=0
        message=''
        if(.not.all(ieee_is_finite(mean6)) .or. &
           .not.all(ieee_is_finite(covariance))) then
            call fail(-1,'ADS mean and covariance must be finite')
            return
        end if
        if(.not.ieee_is_finite(domain_sigma) .or. domain_sigma<=0.0_DP) then
            call fail(-2,'ADS domain_sigma must be positive')
            return
        end if
        if(.not.ieee_is_finite(srp_sigma) .or. srp_sigma<0.0_DP) then
            call fail(-3,'ADS srp_sigma must be non-negative')
            return
        end if
        if(mode/=ADS_COORD_COMPONENT .and. mode/=ADS_COORD_WHITENED) then
            call fail(-4,'unknown ADS coordinate mode')
            return
        end if
        resolved=ADS_PARAM_AUTO
        if(present(parameterization)) resolved=parameterization
        if(resolved==ADS_PARAM_AUTO) resolved=merge(ADS_PARAM_SRP_ETA,ADS_PARAM_NONE,srp_sigma>0.0_DP)
        select case(resolved)
        case(ADS_PARAM_NONE,ADS_PARAM_SRP_MODES)
            if(srp_sigma/=0.0_DP) then
                call fail(-8,'NONE/MODES cannot include eta uncertainty'); return
            end if
        case(ADS_PARAM_SRP_ETA)
            if(srp_sigma<=0.0_DP) then
                call fail(-8,'ETA requires positive srp_sigma'); return
            end if
        case default
            call fail(-8,'unknown ADS parameterization'); return
        end select

        scale=max(1.0_DP,maxval(abs(covariance)))
        symmetry_tolerance=1000.0_DP*epsilon(1.0_DP)*scale
        if(maxval(abs(covariance-transpose(covariance)))>symmetry_tolerance) then
            call fail(-5,'ADS covariance must be symmetric')
            return
        end if

        map%mode=mode
        map%n_variables=merge(7,6,srp_sigma>0.0_DP)
        map%parameterization=resolved
        if(resolved==ADS_PARAM_SRP_MODES) map%n_variables=8
        map%domain_sigma=domain_sigma
        map%srp_sigma=srp_sigma
        map%mean6=mean6
        map%basis6=0.0_DP
        symmetric=0.5_DP*(covariance+transpose(covariance))

        select case(mode)
        case(ADS_COORD_COMPONENT)
            do i=1,6
                if(symmetric(i,i)<-symmetry_tolerance) then
                    call fail(-6,'ADS covariance has a negative diagonal')
                    return
                end if
                map%basis6(i,i)=domain_sigma*sqrt(max(symmetric(i,i),0.0_DP))
            end do
        case(ADS_COORD_WHITENED)
            map%basis6=symmetric
            call dpotrf('L',6,map%basis6,6,info)
            if(info/=0) then
                map%basis6=0.0_DP
                call fail(-7,'whitened ADS coordinates require positive-definite covariance')
                return
            end if
            do j=2,6
                do i=1,j-1
                    map%basis6(i,j)=0.0_DP
                end do
            end do
            map%basis6=domain_sigma*map%basis6
        end select

    contains
        subroutine fail(code,text)
            integer,intent(in) :: code
            character(len=*),intent(in) :: text
            status=code
            message=text
        end subroutine fail
    end subroutine ads_build_coordinate_map

    subroutine physical_to_unit_eta(map,deviation,eta_srp,unit_point,status)
        type(ads_coordinate_map_type), intent(in) :: map
        real(DP), intent(in) :: deviation(6),eta_srp
        real(DP), allocatable, intent(out) :: unit_point(:)
        integer, intent(out) :: status
        real(DP) :: fixed_tolerance,srp_scale
        integer :: i,j

        status=0
        allocate(unit_point(map%n_variables))
        unit_point=0.0_DP
        if(map%parameterization==ADS_PARAM_SRP_MODES) then
            status=-7; return
        end if
        if(.not.all(ieee_is_finite(deviation)) .or. &
           .not.ieee_is_finite(eta_srp)) then
            status=-1
            return
        end if
        select case(map%mode)
        case(ADS_COORD_COMPONENT)
            do i=1,6
                if(map%basis6(i,i)>0.0_DP) then
                    unit_point(i)=deviation(i)/map%basis6(i,i)
                else
                    fixed_tolerance=1000.0_DP*epsilon(1.0_DP)* &
                        max(1.0_DP,abs(map%mean6(i)))
                    if(abs(deviation(i))>fixed_tolerance) then
                        status=-2
                        return
                    end if
                end if
            end do
        case(ADS_COORD_WHITENED)
            do i=1,6
                if(map%basis6(i,i)<=0.0_DP) then
                    status=-3
                    return
                end if
                unit_point(i)=deviation(i)
                do j=1,i-1
                    unit_point(i)=unit_point(i)-map%basis6(i,j)*unit_point(j)
                end do
                unit_point(i)=unit_point(i)/map%basis6(i,i)
            end do
        case default
            status=-4
            return
        end select
        if(map%n_variables==7) then
            srp_scale=map%domain_sigma*map%srp_sigma
            if(srp_scale<=0.0_DP) then
                status=-5
                return
            end if
            unit_point(7)=eta_srp/srp_scale
        else if(abs(eta_srp)>1000.0_DP*epsilon(1.0_DP)) then
            status=-6
        end if
    end subroutine physical_to_unit_eta

    subroutine unit_to_physical_eta(map,unit_point,deviation,eta_srp,status)
        type(ads_coordinate_map_type), intent(in) :: map
        real(DP), intent(in) :: unit_point(:)
        real(DP), intent(out) :: deviation(6),eta_srp
        integer, intent(out) :: status

        status=0
        deviation=0.0_DP
        eta_srp=0.0_DP
        if(map%parameterization==ADS_PARAM_SRP_MODES) then
            status=-3; return
        end if
        if(size(unit_point)/=map%n_variables) then
            status=-1
            return
        end if
        if(.not.all(ieee_is_finite(unit_point))) then
            status=-2
            return
        end if
        deviation=matmul(map%basis6,unit_point(1:6))
        if(map%n_variables==7) then
            eta_srp=map%domain_sigma*map%srp_sigma*unit_point(7)
        end if
    end subroutine unit_to_physical_eta

    subroutine physical_to_unit_parameters(map,deviation,parameters,unit_point,status)
        type(ads_coordinate_map_type),intent(in) :: map
        real(DP),intent(in) :: deviation(6)
        type(ads_physical_parameters_type),intent(in) :: parameters
        real(DP),allocatable,intent(out) :: unit_point(:)
        integer,intent(out) :: status
        type(ads_coordinate_map_type) :: state_map
        real(DP),allocatable :: state_point(:)
        real(DP) :: tolerance,s,u
        if(map%parameterization/=ADS_PARAM_SRP_MODES) then
            call physical_to_unit_eta(map,deviation,parameters%eta,unit_point,status); return
        end if
        allocate(unit_point(8)); unit_point=0; status=-1
        tolerance=1000*epsilon(1.0_DP)
        if(.not.all(ieee_is_finite([parameters%s,parameters%u,parameters%eta]))) return
        if(abs(parameters%eta)>tolerance) return
        if(parameters%s< -tolerance .or. parameters%s>1+tolerance .or. &
           parameters%u< -tolerance .or. parameters%u>1+tolerance) return
        state_map=map; state_map%n_variables=6; state_map%parameterization=ADS_PARAM_NONE
        call physical_to_unit_eta(state_map,deviation,0.0_DP,state_point,status)
        if(status/=0) return
        if(any(abs(state_point)>1+tolerance)) then
            status=-2; return
        end if
        s=max(0.0_DP,min(1.0_DP,parameters%s)); u=max(0.0_DP,min(1.0_DP,parameters%u))
        unit_point(1:6)=max(-1.0_DP,min(1.0_DP,state_point))
        unit_point(7:8)=2*[s,u]-1
    end subroutine

    subroutine unit_to_physical_parameters(map,unit_point,deviation,parameters,status)
        type(ads_coordinate_map_type),intent(in) :: map
        real(DP),intent(in) :: unit_point(:)
        real(DP),intent(out) :: deviation(6)
        type(ads_physical_parameters_type),intent(out) :: parameters
        integer,intent(out) :: status
        real(DP) :: q(8),tolerance
        parameters=ads_physical_parameters_type(); deviation=0
        if(map%parameterization/=ADS_PARAM_SRP_MODES) then
            call unit_to_physical_eta(map,unit_point,deviation,parameters%eta,status); return
        end if
        status=-1
        if(size(unit_point)/=8) return
        if(.not.all(ieee_is_finite(unit_point))) return
        tolerance=1000*epsilon(1.0_DP)
        if(any(abs(unit_point)>1+tolerance)) return
        q=max(-1.0_DP,min(1.0_DP,unit_point))
        deviation=matmul(map%basis6,q(1:6))
        parameters%s=0.5_DP*(1+q(7)); parameters%u=0.5_DP*(1+q(8)); status=0
    end subroutine

    subroutine ads_patch_mode_parameters(map,center,width,s0,ds,u0,du,status)
        type(ads_coordinate_map_type),intent(in) :: map
        real(DP),intent(in) :: center(:),width(:)
        real(DP),intent(out) :: s0,ds,u0,du
        integer,intent(out) :: status
        real(DP) :: tolerance
        s0=0; ds=0; u0=0; du=0; status=-1
        if(map%parameterization/=ADS_PARAM_SRP_MODES .or. map%n_variables/=8) return
        if(size(center)/=8 .or. size(width)/=8) return
        if(.not.all(ieee_is_finite(center)) .or. .not.all(ieee_is_finite(width))) return
        tolerance=1000*epsilon(1.0_DP)
        if(any(width<=0) .or. any(abs(center)+0.5_DP*width>1+tolerance)) return
        s0=0.5_DP*(1+center(7)); ds=0.25_DP*width(7)
        u0=0.5_DP*(1+center(8)); du=0.25_DP*width(8); status=0
    end subroutine

end module pod_uq_ads_coordinates_module
