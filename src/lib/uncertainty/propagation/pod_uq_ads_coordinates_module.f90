!> Coordinate maps between physical uncertainty and the ADS unit box.
module pod_uq_ads_coordinates_module
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use pod_global, only: DP
    use pod_basicmath_module, only: dpotrf
    implicit none
    private

    integer, parameter, public :: ADS_COORD_COMPONENT = 1
    integer, parameter, public :: ADS_COORD_WHITENED = 2

    type, public :: ads_coordinate_map_type
        integer :: mode = ADS_COORD_COMPONENT
        integer :: n_variables = 6
        real(DP) :: domain_sigma = 3.0_DP
        real(DP) :: srp_sigma = 0.0_DP
        real(DP) :: mean6(6) = 0.0_DP
        real(DP) :: basis6(6,6) = 0.0_DP
    end type ads_coordinate_map_type

    public :: ads_build_coordinate_map
    public :: ads_physical_to_unit
    public :: ads_unit_to_physical

contains

    subroutine ads_build_coordinate_map(mean6,covariance,mode,domain_sigma, &
                                        srp_sigma,map,status,message)
        real(DP), intent(in) :: mean6(6),covariance(6,6)
        integer, intent(in) :: mode
        real(DP), intent(in) :: domain_sigma,srp_sigma
        type(ads_coordinate_map_type), intent(out) :: map
        integer, intent(out) :: status
        character(len=*), intent(out) :: message
        real(DP) :: symmetric(6,6),scale,symmetry_tolerance
        integer :: i,j,info

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

        scale=max(1.0_DP,maxval(abs(covariance)))
        symmetry_tolerance=1000.0_DP*epsilon(1.0_DP)*scale
        if(maxval(abs(covariance-transpose(covariance)))>symmetry_tolerance) then
            call fail(-5,'ADS covariance must be symmetric')
            return
        end if

        map%mode=mode
        map%n_variables=merge(7,6,srp_sigma>0.0_DP)
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

    subroutine ads_physical_to_unit(map,deviation,eta_srp,unit_point,status)
        type(ads_coordinate_map_type), intent(in) :: map
        real(DP), intent(in) :: deviation(6),eta_srp
        real(DP), allocatable, intent(out) :: unit_point(:)
        integer, intent(out) :: status
        real(DP) :: fixed_tolerance,srp_scale
        integer :: i,j

        status=0
        allocate(unit_point(map%n_variables))
        unit_point=0.0_DP
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
    end subroutine ads_physical_to_unit

    subroutine ads_unit_to_physical(map,unit_point,deviation,eta_srp,status)
        type(ads_coordinate_map_type), intent(in) :: map
        real(DP), intent(in) :: unit_point(:)
        real(DP), intent(out) :: deviation(6),eta_srp
        integer, intent(out) :: status

        status=0
        deviation=0.0_DP
        eta_srp=0.0_DP
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
    end subroutine ads_unit_to_physical

end module pod_uq_ads_coordinates_module
