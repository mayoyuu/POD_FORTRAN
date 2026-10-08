!> Statistics for fixed, point-paired SRP Monte Carlo evaluations.
module pod_srp_mc_statistics_module
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use pod_global, only: DP
    use pod_uncertainty_diagnostics_module, only: compute_covariance_axes
    implicit none
    private

    integer, parameter :: STATE_DIM=6, QUANTILE_COUNT=4
    real(DP), parameter :: QUANTILE_PROBABILITIES(QUANTILE_COUNT)= &
        [0.5_DP,0.95_DP,0.99_DP,0.999_DP]

    type, public :: srp_mc_group_stats_type
        integer :: sample_count=0
        real(DP) :: mean_error(STATE_DIM)=0.0_DP
        real(DP) :: covariance(STATE_DIM,STATE_DIM)=0.0_DP
        real(DP) :: standard_deviation(STATE_DIM)=0.0_DP
        real(DP) :: skewness(STATE_DIM)=0.0_DP
        real(DP) :: excess_kurtosis(STATE_DIM)=0.0_DP
        real(DP) :: component_min(STATE_DIM)=0.0_DP
        real(DP) :: component_max(STATE_DIM)=0.0_DP
        real(DP) :: position_axes_rms(3)=0.0_DP
        real(DP) :: velocity_axes_rms(3)=0.0_DP
        real(DP) :: position_norm_rms=0.0_DP
        real(DP) :: velocity_norm_rms=0.0_DP
        real(DP) :: position_norm_quantiles(QUANTILE_COUNT)=0.0_DP
        real(DP) :: velocity_norm_quantiles(QUANTILE_COUNT)=0.0_DP
        real(DP) :: position_norm_max=0.0_DP
        real(DP) :: velocity_norm_max=0.0_DP
        integer :: position_worst_id=0
        integer :: velocity_worst_id=0
        real(DP) :: position_max_to_longest_axis=0.0_DP
        real(DP) :: velocity_max_to_longest_axis=0.0_DP
    end type srp_mc_group_stats_type

    type, public :: srp_mc_pair_stats_type
        integer :: sample_count=0
        real(DP) :: mean_delta(STATE_DIM)=0.0_DP
        real(DP) :: covariance(STATE_DIM,STATE_DIM)=0.0_DP
        real(DP) :: standard_deviation(STATE_DIM)=0.0_DP
        real(DP) :: component_delta_min(STATE_DIM)=0.0_DP
        real(DP) :: component_delta_max(STATE_DIM)=0.0_DP
        real(DP) :: position_delta_rms=0.0_DP
        real(DP) :: velocity_delta_rms=0.0_DP
        real(DP) :: position_delta_quantiles(QUANTILE_COUNT)=0.0_DP
        real(DP) :: velocity_delta_quantiles(QUANTILE_COUNT)=0.0_DP
        real(DP) :: position_delta_max=0.0_DP
        real(DP) :: velocity_delta_max=0.0_DP
        integer :: position_worst_id=0
        integer :: velocity_worst_id=0
        real(DP) :: position_max_increment=0.0_DP
        real(DP) :: velocity_max_increment=0.0_DP
        real(DP) :: position_max_amplification=0.0_DP
        real(DP) :: velocity_max_amplification=0.0_DP
        real(DP) :: longest_position_axis_increment=0.0_DP
        real(DP) :: longest_velocity_axis_increment=0.0_DP
        real(DP) :: longest_position_axis_amplification=0.0_DP
        real(DP) :: longest_velocity_axis_amplification=0.0_DP
        real(DP) :: delta_position_max_to_treatment_axis=0.0_DP
        real(DP) :: delta_position_max_to_control_axis=0.0_DP
        real(DP) :: delta_velocity_max_to_treatment_axis=0.0_DP
        real(DP) :: delta_velocity_max_to_control_axis=0.0_DP
    end type srp_mc_pair_stats_type

    public :: compute_group_stats, compute_pair_stats

contains

    subroutine compute_group_stats(states,nominal,stats,status)
        real(DP),intent(in) :: states(:,:),nominal(:)
        type(srp_mc_group_stats_type),intent(out) :: stats
        integer,intent(out) :: status

        real(DP),allocatable :: errors(:,:),position_norm(:),velocity_norm(:)
        real(DP) :: eigenvalues(3),eigenvectors(3,3)
        integer :: n,axis_status

        stats=srp_mc_group_stats_type()
        status=0
        if(size(states,1)/=STATE_DIM .or. size(nominal)/=STATE_DIM) then
            status=-1
            return
        end if
        n=size(states,2)
        if(n<2) then
            status=-2
            return
        end if
        if(any(.not.ieee_is_finite(states)) .or. any(.not.ieee_is_finite(nominal))) then
            status=-3
            return
        end if

        allocate(errors(STATE_DIM,n),position_norm(n),velocity_norm(n))
        errors=states-spread(nominal,dim=2,ncopies=n)
        stats%sample_count=n
        call component_moments(errors,stats%mean_error,stats%covariance, &
            stats%skewness,stats%excess_kurtosis)
        stats%standard_deviation=sqrt(max([(stats%covariance(axis_status,axis_status), &
            axis_status=1,STATE_DIM)],0.0_DP))
        stats%component_min=minval(errors,dim=2)
        stats%component_max=maxval(errors,dim=2)

        call compute_covariance_axes(stats%covariance(1:3,1:3),eigenvalues, &
            eigenvectors,stats%position_axes_rms,axis_status)
        if(axis_status/=0) then
            status=-4
            deallocate(errors,position_norm,velocity_norm)
            return
        end if
        call compute_covariance_axes(stats%covariance(4:6,4:6),eigenvalues, &
            eigenvectors,stats%velocity_axes_rms,axis_status)
        if(axis_status/=0) then
            status=-5
            deallocate(errors,position_norm,velocity_norm)
            return
        end if

        position_norm=sqrt(sum(errors(1:3,:)**2,dim=1))
        velocity_norm=sqrt(sum(errors(4:6,:)**2,dim=1))
        call fill_norm_stats(position_norm,stats%position_norm_rms, &
            stats%position_norm_quantiles,stats%position_norm_max, &
            stats%position_worst_id)
        call fill_norm_stats(velocity_norm,stats%velocity_norm_rms, &
            stats%velocity_norm_quantiles,stats%velocity_norm_max, &
            stats%velocity_worst_id)
        stats%position_max_to_longest_axis=safe_ratio(stats%position_norm_max, &
            stats%position_axes_rms(1))
        stats%velocity_max_to_longest_axis=safe_ratio(stats%velocity_norm_max, &
            stats%velocity_axes_rms(1))
        deallocate(errors,position_norm,velocity_norm)
    end subroutine compute_group_stats

    subroutine compute_pair_stats(treatment,control,treatment_stats,control_stats, &
                                  stats,status)
        real(DP),intent(in) :: treatment(:,:),control(:,:)
        type(srp_mc_group_stats_type),intent(in) :: treatment_stats,control_stats
        type(srp_mc_pair_stats_type),intent(out) :: stats
        integer,intent(out) :: status

        real(DP),allocatable :: delta(:,:),position_norm(:),velocity_norm(:)
        real(DP) :: correction(STATE_DIM),dummy_skew(STATE_DIM),dummy_kurtosis(STATE_DIM)
        integer :: n,i,j

        stats=srp_mc_pair_stats_type()
        status=0
        if(size(treatment,1)/=STATE_DIM .or. size(control,1)/=STATE_DIM .or. &
           size(treatment,2)/=size(control,2)) then
            status=-1
            return
        end if
        n=size(treatment,2)
        if(n<2 .or. treatment_stats%sample_count/=n .or. &
           control_stats%sample_count/=n) then
            status=-2
            return
        end if
        if(any(.not.ieee_is_finite(treatment)) .or. &
           any(.not.ieee_is_finite(control))) then
            status=-3
            return
        end if

        allocate(delta(STATE_DIM,n),position_norm(n),velocity_norm(n))
        delta=treatment-control
        stats%sample_count=n
        stats%mean_delta=0.0_DP
        correction=0.0_DP
        do j=1,n
            do i=1,STATE_DIM
                call kahan_add(stats%mean_delta(i),correction(i),delta(i,j))
            end do
        end do
        stats%mean_delta=stats%mean_delta/real(n,DP)
        call component_moments(delta,stats%mean_delta,stats%covariance, &
            dummy_skew,dummy_kurtosis)
        stats%standard_deviation=sqrt(max([(stats%covariance(i,i),i=1,STATE_DIM)],0.0_DP))
        stats%component_delta_min=minval(delta,dim=2)
        stats%component_delta_max=maxval(delta,dim=2)

        position_norm=sqrt(sum(delta(1:3,:)**2,dim=1))
        velocity_norm=sqrt(sum(delta(4:6,:)**2,dim=1))
        call fill_norm_stats(position_norm,stats%position_delta_rms, &
            stats%position_delta_quantiles,stats%position_delta_max, &
            stats%position_worst_id)
        call fill_norm_stats(velocity_norm,stats%velocity_delta_rms, &
            stats%velocity_delta_quantiles,stats%velocity_delta_max, &
            stats%velocity_worst_id)

        stats%position_max_increment=treatment_stats%position_norm_max- &
            control_stats%position_norm_max
        stats%velocity_max_increment=treatment_stats%velocity_norm_max- &
            control_stats%velocity_norm_max
        stats%longest_position_axis_increment=treatment_stats%position_axes_rms(1)- &
            control_stats%position_axes_rms(1)
        stats%longest_velocity_axis_increment=treatment_stats%velocity_axes_rms(1)- &
            control_stats%velocity_axes_rms(1)
        stats%position_max_amplification=safe_ratio(treatment_stats%position_norm_max, &
            control_stats%position_norm_max)
        stats%velocity_max_amplification=safe_ratio(treatment_stats%velocity_norm_max, &
            control_stats%velocity_norm_max)
        stats%longest_position_axis_amplification=safe_ratio( &
            treatment_stats%position_axes_rms(1),control_stats%position_axes_rms(1))
        stats%longest_velocity_axis_amplification=safe_ratio( &
            treatment_stats%velocity_axes_rms(1),control_stats%velocity_axes_rms(1))
        stats%delta_position_max_to_treatment_axis=safe_ratio( &
            stats%position_delta_max,treatment_stats%position_axes_rms(1))
        stats%delta_position_max_to_control_axis=safe_ratio( &
            stats%position_delta_max,control_stats%position_axes_rms(1))
        stats%delta_velocity_max_to_treatment_axis=safe_ratio( &
            stats%velocity_delta_max,treatment_stats%velocity_axes_rms(1))
        stats%delta_velocity_max_to_control_axis=safe_ratio( &
            stats%velocity_delta_max,control_stats%velocity_axes_rms(1))
        deallocate(delta,position_norm,velocity_norm)
    end subroutine compute_pair_stats

    subroutine component_moments(values,mean_value,covariance,skewness,kurtosis)
        real(DP),intent(in) :: values(:,:)
        real(DP),intent(out) :: mean_value(:),covariance(:,:),skewness(:),kurtosis(:)
        real(DP) :: mean_correction(STATE_DIM),cov_correction(STATE_DIM,STATE_DIM)
        real(DP) :: third(STATE_DIM),fourth(STATE_DIM)
        real(DP) :: third_correction(STATE_DIM),fourth_correction(STATE_DIM)
        real(DP) :: deviation(STATE_DIM),variance,sigma,scale,m2,g1,g2
        integer :: n,i,k,j

        n=size(values,2)
        mean_value=0.0_DP
        mean_correction=0.0_DP
        do j=1,n
            do i=1,STATE_DIM
                call kahan_add(mean_value(i),mean_correction(i),values(i,j))
            end do
        end do
        mean_value=mean_value/real(n,DP)

        covariance=0.0_DP
        cov_correction=0.0_DP
        third=0.0_DP
        fourth=0.0_DP
        third_correction=0.0_DP
        fourth_correction=0.0_DP
        do j=1,n
            deviation=values(:,j)-mean_value
            do i=1,STATE_DIM
                do k=1,STATE_DIM
                    call kahan_add(covariance(i,k),cov_correction(i,k), &
                        deviation(i)*deviation(k))
                end do
                call kahan_add(third(i),third_correction(i),deviation(i)**3)
                call kahan_add(fourth(i),fourth_correction(i),deviation(i)**4)
            end do
        end do
        covariance=covariance/real(n-1,DP)
        covariance=0.5_DP*(covariance+transpose(covariance))
        skewness=0.0_DP
        kurtosis=0.0_DP
        do i=1,STATE_DIM
            variance=max(covariance(i,i),0.0_DP)
            m2=variance*real(n-1,DP)/real(n,DP)
            sigma=sqrt(m2)
            scale=max(1.0_DP,maxval(abs(values(i,:))),abs(mean_value(i)))
            if(sigma>100.0_DP*epsilon(1.0_DP)*scale) then
                g1=(third(i)/real(n,DP))/sigma**3
                g2=(fourth(i)/real(n,DP))/m2**2-3.0_DP
                if(n>2) skewness(i)=sqrt(real(n,DP)*real(n-1,DP))/ &
                    real(n-2,DP)*g1
                if(n>3) kurtosis(i)=real(n-1,DP)/ &
                    (real(n-2,DP)*real(n-3,DP))* &
                    (real(n+1,DP)*g2+6.0_DP)
            end if
        end do
    end subroutine component_moments

    subroutine fill_norm_stats(values,rms,quantiles,maximum,worst_id)
        real(DP),intent(in) :: values(:)
        real(DP),intent(out) :: rms,quantiles(:),maximum
        integer,intent(out) :: worst_id
        real(DP),allocatable :: ordered(:)
        real(DP) :: sum_square,correction
        integer :: i

        allocate(ordered(size(values)))
        ordered=values
        sum_square=0.0_DP
        correction=0.0_DP
        do i=1,size(values)
            call kahan_add(sum_square,correction,values(i)*values(i))
        end do
        rms=sqrt(sum_square/real(size(values),DP))
        maximum=maxval(values)
        worst_id=maxloc(values,dim=1)
        call quicksort(ordered,1,size(ordered))
        do i=1,size(quantiles)
            quantiles(i)=linear_quantile(ordered,QUANTILE_PROBABILITIES(i))
        end do
        deallocate(ordered)
    end subroutine fill_norm_stats

    pure function linear_quantile(ordered,probability) result(value)
        real(DP),intent(in) :: ordered(:),probability
        real(DP) :: value,index_value,fraction
        integer :: lower_index,upper_index

        index_value=probability*real(size(ordered)-1,DP)
        lower_index=1+floor(index_value)
        upper_index=min(lower_index+1,size(ordered))
        fraction=index_value-floor(index_value)
        value=(1.0_DP-fraction)*ordered(lower_index)+fraction*ordered(upper_index)
    end function linear_quantile

    recursive subroutine quicksort(values,left,right)
        real(DP),intent(inout) :: values(:)
        integer,intent(in) :: left,right
        integer :: i,j
        real(DP) :: pivot,temp

        if(left>=right) return
        i=left
        j=right
        pivot=values((left+right)/2)
        do
            do while(values(i)<pivot)
                i=i+1
            end do
            do while(values(j)>pivot)
                j=j-1
            end do
            if(i<=j) then
                temp=values(i)
                values(i)=values(j)
                values(j)=temp
                i=i+1
                j=j-1
            end if
            if(i>j) exit
        end do
        if(left<j) call quicksort(values,left,j)
        if(i<right) call quicksort(values,i,right)
    end subroutine quicksort

    pure subroutine kahan_add(total,correction,value)
        real(DP),intent(inout) :: total,correction
        real(DP),intent(in) :: value
        real(DP) :: adjusted,updated

        adjusted=value-correction
        updated=total+adjusted
        correction=(updated-total)-adjusted
        total=updated
    end subroutine kahan_add

    pure function safe_ratio(numerator,denominator) result(value)
        real(DP),intent(in) :: numerator,denominator
        real(DP) :: value

        if(abs(denominator)>tiny(1.0_DP)) then
            value=numerator/denominator
        else if(abs(numerator)<=tiny(1.0_DP)) then
            value=1.0_DP
        else
            value=huge(1.0_DP)
        end if
    end function safe_ratio

end module pod_srp_mc_statistics_module
