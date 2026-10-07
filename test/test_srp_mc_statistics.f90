program test_srp_mc_statistics
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use pod_global, only: DP
    use pod_srp_mc_statistics_module, only: srp_mc_group_stats_type, &
        srp_mc_pair_stats_type, compute_group_stats, compute_pair_stats
    implicit none

    type(srp_mc_group_stats_type) :: treatment_stats, control_stats
    type(srp_mc_pair_stats_type) :: pair_stats
    real(DP) :: treatment(6,5), control(6,5), nominal(6)
    real(DP), parameter :: SQRT3=sqrt(3.0_DP)
    integer :: i, status

    nominal=0.0_DP
    treatment=0.0_DP
    control=0.0_DP
    treatment(:,5)=5.0_DP
    control(:,5)=2.5_DP

    call compute_group_stats(treatment,nominal,treatment_stats,status)
    call assert_true(status==0,'treatment status')
    call assert_true(treatment_stats%sample_count==5,'sample count')
    call assert_vec(treatment_stats%mean_error,[1.0_DP,1.0_DP,1.0_DP, &
        1.0_DP,1.0_DP,1.0_DP],1.0e-13_DP,'mean')
    call assert_vec(treatment_stats%component_min,[0.0_DP,0.0_DP,0.0_DP, &
        0.0_DP,0.0_DP,0.0_DP],1.0e-13_DP,'minimum')
    call assert_vec(treatment_stats%component_max,[5.0_DP,5.0_DP,5.0_DP, &
        5.0_DP,5.0_DP,5.0_DP],1.0e-13_DP,'maximum')
    do i=1,6
        call assert_close(treatment_stats%covariance(i,i),5.0_DP,1.0e-13_DP,'variance')
        call assert_close(treatment_stats%skewness(i),sqrt(5.0_DP),1.0e-13_DP,'skewness')
        call assert_close(treatment_stats%excess_kurtosis(i),5.0_DP,1.0e-13_DP,'kurtosis')
    end do
    call assert_vec(treatment_stats%position_axes_rms,[sqrt(15.0_DP),0.0_DP,0.0_DP], &
        1.0e-12_DP,'position axes')
    call assert_vec(treatment_stats%velocity_axes_rms,[sqrt(15.0_DP),0.0_DP,0.0_DP], &
        1.0e-12_DP,'velocity axes')
    call assert_close(treatment_stats%position_norm_rms,sqrt(15.0_DP),1.0e-13_DP,'position RMS')
    call assert_close(treatment_stats%position_norm_quantiles(1),0.0_DP,1.0e-13_DP,'p50')
    call assert_close(treatment_stats%position_norm_quantiles(2),4.0_DP*SQRT3,1.0e-12_DP,'p95')
    call assert_close(treatment_stats%position_norm_quantiles(3),4.8_DP*SQRT3,1.0e-12_DP,'p99')
    call assert_close(treatment_stats%position_norm_quantiles(4),4.98_DP*SQRT3,1.0e-12_DP,'p999')
    call assert_close(treatment_stats%position_norm_max,5.0_DP*SQRT3,1.0e-12_DP,'position max')
    call assert_true(treatment_stats%position_worst_id==5,'position worst id')
    call assert_true(treatment_stats%velocity_worst_id==5,'velocity worst id')

    call compute_group_stats(control,nominal,control_stats,status)
    call assert_true(status==0,'control status')
    call compute_pair_stats(treatment,control,treatment_stats,control_stats,pair_stats,status)
    call assert_true(status==0,'pair status')
    call assert_vec(pair_stats%mean_delta,[0.5_DP,0.5_DP,0.5_DP, &
        0.5_DP,0.5_DP,0.5_DP],1.0e-13_DP,'mean delta')
    call assert_close(pair_stats%covariance(1,1),1.25_DP,1.0e-13_DP,'delta covariance')
    call assert_close(pair_stats%position_delta_rms,sqrt(3.75_DP),1.0e-13_DP,'delta r rms')
    call assert_close(pair_stats%position_delta_quantiles(2),2.0_DP*SQRT3,1.0e-12_DP,'delta p95')
    call assert_close(pair_stats%position_delta_max,2.5_DP*SQRT3,1.0e-12_DP,'delta max')
    call assert_true(pair_stats%position_worst_id==5,'delta worst id')

    treatment(1,1)=ieee_value(0.0_DP,ieee_quiet_nan)
    call compute_group_stats(treatment,nominal,treatment_stats,status)
    call assert_true(status/=0,'nonfinite rejected')

    write(*,'(a)') 'PASS: SRP MC statistics.'

contains

    subroutine assert_true(condition,label)
        logical,intent(in) :: condition
        character(len=*),intent(in) :: label
        if(.not.condition) then
            write(*,'(a)') 'FAIL: '//trim(label)
            error stop 1
        end if
    end subroutine assert_true

    subroutine assert_close(actual,expected,tolerance,label)
        real(DP),intent(in) :: actual,expected,tolerance
        character(len=*),intent(in) :: label
        if(abs(actual-expected)>tolerance) then
            write(*,'(a,3(1x,es24.16))') 'FAIL: '//trim(label),actual,expected,tolerance
            error stop 1
        end if
    end subroutine assert_close

    subroutine assert_vec(actual,expected,tolerance,label)
        real(DP),intent(in) :: actual(:),expected(:),tolerance
        character(len=*),intent(in) :: label
        if(size(actual)/=size(expected).or.maxval(abs(actual-expected))>tolerance) then
            write(*,'(a,1x,es24.16)') 'FAIL: '//trim(label),maxval(abs(actual-expected))
            error stop 1
        end if
    end subroutine assert_vec

end program test_srp_mc_statistics
