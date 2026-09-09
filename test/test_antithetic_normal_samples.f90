!> @file test_antithetic_normal_samples.f90
!! @brief Specify a reproducible antithetic standard-normal sample generator.
!!
!! Six-dimensional tensor Gauss-Hermite rules are too expensive for the
!! 15-day SRP study.  This test requires a deterministic local generator that
!! does not touch Fortran's global random-number state and exactly pairs each
    !! sample with its negative.  The finite ensemble is also moment matched,
    !! so every marginal variance must be unity to roundoff.
program test_antithetic_normal_samples
    use pod_global, only: DP
    use pod_uncertainty_diagnostics_module, only: &
        generate_antithetic_standard_normal_samples, compute_weighted_moments
    implicit none

    integer, parameter :: N_DIM=6, N_SAMPLES=4096
    real(DP) :: samples_a(N_DIM,N_SAMPLES),samples_b(N_DIM,N_SAMPLES)
    real(DP) :: odd_samples(2,3),weights(N_SAMPLES)
    real(DP) :: mean_value(N_DIM),covariance(N_DIM,N_DIM)
    real(DP) :: skewness(N_DIM),kurtosis(N_DIM),off_diagonal(N_DIM,N_DIM)
    integer :: status,i

    call generate_antithetic_standard_normal_samples(samples_a,260907,status)
    call assert_true(status==0,'valid generator status')
    call generate_antithetic_standard_normal_samples(samples_b,260907,status)
    call assert_true(status==0,'repeat generator status')
    call assert_true(all(samples_a==samples_b),'same seed is reproducible')
    call assert_close(maxval(abs(samples_a(:,1:N_SAMPLES/2)+ &
        samples_a(:,N_SAMPLES/2+1:N_SAMPLES))),0.0_DP,0.0_DP, &
        'samples are exact antithetic pairs')

    weights=1.0_DP/real(N_SAMPLES,DP)
    call compute_weighted_moments(samples_a,weights,mean_value,covariance, &
                                  skewness,kurtosis,status)
    call assert_true(status==0,'sample moment status')
    call assert_close(maxval(abs(mean_value)),0.0_DP,1.0e-14_DP, &
                      'antithetic mean')
    do i=1,N_DIM
        call assert_close(covariance(i,i),1.0_DP,1.0e-12_DP, &
                          'unit marginal variance')
    end do
    off_diagonal=covariance
    do i=1,N_DIM
        off_diagonal(i,i)=0.0_DP
    end do
    call assert_true(maxval(abs(off_diagonal))<8.0e-2_DP, &
                     'small off-diagonal covariance')
    call assert_true(maxval(abs(skewness))<1.0e-12_DP, &
                     'antithetic marginal skewness')
    call assert_true(maxval(abs(kurtosis))<2.0e-1_DP, &
                     'normal marginal excess kurtosis')

    call generate_antithetic_standard_normal_samples(odd_samples,260907,status)
    call assert_true(status==-1,'odd sample count is rejected')
    call generate_antithetic_standard_normal_samples(samples_b,0,status)
    call assert_true(status==-2,'non-positive seed is rejected')

    write(*,'(a)') 'PASS: reproducible antithetic standard-normal samples.'

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
            write(*,'(a,3(1x,es24.16))') 'FAIL: '//trim(label), &
                actual,expected,tolerance
            error stop 1
        end if
    end subroutine assert_close

end program test_antithetic_normal_samples
