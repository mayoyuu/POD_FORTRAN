program test_covariance_condition_number
    use pod_global, only: DP
    use pod_error_analysis_module, only: compute_covariance_condition_number
    implicit none

    real(DP) :: cov(6,6)
    real(DP) :: cond
    integer :: info

    cov = 0.0_DP
    cov(1,1) = 1.0_DP
    cov(2,2) = 2.0_DP
    cov(3,3) = 4.0_DP
    cov(4,4) = 8.0_DP
    cov(5,5) = 16.0_DP
    cov(6,6) = 32.0_DP

    call compute_covariance_condition_number(cov, cond, info)

    if (info /= 0) error stop "positive diagonal covariance should be valid"
    if (abs(cond - 32.0_DP) > 1.0e-10_DP) then
        write(*,*) "Expected condition number 32, got ", cond
        error stop "unexpected covariance condition number"
    end if

    cov(6,6) = 0.0_DP
    call compute_covariance_condition_number(cov, cond, info)

    if (info == 0) error stop "singular covariance should be flagged"
    if (cond < huge(1.0_DP) / 2.0_DP) error stop "singular covariance should report huge condition number"

    write(*,*) "test_covariance_condition_number passed"
end program test_covariance_condition_number

