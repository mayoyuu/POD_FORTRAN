program test_zero_covariance_sampling
    use pod_global, only: DP
    use pod_random_module, only: init_random_seed, generate_multivariate_normal
    implicit none

    integer, parameter :: dim = 6, n_samples = 4
    real(DP) :: mean(dim), cov(dim, dim), samples(dim, n_samples)
    real(DP) :: expected(dim, n_samples)

    mean = [1.0_DP, -2.0_DP, 3.0_DP, -4.0_DP, 5.0_DP, -6.0_DP]
    expected = spread(mean, dim=2, ncopies=n_samples)
    cov = 0.0_DP
    samples = huge(1.0_DP)

    call init_random_seed(.true.)
    call generate_multivariate_normal(mean, cov, samples)

    if (any(samples /= expected)) then
        write(*,*) 'FAIL: zero covariance did not return deterministic mean samples'
        error stop 1
    end if

    cov = 0.0_DP
    cov(1,1) = 1.0_DP
    cov(2,2) = 1.0_DP
    cov(3,3) = 1.0_DP
    cov(4,4) = 1.0_DP
    cov(5,5) = 1.0_DP
    cov(6,6) = 1.0_DP
    samples = expected

    call generate_multivariate_normal(mean, cov, samples)

    if (.not. any(samples /= expected)) then
        write(*,*) 'FAIL: positive-definite covariance did not sample deviations'
        error stop 1
    end if

    write(*,*) 'test_zero_covariance_sampling passed'
end program test_zero_covariance_sampling
