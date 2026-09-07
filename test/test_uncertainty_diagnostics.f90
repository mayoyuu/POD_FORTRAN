!> @file test_uncertainty_diagnostics.f90
!! @brief Analytic tests for deterministic quadrature and weighted moments.
!!
!! These tests contain no orbit dynamics.  They establish that the statistics
!! later written by the 15-day SRP test have the intended probability measure,
!! covariance convention, higher moments, ordered axes and effective rank.
program test_uncertainty_diagnostics
    use pod_global, only: DP
    use pod_uncertainty_diagnostics_module, only: &
        gauss_legendre_probability_rule, gauss_normal_probability_rule, &
        compute_weighted_moments, compute_covariance_axes, compute_effective_rank
    implicit none

    integer, parameter :: N = 9
    real(DP), parameter :: TOL = 2.0e-13_DP
    real(DP) :: nodes(N), weights(N), samples(2,N)
    real(DP) :: mean(2), covariance(2,2), skewness(2), excess_kurtosis(2)
    real(DP) :: p3(3,3), eigenvalues(3), eigenvectors(3,3), axis_sigma(3)
    real(DP) :: m1, m2, m4
    integer :: status, rank_value, i

    ! Nine-point probability Gauss-Legendre integrates these moments exactly.
    call gauss_legendre_probability_rule(nodes, weights, status)
    call assert_true(status == 0, 'Gauss-Legendre rule status')
    m1 = sum(weights*nodes)
    m2 = sum(weights*nodes**2)
    m4 = sum(weights*nodes**4)
    call assert_close(sum(weights), 1.0_DP, TOL, 'uniform weight sum')
    call assert_close(m1, 0.0_DP, TOL, 'uniform first moment')
    call assert_close(m2, 1.0_DP/3.0_DP, TOL, 'uniform second moment')
    call assert_close(m4, 1.0_DP/5.0_DP, TOL, 'uniform fourth moment')

    ! For q~N(0,1/3), E(q^4)=3 sigma^4=1/3.
    call gauss_normal_probability_rule(1.0_DP/sqrt(3.0_DP), nodes, weights, status)
    call assert_true(status == 0, 'Gauss-Hermite rule status')
    m1 = sum(weights*nodes)
    m2 = sum(weights*nodes**2)
    m4 = sum(weights*nodes**4)
    call assert_close(sum(weights), 1.0_DP, TOL, 'normal weight sum')
    call assert_close(m1, 0.0_DP, TOL, 'normal first moment')
    call assert_close(m2, 1.0_DP/3.0_DP, TOL, 'normal second moment')
    call assert_close(m4, 1.0_DP/3.0_DP, 5.0e-13_DP, 'normal fourth moment')

    ! The affine map x=[2q,-q] has an analytic rank-one covariance.
    call gauss_legendre_probability_rule(nodes, weights, status)
    do i = 1, N
        samples(:,i) = [2.0_DP*nodes(i), -nodes(i)]
    end do
    call compute_weighted_moments(samples, weights, mean, covariance, &
                                  skewness, excess_kurtosis, status)
    call assert_true(status == 0, 'weighted moment status')
    call assert_vector2_close(mean, [0.0_DP,0.0_DP], TOL, 'affine mean')
    call assert_matrix2_close(covariance, reshape([4.0_DP/3.0_DP, -2.0_DP/3.0_DP, &
                                                  -2.0_DP/3.0_DP, 1.0_DP/3.0_DP], [2,2]), &
                              5.0e-13_DP, 'affine covariance')
    call assert_vector2_close(skewness, [0.0_DP,0.0_DP], 2.0e-12_DP, 'uniform skewness')
    call assert_vector2_close(excess_kurtosis, [-1.2_DP,-1.2_DP], 2.0e-12_DP, &
                              'uniform excess kurtosis')
    call compute_effective_rank(covariance, 1.0e-12_DP, rank_value, status)
    call assert_true(status == 0 .and. rank_value == 1, 'rank-one covariance')

    ! Covariance axes are sorted from longest to shortest.
    p3 = 0.0_DP
    p3(1,1) = 1.0_DP
    p3(2,2) = 9.0_DP
    p3(3,3) = 4.0_DP
    call compute_covariance_axes(p3, eigenvalues, eigenvectors, axis_sigma, status)
    call assert_true(status == 0, 'covariance axes status')
    call assert_vector3_close(eigenvalues, [9.0_DP,4.0_DP,1.0_DP], TOL, &
                              'descending covariance eigenvalues')
    call assert_vector3_close(axis_sigma, [3.0_DP,2.0_DP,1.0_DP], TOL, &
                              'one-sigma semi-axis scales')
    call assert_matrix3_close(matmul(transpose(eigenvectors),eigenvectors), &
                              identity3(), 5.0e-13_DP, 'eigenvector orthonormality')

    write(*,'(a)') 'PASS: deterministic quadrature and uncertainty diagnostics.'

contains

    pure function identity3() result(identity)
        real(DP) :: identity(3,3)
        identity = 0.0_DP
        identity(1,1)=1.0_DP
        identity(2,2)=1.0_DP
        identity(3,3)=1.0_DP
    end function identity3

    subroutine assert_true(condition, label)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: label
        if (.not. condition) then
            write(*,'(a)') 'FAIL: '//trim(label)
            error stop 1
        end if
    end subroutine assert_true

    subroutine assert_close(actual, expected, tolerance, label)
        real(DP), intent(in) :: actual, expected, tolerance
        character(len=*), intent(in) :: label
        if (abs(actual-expected) > tolerance) then
            write(*,'(a,3(1x,es24.16))') 'FAIL: '//trim(label), actual, expected, tolerance
            error stop 1
        end if
    end subroutine assert_close

    subroutine assert_vector2_close(actual, expected, tolerance, label)
        real(DP), intent(in) :: actual(2), expected(2), tolerance
        character(len=*), intent(in) :: label
        if (maxval(abs(actual-expected)) > tolerance) then
            write(*,'(a,1x,es24.16)') 'FAIL: '//trim(label), maxval(abs(actual-expected))
            error stop 1
        end if
    end subroutine assert_vector2_close

    subroutine assert_vector3_close(actual, expected, tolerance, label)
        real(DP), intent(in) :: actual(3), expected(3), tolerance
        character(len=*), intent(in) :: label
        if (maxval(abs(actual-expected)) > tolerance) then
            write(*,'(a,1x,es24.16)') 'FAIL: '//trim(label), maxval(abs(actual-expected))
            error stop 1
        end if
    end subroutine assert_vector3_close

    subroutine assert_matrix2_close(actual, expected, tolerance, label)
        real(DP), intent(in) :: actual(2,2), expected(2,2), tolerance
        character(len=*), intent(in) :: label
        if (maxval(abs(actual-expected)) > tolerance) then
            write(*,'(a,1x,es24.16)') 'FAIL: '//trim(label), maxval(abs(actual-expected))
            error stop 1
        end if
    end subroutine assert_matrix2_close

    subroutine assert_matrix3_close(actual, expected, tolerance, label)
        real(DP), intent(in) :: actual(3,3), expected(3,3), tolerance
        character(len=*), intent(in) :: label
        if (maxval(abs(actual-expected)) > tolerance) then
            write(*,'(a,1x,es24.16)') 'FAIL: '//trim(label), maxval(abs(actual-expected))
            error stop 1
        end if
    end subroutine assert_matrix3_close

end program test_uncertainty_diagnostics

