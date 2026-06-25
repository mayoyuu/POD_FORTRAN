program test_gmm_weight_update_formula
    use pod_global, only: DP
    use pod_gmm_math_module, only: compute_gmm_posterior_weights
    implicit none

    call test_eq34_normalization()
    call test_normalization_sums_to_one()

    write(*,*) 'test_gmm_weight_update_formula passed'

contains

    subroutine assert_close(actual, expected, tol, message)
        real(DP), intent(in) :: actual(:)
        real(DP), intent(in) :: expected(:)
        real(DP), intent(in) :: tol
        character(len=*), intent(in) :: message

        if (maxval(abs(actual - expected)) > tol) then
            write(*,*) trim(message)
            write(*,*) 'actual  = ', actual
            write(*,*) 'expected= ', expected
            error stop 1
        end if
    end subroutine assert_close

    subroutine assert_scalar_close(actual, expected, tol, message)
        real(DP), intent(in) :: actual
        real(DP), intent(in) :: expected
        real(DP), intent(in) :: tol
        character(len=*), intent(in) :: message

        if (abs(actual - expected) > tol) then
            write(*,*) trim(message), actual, expected
            error stop 1
        end if
    end subroutine assert_scalar_close

    subroutine test_eq34_normalization()
        real(DP) :: prior_weights(3), likelihoods(3), posterior(3)
        real(DP), parameter :: expected(3) = [0.06122448979591837_DP, &
            0.8163265306122449_DP, 0.12244897959183673_DP]
        real(DP), parameter :: tol = 1.0e-14_DP

        prior_weights = [0.3_DP, 0.4_DP, 0.3_DP]
        likelihoods   = [1.0_DP, 10.0_DP, 2.0_DP]

        call compute_gmm_posterior_weights(prior_weights, likelihoods, posterior)

        call assert_close(posterior, expected, tol, &
            'Eq. (34) normalization with known priors and likelihoods')
        call assert_scalar_close(sum(posterior), 1.0_DP, tol, &
            'posterior weights must sum to one')
    end subroutine test_eq34_normalization

    subroutine test_normalization_sums_to_one()
        real(DP) :: prior_weights(5), likelihoods(5), posterior(5)
        real(DP), parameter :: tol = 1.0e-14_DP
        integer :: i

        prior_weights = [0.1_DP, 0.2_DP, 0.1_DP, 0.3_DP, 0.3_DP]
        likelihoods   = [1.0e-10_DP, 1.0e-20_DP, 1.0_DP, 1.0e-5_DP, 1.0e-15_DP]

        call compute_gmm_posterior_weights(prior_weights, likelihoods, posterior)

        call assert_scalar_close(sum(posterior), 1.0_DP, tol, &
            'posterior weights sum to one with extreme likelihoods')
        do i = 1, 5
            if (posterior(i) < 0.0_DP .or. posterior(i) > 1.0_DP) then
                write(*,*) 'weight out of range: component', i, posterior(i)
                error stop 1
            end if
        end do
    end subroutine test_normalization_sums_to_one

end program test_gmm_weight_update_formula
