program test_gmm_weight_floor
    use pod_global, only: DP
    use pod_gmm_math_module, only: apply_gmm_weight_floor
    implicit none

    call test_floor_lifts_zero_weights()
    call test_zero_floor_only_normalizes()

    write(*,*) 'test_gmm_weight_floor passed'

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

    subroutine test_floor_lifts_zero_weights()
        real(DP) :: weights(3)
        real(DP) :: expected(3)
        real(DP), parameter :: floor = 1.0e-4_DP
        real(DP), parameter :: tol = 1.0e-12_DP

        weights = [1.0_DP, 0.0_DP, 0.0_DP]
        expected = [1.0_DP, floor, floor] / (1.0_DP + 2.0_DP * floor)

        call apply_gmm_weight_floor(weights, floor)

        call assert_close(weights, expected, tol, 'floor should lift zero weights and renormalize')
        call assert_scalar_close(sum(weights), 1.0_DP, tol, 'floored weights should sum to one')
    end subroutine test_floor_lifts_zero_weights

    subroutine test_zero_floor_only_normalizes()
        real(DP) :: weights(3)
        real(DP) :: expected(3)
        real(DP), parameter :: tol = 1.0e-12_DP

        weights = [2.0_DP, 3.0_DP, 5.0_DP]
        expected = [0.2_DP, 0.3_DP, 0.5_DP]

        call apply_gmm_weight_floor(weights, 0.0_DP)

        call assert_close(weights, expected, tol, 'zero floor should preserve weights up to normalization')
    end subroutine test_zero_floor_only_normalizes

end program test_gmm_weight_floor