! Regression for Real PlmX leaving the P(n,n+1) sentinel undefined.
program test_gravity_legendre_tail
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_finite
    use pod_global, only: DP
    use pod_gravity_model_module, only: PlmX, PlmX_da, GravityTempPool
    use pod_dace_classes, only: DA, dace_initialize, assignment(=)
    implicit none
    integer :: i
    real(DP), parameter :: latitudes(3) = [-0.7_DP, 0.0_DP, 0.6_DP]

    call dace_initialize(2, 1)
    do i = 1, size(latitudes)
        call check_degree(2, latitudes(i))
        call check_degree(10, latitudes(i))
    end do
    print *, 'PASS: Real Legendre output including the derivative sentinel matches DA'
contains
    subroutine check_degree(n, zr)
        integer, intent(in) :: n
        real(DP), intent(in) :: zr
        real(DP) :: actual(n,n+1), expected
        type(DA) :: zr_da, expected_da(n,n+1)
        type(GravityTempPool) :: pool
        integer :: row, col

        actual = ieee_value(0.0_DP, ieee_quiet_nan)
        call PlmX(n, zr, actual)
        if (.not. all(ieee_is_finite(actual))) then
            print *, 'FAIL: Real Legendre output contains an undefined sentinel; degree=', n
            error stop 1
        end if
        if (any(actual(:,n+1) /= 0.0_DP)) error stop 'FAIL: P(l,n+1) must be zero'

        call pool%init(3)
        call zr_da%init()
        zr_da = zr
        call PlmX_da(n, zr_da, expected_da, pool)
        do col = 1, n+1
            do row = 1, n
                expected = expected_da(row,col)%cons()
                if (abs(actual(row,col)-expected) > 1.0e-13_DP*max(1.0_DP,abs(expected))) &
                    error stop 'FAIL: Real and DA Legendre values differ'
                call expected_da(row,col)%destroy()
            end do
        end do
        call zr_da%destroy()
        call pool%destroy()
    end subroutine check_degree
end program test_gravity_legendre_tail
