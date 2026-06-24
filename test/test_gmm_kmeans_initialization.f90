program test_gmm_kmeans_initialization
    use pod_global, only: DP
    use pod_gmm_math_module, only: fit_gmm_to_particles
    use pod_uq_gmm_state_module, only: uq_gmm_state_type
    implicit none

    integer, parameter :: n_particles = 34
    real(DP) :: particles(1, n_particles)
    real(DP) :: ordered(1, n_particles)
    type(uq_gmm_state_type) :: gmm
    real(DP) :: centers(3)
    integer :: i

    do i = 1, 30
        ordered(1, i) = -1.0_DP + 2.0_DP * real(i - 1, DP) / 29.0_DP
    end do
    ordered(1, 31:34) = [5.0_DP, 5.2_DP, 50.0_DP, 50.2_DP]

    do i = 1, n_particles
        particles(:, i) = ordered(:, n_particles + 1 - i)
    end do

    call gmm%allocate_components(3, 1)
    call fit_gmm_to_particles(particles, gmm, 0, 1.0e-8_DP)

    do i = 1, 3
        centers(i) = gmm%components(i)%mean(1)
    end do
    call sort3(centers)

    if (abs(centers(1)) > 1.0_DP) error stop "missing large cluster near 0"
    if (abs(centers(2) - 5.1_DP) > 0.5_DP) error stop "missing middle cluster near 5.1"
    if (abs(centers(3) - 50.1_DP) > 0.5_DP) error stop "missing far cluster near 50.1"

    write(*,*) "test_gmm_kmeans_initialization passed"

contains
    subroutine sort3(vals)
        real(DP), intent(inout) :: vals(3)
        real(DP) :: tmp
        integer :: a, b

        do a = 1, 2
            do b = a + 1, 3
                if (vals(b) < vals(a)) then
                    tmp = vals(a)
                    vals(a) = vals(b)
                    vals(b) = tmp
                end if
            end do
        end do
    end subroutine sort3
end program test_gmm_kmeans_initialization
