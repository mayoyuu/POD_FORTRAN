module pod_uq_crtbp_mc_module
    use pod_global, only: DP, MAX_STRING_LEN
    use pod_crtbp_module, only: set_crtbp_mu
    use pod_crtbp_integrator_module, only: adaptive_integrate_crtbp
    use pod_random_module, only: generate_multivariate_normal, init_random_seed
    implicit none
    private

    public :: crtbp_mc_propagate

contains

    subroutine crtbp_mc_propagate(nominal_state, cov, mu, t_end, &
            n_particles, rel_tol, abs_tol, dt_min, dt_max, max_steps, &
            final_samples, final_mean, final_cov, verbose)
        real(DP), intent(in) :: nominal_state(6)
        real(DP), intent(in) :: cov(6,6)
        real(DP), intent(in) :: mu, t_end
        integer,  intent(in) :: n_particles
        real(DP), intent(in) :: rel_tol, abs_tol, dt_min, dt_max
        integer,  intent(in) :: max_steps
        real(DP), allocatable, intent(out) :: final_samples(:,:)
        real(DP), intent(out) :: final_mean(6)
        real(DP), intent(out) :: final_cov(6,6)
        logical, intent(in) :: verbose

        integer :: dim, i, n_steps, j
        real(DP), allocatable :: temp_times(:), temp_states(:,:)
        real(DP), dimension(6) :: current_state

        dim = 6

        allocate(final_samples(dim, n_particles))

        call set_crtbp_mu(mu)
        call init_random_seed(.true.)
        call generate_multivariate_normal(nominal_state, cov, final_samples)

        if (verbose) write(*, '(A,I0)') '[MC CRTBP] Propagating ', n_particles

        !$omp parallel do default(none) &
        !$omp private(i, current_state, temp_times, temp_states, n_steps) &
        !$omp shared(n_particles, final_samples, t_end, rel_tol, abs_tol, dt_min, dt_max, max_steps)
        do i = 1, n_particles
            current_state = final_samples(:, i)

            call adaptive_integrate_crtbp(current_state, 0.0_DP, t_end, &
                rel_tol, abs_tol, dt_min, dt_max, max_steps, &
                temp_times, temp_states, n_steps)

            final_samples(:, i) = temp_states(n_steps, :)

            if (allocated(temp_times)) deallocate(temp_times)
            if (allocated(temp_states)) deallocate(temp_states)
        end do
        !$omp end parallel do

        ! Compute mean
        final_mean = 0.0_DP
        do i = 1, n_particles
            final_mean = final_mean + final_samples(:, i)
        end do
        final_mean = final_mean / real(n_particles, DP)

        ! Compute covariance
        final_cov = 0.0_DP
        do i = 1, n_particles
            do j = 1, dim
                final_cov(:, j) = final_cov(:, j) + &
                    (final_samples(:, i) - final_mean) * (final_samples(j, i) - final_mean(j))
            end do
        end do
        final_cov = final_cov / real(n_particles - 1, DP)

        if (verbose) write(*, '(A)') '[MC CRTBP] Propagation complete.'
    end subroutine crtbp_mc_propagate

end module pod_uq_crtbp_mc_module
