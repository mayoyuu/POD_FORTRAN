module pod_uq_crtbp_da_module
    use pod_global, only: DP
    use pod_crtbp_module, only: set_crtbp_mu
    use pod_crtbp_integrator_module, only: da_adaptive_integrate_crtbp
    use pod_dace_classes
    use pod_random_module, only: generate_multivariate_normal, init_random_seed
    implicit none
    private

    public :: crtbp_da_propagate

contains

    subroutine crtbp_da_propagate(nominal_state, cov, mu, t_end, &
            da_order, n_particles, rel_tol, abs_tol, dt_min, dt_max, max_steps, &
            final_samples, final_mean, final_cov, propagated_ref, verbose)
        real(DP), intent(in) :: nominal_state(6)
        real(DP), intent(in) :: cov(6,6)
        real(DP), intent(in) :: mu, t_end
        integer,  intent(in) :: da_order, n_particles
        real(DP), intent(in) :: rel_tol, abs_tol, dt_min, dt_max
        integer,  intent(in) :: max_steps
        real(DP), allocatable, intent(out) :: final_samples(:,:)
        real(DP), intent(out) :: final_mean(6)
        real(DP), intent(out) :: final_cov(6,6)
        real(DP), intent(out) :: propagated_ref(6)
        logical, intent(in) :: verbose

        type(AlgebraicVector) :: state_da_0, state_da_f
        type(CompiledDA) :: compiled
        integer :: i, j, dim
        real(DP) :: eval_inputs(6), eval_results(6)
        real(DP), allocatable :: particles(:,:)

        dim = 6

        ! 1. Generate particles (also used as deviation samples for DA eval)
        allocate(particles(dim, n_particles))
        call init_random_seed(.true.)
        call generate_multivariate_normal(nominal_state, cov, particles)

        ! 2. Set mu and push DA order
        call set_crtbp_mu(mu)
        call dace_push_to(da_order)

        ! 3. Build initial DA state: nominal + da_var(i)
        call state_da_0%init(dim)
        do i = 1, dim
            state_da_0%elements(i) = nominal_state(i) + da_var(i)
        end do

        ! 4. Single DA adaptive integration (returns final state only)
        call state_da_f%init(dim)
        if (verbose) write(*, '(A)') '[DA CRTBP] Starting DA integration...'

        call da_adaptive_integrate_crtbp(state_da_0, 0.0_DP, t_end, &
            rel_tol, abs_tol, dt_min, dt_max, max_steps, &
            state_da_f)

        ! 5. Extract constant part (propagated reference)
        propagated_ref = state_da_f%cons()

        ! 6. Compile polynomial
        compiled = state_da_f%compile()

        ! 7. Evaluate all particles
        allocate(final_samples(dim, n_particles))
        eval_inputs = 0.0_DP

        !$omp parallel do default(none) &
        !$omp private(i, eval_inputs, eval_results) &
        !$omp shared(n_particles, particles, nominal_state, compiled, final_samples, dim)
        do i = 1, n_particles
            eval_inputs(:) = particles(:, i) - nominal_state(:)
            eval_results = compiled%eval(eval_inputs)
            final_samples(:, i) = eval_results
        end do
        !$omp end parallel do

        ! 8. Compute moments
        final_mean = 0.0_DP
        do i = 1, n_particles
            final_mean = final_mean + final_samples(:, i)
        end do
        final_mean = final_mean / real(n_particles, DP)

        final_cov = 0.0_DP
        do i = 1, n_particles
            do j = 1, dim
                final_cov(:, j) = final_cov(:, j) + &
                    (final_samples(:, i) - final_mean) * (final_samples(j, i) - final_mean(j))
            end do
        end do
        final_cov = final_cov / real(n_particles - 1, DP)

        ! 9. Cleanup
        call compiled%destroy()
        call state_da_0%destroy()
        call state_da_f%destroy()
        deallocate(particles)
        call dace_pop_to()

        if (verbose) write(*, '(A)') '[DA CRTBP] Propagation complete.'
    end subroutine crtbp_da_propagate

end module pod_uq_crtbp_da_module
