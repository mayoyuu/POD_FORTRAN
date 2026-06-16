!> CRTBP ADS uncertainty propagation module
!> Provides CRTBP-specific ADS propagation with explicit DA integrator calls.
module pod_uq_crtbp_ads_module
    use pod_global, only: DP
    use pod_ads_split_module, only: splitting_history_type, patch_type, &
        manifold_type, sh_push, sh_pop, sh_count, sh_center, sh_width, &
        sh_contain, sh_map_point, patch_init, patch_destroy, &
        patch_get_trunc_err, patch_get_split_dir, patch_split, &
        mf_init, mf_destroy, mf_push, mf_pop_front
    use pod_crtbp_integrator_module, only: da_adaptive_integrate_crtbp
    use pod_crtbp_module, only: set_crtbp_mu
    use pod_dace_classes, only: AlgebraicVector, DA, da_var, &
        dace_initialize, dace_push_to, dace_pop_to, &
        operator(+), operator(*), assignment(=)
    use pod_random_module, only: generate_multivariate_normal, init_random_seed
    implicit none
    private

    public :: crtbp_ads_propagate, crtbp_ads_propagate_deviates

contains

    ! =========================================================================
    ! ads_get_split_domain_crtbp -- CRTBP-specific BFS domain splitting
    !   Calls da_adaptive_integrate_crtbp for flow evaluation.
    !   Corresponds to C++ Manifold::getSplitDomain
    ! =========================================================================
    subroutine ads_get_split_domain_crtbp(initial_da_vec, err_toll, n_split_max, &
                                           mu, t_end, abs_tol, rel_tol, dt_min, dt_max, &
                                           max_steps, da_order, result_manifold, verbose)
        type(AlgebraicVector), intent(in) :: initial_da_vec
        real(DP), intent(in) :: err_toll(6)
        integer, intent(in) :: n_split_max
        real(DP), intent(in) :: mu, t_end
        real(DP), intent(in) :: abs_tol, rel_tol, dt_min, dt_max
        integer, intent(in) :: max_steps, da_order
        type(manifold_type), intent(out) :: result_manifold
        logical, intent(in) :: verbose

        type(manifold_type) :: queue
        type(patch_type) :: p, f, left, right
        type(AlgebraicVector) :: f_da_vec, initial_copy
        real(DP) :: errors(6), rel_errors(6)
        integer :: pos(1), dir, iter_count, max_queue_size, i

        call mf_init(result_manifold)
        call mf_init(queue)

        ! Deep-copy initial_da_vec (intent(in)) into a local then move into patch
        call initial_copy%init(6)
        do i = 1, 6
            initial_copy%elements(i) = initial_da_vec%elements(i)
        end do
        call patch_init(p, initial_copy)
        call initial_copy%destroy()
        call mf_push(queue, p)

        iter_count = 0
        max_queue_size = 0

        ! Set CRTBP mu before any DA integrator calls
        call set_crtbp_mu(mu)

        do while (queue%n_patches > 0)
            iter_count = iter_count + 1
            max_queue_size = max(max_queue_size, queue%n_patches)

            call mf_pop_front(queue, p)

            ! --- CRTBP DA integration ---
            call f_da_vec%init(6)
            call da_adaptive_integrate_crtbp( &
                p%da_vec, 0.0_DP, t_end, &
                rel_tol, abs_tol, dt_min, dt_max, max_steps, &
                f_da_vec)

            call patch_init(f, f_da_vec, p%history)
            call f_da_vec%destroy()

            ! --- Truncation error ---
            call patch_get_trunc_err(f, da_order, errors)
            rel_errors = max(0.0_DP, errors - err_toll)

            if (maxval(rel_errors) <= 0.0_DP .or. &
                sh_count(p%history, 0) >= n_split_max) then
                call mf_push(result_manifold, f)
                call patch_destroy(p)
            else
                pos = maxloc(rel_errors)
                dir = patch_get_split_dir(f, pos(1), da_order)
                call patch_split(p, dir, left, right)
                call mf_push(queue, left)
                call mf_push(queue, right)
                call patch_destroy(p)
                call patch_destroy(f)
            end if
        end do

        call mf_destroy(queue)

        if (verbose) then
            write(*,'(A,I0,A,I0,A,I0)') '[ADS CRTBP] splitting done: ', &
                result_manifold%n_patches, ' patches, ', &
                iter_count, ' iterations, max queue size ', max_queue_size
        end if
    end subroutine ads_get_split_domain_crtbp

    ! =========================================================================
    ! crtbp_ads_propagate_deviates -- propagate given deviates via ADS
    ! =========================================================================
    subroutine crtbp_ads_propagate_deviates(nominal_state, deviates, mu, t_end, &
            da_order, n_split_max, err_toll, rel_tol, abs_tol, dt_min, dt_max, max_steps, &
            final_samples, final_mean, final_cov, n_patches_out, verbose)
        real(DP), intent(in) :: nominal_state(6)
        real(DP), intent(in) :: deviates(:,:)
        real(DP), intent(in) :: mu, t_end
        integer,  intent(in) :: da_order, n_split_max, max_steps
        real(DP), intent(in) :: err_toll(6)
        real(DP), intent(in) :: rel_tol, abs_tol, dt_min, dt_max
        real(DP), allocatable, intent(out) :: final_samples(:,:)
        real(DP), intent(out) :: final_mean(6), final_cov(6,6)
        integer,  intent(out) :: n_patches_out
        logical, intent(in) :: verbose

        type(AlgebraicVector) :: state_da_0
        type(manifold_type) :: manifold
        real(DP) :: pt_local(6)
        real(DP), allocatable :: eval_res(:)
        integer :: dim, n_dev, i, j, k
        logical :: found

        dim = 6
        n_dev = size(deviates, 2)

        call dace_initialize(da_order, dim)

        call state_da_0%init(dim)
        do i = 1, dim
            state_da_0%elements(i) = nominal_state(i) + da_var(i)
        end do

        if (verbose) write(*,'(A)') '[ADS CRTBP] Starting domain splitting...'
        call ads_get_split_domain_crtbp(state_da_0, err_toll, n_split_max, &
            mu, t_end, abs_tol, rel_tol, dt_min, dt_max, max_steps, &
            da_order, manifold, verbose)
        n_patches_out = manifold%n_patches

        allocate(final_samples(dim, n_dev))
        final_samples = 0.0_DP

        do k = 1, n_dev
            found = .false.
            do i = 1, manifold%n_patches
                if (sh_contain(manifold%patches(i)%history, deviates(:, k))) then
                    pt_local = deviates(:, k)
                    call sh_map_point(manifold%patches(i)%history, pt_local)
                    eval_res = manifold%patches(i)%da_vec%eval(pt_local)
                    final_samples(:, k) = eval_res
                    found = .true.
                    deallocate(eval_res)
                    exit
                end if
            end do
            if (.not. found) then
                write(*,'(A,I0,A)') '[ADS CRTBP] WARNING: deviate ', k, ' not in any patch!'
                final_samples(:, k) = nominal_state
            end if
        end do

        final_mean = 0.0_DP
        do k = 1, n_dev
            final_mean = final_mean + final_samples(:, k)
        end do
        final_mean = final_mean / real(n_dev, DP)

        final_cov = 0.0_DP
        do k = 1, n_dev
            do j = 1, dim
                final_cov(:, j) = final_cov(:, j) + &
                    (final_samples(:, k) - final_mean) * (final_samples(j, k) - final_mean(j))
            end do
        end do
        final_cov = final_cov / real(n_dev - 1, DP)

        call state_da_0%destroy()
        call mf_destroy(manifold)

        if (verbose) write(*,'(A)') '[ADS CRTBP] Propagation complete.'
    end subroutine crtbp_ads_propagate_deviates

    ! =========================================================================
    ! crtbp_ads_propagate -- full cov => samples => propagate pipeline
    ! =========================================================================
    subroutine crtbp_ads_propagate(nominal_state, cov, mu, t_end, &
            da_order, n_split_max, err_toll, rel_tol, abs_tol, dt_min, dt_max, max_steps, &
            n_particles, final_samples, final_mean, final_cov, n_patches_out, verbose)
        real(DP), intent(in) :: nominal_state(6), cov(6,6), mu, t_end
        integer,  intent(in) :: da_order, n_split_max, max_steps, n_particles
        real(DP), intent(in) :: err_toll(6)
        real(DP), intent(in) :: rel_tol, abs_tol, dt_min, dt_max
        real(DP), allocatable, intent(out) :: final_samples(:,:)
        real(DP), intent(out) :: final_mean(6), final_cov(6,6)
        integer,  intent(out) :: n_patches_out
        logical, intent(in) :: verbose

        real(DP), allocatable :: particles(:,:), deviates(:,:)
        integer :: i, dim

        dim = 6

        allocate(particles(dim, n_particles))
        call init_random_seed(.true.)
        call generate_multivariate_normal(nominal_state, cov, particles)

        allocate(deviates(dim, n_particles))
        do i = 1, n_particles
            deviates(:, i) = particles(:, i) - nominal_state
        end do

        call crtbp_ads_propagate_deviates(nominal_state, deviates, mu, t_end, &
            da_order, n_split_max, err_toll, rel_tol, abs_tol, dt_min, dt_max, max_steps, &
            final_samples, final_mean, final_cov, n_patches_out, verbose)

        deallocate(particles, deviates)
    end subroutine crtbp_ads_propagate

end module pod_uq_crtbp_ads_module
