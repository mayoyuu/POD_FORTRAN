!> ADS uncertainty propagation module
!> Implements Automatic Domain Splitting based nonlinear uncertainty propagation.
!> Extends uq_propagator_base.
module pod_uq_prop_ads_module
    use pod_global, only: DP, MAX_STRING_LEN
    use pod_uq_base_module, only: uq_propagator_base
    use pod_uq_state_module, only: uq_state_type
    use pod_ads_split_module, only: manifold_type
    use pod_uq_crtbp_ads_module, only: crtbp_ads_propagate
    implicit none
    private

    public :: uq_ads_propagator

    type, extends(uq_propagator_base), public :: uq_ads_propagator
        integer  :: da_order = 2
        integer  :: n_split_max = 12
        real(DP) :: err_toll(6) = 1.0d-4
        real(DP) :: mu = 0.012153614091892_DP
        real(DP) :: ads_abs_tol = 1.0d-12
        real(DP) :: ads_rel_tol = 1.0d-12
        real(DP) :: ads_dt_min = 1.0d-6
        real(DP) :: ads_dt_max = 3600.0_DP
        integer  :: ads_max_steps = 100000
        type(manifold_type) :: domain
        integer  :: n_patches = 0
    contains
        procedure :: propagate       => ads_propagate
        procedure :: get_method_name => ads_get_method_name
        procedure :: set_ads_order
        procedure :: set_ads_tolerances
        procedure :: set_ads_mu
    end type uq_ads_propagator

contains

    function ads_get_method_name(this) result(name)
        class(uq_ads_propagator), intent(in) :: this
        character(len=MAX_STRING_LEN) :: name
        write(name, '(A,I0)') 'ADS-', this%da_order
    end function ads_get_method_name

    subroutine set_ads_order(this, order)
        class(uq_ads_propagator), intent(inout) :: this
        integer, intent(in) :: order
        this%da_order = order
    end subroutine set_ads_order

    subroutine set_ads_tolerances(this, abs_tol, rel_tol, dt_min, dt_max, &
                                   max_steps, err_toll, n_split_max)
        class(uq_ads_propagator), intent(inout) :: this
        real(DP), optional, intent(in) :: abs_tol, rel_tol, dt_min, dt_max
        integer, optional, intent(in) :: max_steps, n_split_max
        real(DP), optional, intent(in) :: err_toll(6)
        if (present(abs_tol)) this%ads_abs_tol = abs_tol
        if (present(rel_tol)) this%ads_rel_tol = rel_tol
        if (present(dt_min)) this%ads_dt_min = dt_min
        if (present(dt_max)) this%ads_dt_max = dt_max
        if (present(max_steps)) this%ads_max_steps = max_steps
        if (present(err_toll)) this%err_toll = err_toll
        if (present(n_split_max)) this%n_split_max = n_split_max
    end subroutine set_ads_tolerances

    subroutine set_ads_mu(this, mu)
        class(uq_ads_propagator), intent(inout) :: this
        real(DP), intent(in) :: mu
        this%mu = mu
    end subroutine set_ads_mu

    subroutine ads_propagate(this, t_start, t_end, input_state, output_state)
        class(uq_ads_propagator), intent(inout) :: this
        real(DP), intent(in) :: t_start, t_end
        type(uq_state_type), intent(in) :: input_state
        type(uq_state_type), intent(inout) :: output_state

        real(DP), allocatable :: final_samples(:,:)
        real(DP) :: final_mean(6), final_cov(6,6)
        integer :: n_particles, n_patches_out

        n_particles = 100000

        if (.not. allocated(input_state%mean)) then
            write(*,*) '[ADS] ERROR: input_state%mean not allocated'
            return
        end if

        call crtbp_ads_propagate( &
            input_state%mean(1:6), input_state%cov(1:6,1:6), &
            this%mu, t_end - t_start, this%da_order, &
            this%n_split_max, this%err_toll, &
            this%ads_rel_tol, this%ads_abs_tol, &
            this%ads_dt_min, this%ads_dt_max, this%ads_max_steps, &
            n_particles, &
            final_samples, final_mean, final_cov, n_patches_out, this%verbose)

        this%n_patches = n_patches_out

        call output_state%deallocate_memory()
        if (allocated(output_state%mean)) deallocate(output_state%mean)
        if (allocated(output_state%cov)) deallocate(output_state%cov)
        allocate(output_state%mean(6), output_state%cov(6,6))
        output_state%mean = final_mean
        output_state%cov = final_cov

        if (allocated(final_samples)) deallocate(final_samples)
    end subroutine ads_propagate

end module pod_uq_prop_ads_module
