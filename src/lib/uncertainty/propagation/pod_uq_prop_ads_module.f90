!> Object-oriented adapter for HFEM Automatic Domain Splitting.
module pod_uq_ads_module
    use pod_global, only: DP, MAX_STRING_LEN
    use pod_uq_base_module, only: uq_propagator_base
    use pod_uq_state_module, only: uq_state_type
    use pod_uq_hfem_ads_module, only: hfem_ads_options_type, &
        hfem_ads_stats_type, hfem_ads_propagate
    implicit none
    private

    type, extends(uq_propagator_base), public :: uq_ads_propagator
        type(hfem_ads_options_type) :: options
        type(hfem_ads_stats_type) :: last_stats
        integer :: last_status = 0
        character(len=MAX_STRING_LEN) :: last_message = ''
    contains
        procedure :: propagate => ads_propagate
        procedure :: get_method_name => ads_get_method_name
        procedure :: set_ads_order
        procedure :: set_ads_options
    end type uq_ads_propagator

contains

    function ads_get_method_name(this) result(name)
        class(uq_ads_propagator), intent(in) :: this
        character(len=MAX_STRING_LEN) :: name

        write(name,'(A,I0)') 'HFEM-ADS-', this%options%da_order
    end function ads_get_method_name

    subroutine set_ads_order(this, order)
        class(uq_ads_propagator), intent(inout) :: this
        integer, intent(in) :: order

        this%options%da_order = order
    end subroutine set_ads_order

    subroutine set_ads_options(this, options)
        class(uq_ads_propagator), intent(inout) :: this
        type(hfem_ads_options_type), intent(in) :: options

        this%options = options
    end subroutine set_ads_options

    subroutine ads_propagate(this, t_start, t_end, input_state, output_state)
        class(uq_ads_propagator), intent(inout) :: this
        real(DP), intent(in) :: t_start, t_end
        type(uq_state_type), intent(in) :: input_state
        type(uq_state_type), intent(inout) :: output_state
        real(DP), allocatable :: propagated_samples(:,:)

        this%last_status = 0
        this%last_message = ''
        call output_state%deallocate_memory()
        if (.not. allocated(input_state%mean) .or. &
            .not. allocated(input_state%cov) .or. &
            .not. allocated(input_state%samples)) then
            this%last_status = -1
            this%last_message = 'ADS input state requires mean, covariance, and samples'
            return
        end if
        if (size(input_state%mean) < 6 .or. size(input_state%cov,1) < 6 .or. &
            size(input_state%cov,2) < 6) then
            this%last_status = -2
            this%last_message = 'ADS input state has fewer than six orbit components'
            return
        end if

        call hfem_ads_propagate(input_state%mean(1:6), &
            input_state%cov(1:6,1:6), this%epoch0, t_start, t_end, &
            input_state%samples, this%options, propagated_samples, &
            this%last_stats, this%last_status, this%last_message)
        if (this%last_status /= 0) return
        call output_state%allocate_memory(size(propagated_samples,1), &
            size(propagated_samples,2))
        output_state%samples = propagated_samples
        call output_state%compute_moments()
    end subroutine ads_propagate

end module pod_uq_ads_module
