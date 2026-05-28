!> @file pod_uq_prop_ut_module.f90
!> @brief 基于无迹变换 (Unscented Transform) 的轨道不确定性传播器
!> @author Song Yu
!> @date 2026-05-26
!> 继承自 uq_propagator_base，复用 pod_uq_transform_ut_module 的 sigma 点生成
!> 和矩重建函数。通过 2n+1 个 sigma 点的确定性传播，以远少于 MC 的积分次数
!> 获得传播后的均值和协方差。
module pod_uq_ut_module
    use pod_global, only: DP, MAX_STRING_LEN
    use pod_uq_base_module, only: uq_propagator_base
    use pod_uq_state_module, only: uq_state_type
    use pod_uq_transform_ut_module, only: generate_sigma_points, reconstruct_ut_moments
    use pod_integrator_module, only: adaptive_step_integrate, METHOD_RKF45, METHOD_RKF78
    use pod_force_model_module, only: set_propagation_epoch
    use pod_config, only: config
    implicit none
    private
    public :: uq_ut_propagator

    type, extends(uq_propagator_base) :: uq_ut_propagator
        real(DP) :: kappa = 0.5_DP
    contains
        procedure, pass :: propagate => ut_propagate
        procedure, pass :: get_method_name => ut_get_method_name
        procedure, pass :: set_kappa => ut_set_kappa
    end type uq_ut_propagator

contains

    function ut_get_method_name(this) result(name)
        class(uq_ut_propagator), intent(in) :: this
        character(len=MAX_STRING_LEN)       :: name
        if (this%integrator_type == -1) continue
        name = "Unscented Transform (UT) Propagator"
    end function ut_get_method_name

    subroutine ut_set_kappa(this, kappa_val)
        class(uq_ut_propagator), intent(inout) :: this
        real(DP), intent(in) :: kappa_val
        this%kappa = kappa_val
    end subroutine ut_set_kappa

    subroutine ut_propagate(this, t_start, t_end, input_state, output_state)
        class(uq_ut_propagator), intent(inout) :: this
        real(DP), intent(in)                   :: t_start, t_end
        type(uq_state_type), intent(in)        :: input_state
        type(uq_state_type), intent(inout)     :: output_state

        integer :: dim, n_sigma, i
        real(DP), allocatable :: sigmas_x(:,:), propagated_sigmas(:,:), weights(:)
        real(DP), allocatable :: mean_y(:), P_yy(:,:), P_xy(:,:)
        real(DP), dimension(6) :: current_state
        real(DP) :: t_start_nondim, t_end_nondim
        integer  :: n_steps
        real(DP), allocatable :: temp_times(:), temp_states(:,:)

        dim = size(input_state%mean)
        n_sigma = 2 * dim + 1

        call generate_sigma_points(input_state%mean, input_state%cov, this%kappa, &
                                    sigmas_x, weights)

        allocate(propagated_sigmas(dim, n_sigma))

        call set_propagation_epoch(this%epoch0)

        t_start_nondim = t_start / config%TU
        t_end_nondim   = t_end / config%TU

        if (this%verbose) write(*,*) '[UT Propagator] 开始传播 Sigma 点, 共 ', n_sigma, ' 个...'

        !$omp parallel do default(none) &
        !$omp private(i, current_state, temp_times, temp_states, n_steps) &
        !$omp shared(dim, n_sigma, sigmas_x, propagated_sigmas, t_start_nondim, t_end_nondim, this, config)
        do i = 1, n_sigma
            current_state = sigmas_x(:, i)

            current_state(1:3) = current_state(1:3) / config%LU
            current_state(4:6) = current_state(4:6) / config%VU

            select case (this%integrator_type)
                case (METHOD_RKF78)
                    call adaptive_step_integrate(current_state, t_start_nondim, t_end_nondim, &
                                                 METHOD_RKF78, temp_times, temp_states, n_steps)
                    current_state = temp_states(n_steps, :)
                case (METHOD_RKF45)
                    call adaptive_step_integrate(current_state, t_start_nondim, t_end_nondim, &
                                                 METHOD_RKF45, temp_times, temp_states, n_steps)
                    current_state = temp_states(n_steps, :)
                case default
                    write(*,*) '[ERROR] UT Propagator: 未知的积分器类型！'
                    return
            end select

            current_state(1:3) = current_state(1:3) * config%LU
            current_state(4:6) = current_state(4:6) * config%VU

            propagated_sigmas(:, i) = current_state

            if (allocated(temp_times)) deallocate(temp_times)
            if (allocated(temp_states)) deallocate(temp_states)
        end do
        !$omp end parallel do

        call reconstruct_ut_moments(sigmas_x, propagated_sigmas, weights, mean_y, P_yy, P_xy)

        if (allocated(output_state%mean)) deallocate(output_state%mean)
        allocate(output_state%mean(dim))
        output_state%mean = mean_y

        if (allocated(output_state%cov)) deallocate(output_state%cov)
        allocate(output_state%cov(dim, dim))
        output_state%cov = P_yy

        deallocate(sigmas_x, propagated_sigmas, weights)
        deallocate(mean_y, P_yy, P_xy)

        if (this%verbose) write(*,*) '[UT Propagator] UT 传播计算完毕。'
    end subroutine ut_propagate

end module pod_uq_ut_module
