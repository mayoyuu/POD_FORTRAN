module pod_uq_mc_module
    use pod_global, only: DP, MAX_STRING_LEN
    use pod_uq_base_module, only: uq_propagator_base
    use pod_uq_state_module, only: uq_state_type
    use pod_integrator_module, only: rk4_integrate, adaptive_step_integrate, METHOD_RKF45, METHOD_RKF78
    use pod_force_model_module, only: set_propagation_epoch
    use pod_config, only: config
    use pod_integrator_module
    implicit none
    private
    public :: uq_mc_propagator

    type, extends(uq_propagator_base) :: uq_mc_propagator
    contains
        procedure, pass :: propagate => mc_propagate
        procedure, pass :: get_method_name => mc_get_method_name
    end type uq_mc_propagator

contains

    !> 履行基类契约：返回方法名称
    function mc_get_method_name(this) result(name)
        class(uq_mc_propagator), intent(in) :: this
        character(len=MAX_STRING_LEN)       :: name
        if (this%integrator_type == -1) continue
        name = "Monte Carlo (MC) Propagator"
    end function mc_get_method_name

    !> 核心传播实现 (OpenMP 全量并行积分)
    subroutine mc_propagate(this, t_start, t_end, input_state, output_state)
        class(uq_mc_propagator), intent(inout) :: this
        real(DP), intent(in)                   :: t_start, t_end
        type(uq_state_type), intent(in)        :: input_state
        type(uq_state_type), intent(inout)     :: output_state
        
        integer :: n_particles, dim, i
        ! 必须声明在线程内部使用的局部变量 (会在 OpenMP 中设为 private)
        real(DP), dimension(6) :: current_state
        real(DP) :: t_start_nondim, t_end_nondim
        integer  :: n_steps
        real(DP), allocatable :: temp_times(:)
        real(DP), allocatable :: temp_states(:,:)


        dim = size(input_state%samples, 1)
        n_particles = size(input_state%samples, 2)
        
        call output_state%allocate_memory(dim, n_particles)

        ! ========================================================
        ! 核心修改 1：从 this%epoch0 获取基准历元，注入给 Real 物理引擎
        ! ========================================================
        call set_propagation_epoch(this%epoch0)

        ! 相对时间无量纲化
        t_start_nondim = t_start / config%TU
        t_end_nondim   = t_end / config%TU

        if (this%verbose) write(*,*) '[MC Propagator] 开始多线程物理动力学积分...'

        !$omp parallel do default(none) &
        !$omp private(i, current_state, temp_times, temp_states, n_steps) &
        !$omp shared(this, dim, n_particles, input_state, output_state, t_start_nondim, t_end_nondim, config)
        do i = 1, n_particles
            current_state = input_state%samples(:, i)
            
            ! 粒子独立无量纲化
            current_state(1:3) = current_state(1:3) / config%LU
            current_state(4:6) = current_state(4:6) / config%VU
            
            select case (this%integrator_type)
                case (METHOD_RKF78)
                    call adaptive_step_integrate(current_state, t_start_nondim, t_end_nondim, METHOD_RKF78, &
                                                 temp_times, temp_states, n_steps)
                    current_state = temp_states(n_steps, :)
                case (METHOD_RKF45)
                    call adaptive_step_integrate(current_state, t_start_nondim, t_end_nondim, METHOD_RKF45, &
                                                 temp_times, temp_states, n_steps)
                    current_state = temp_states(n_steps, :)
                case default
                    write(*,*) '[ERROR] DA Propagator: 未知的积分器类型！'
                    return
            end select
            
            ! 粒子独立还原量纲
            current_state(1:3) = current_state(1:3) * config%LU
            current_state(4:6) = current_state(4:6) * config%VU
            
            output_state%samples(:, i) = current_state
            
            if (allocated(temp_times)) deallocate(temp_times)
            if (allocated(temp_states)) deallocate(temp_states)
        end do
        !$omp end parallel do

        ! ========================================================
        ! 核心修改 2：调用封装好的 compute_moments 计算后验均值与协方差
        ! ========================================================
        call output_state%compute_moments()

        if (this%verbose) write(*,*) '[MC Propagator] 蒙特卡洛传播计算完毕。'
    end subroutine mc_propagate

end module pod_uq_mc_module