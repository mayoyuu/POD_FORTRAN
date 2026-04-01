module pod_uq_mc_module
    use pod_global, only: DP, MAX_STRING_LEN
    use pod_uq_base_module, only: uq_propagator_base, INTEG_RK4, INTEG_RKF45, INTEG_RKF78
    use pod_uq_state_module, only: uq_state_type
    use pod_integrator_module, only: rk4_integrate, adaptive_step_integrate, METHOD_RKF45, METHOD_RKF78
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
        real(DP), dimension(6) :: current_state, next_state
        real(DP) :: t_current, dt_step
        integer  :: step_idx, n_steps
        real(DP), allocatable :: temp_times(:)
        real(DP), allocatable :: temp_states(:,:)

        dim = size(input_state%samples, 1)
        n_particles = size(input_state%samples, 2)
        
        ! 1. 为输出分布分配内存
        call output_state%allocate_memory(dim, n_particles)

        if (this%verbose) write(*,*) '[MC Propagator] 开始多线程物理动力学积分...'
        if (this%verbose) write(*,*) '[MC Propagator] 并行粒子总数: ', n_particles

        ! 2. 开启 OpenMP 并行域
        ! 注意：default(none) 强制要求我们显式声明每一个变量的共享/私有属性，极其安全
        ! ... [分配 output_state 内存等常规操作] ...

        !$omp parallel do default(none) &
        !$omp private(i, current_state, next_state, t_current, dt_step, n_steps, step_idx, temp_times, temp_states) &
        !$omp shared(this, dim, n_particles, input_state, output_state, t_start, t_end)
        do i = 1, n_particles
            
            current_state = input_state%samples(:, i)
            
            select case (this%integrator_type)
                case (INTEG_RKF45)
                    call adaptive_step_integrate(current_state, t_start, t_end, METHOD_RKF45, &
                                                 temp_times, temp_states, n_steps)
                    current_state = temp_states(n_steps, :)
                    
                case (INTEG_RKF78)
                    call adaptive_step_integrate(current_state, t_start, t_end, METHOD_RKF78, &
                                                 temp_times, temp_states, n_steps)
                    current_state = temp_states(n_steps, :)
                    
                case (INTEG_RK4)
                    ! === RK4 定步长保留 ===
                    n_steps = int((t_end - t_start) / this%dt_initial)
                    dt_step = this%dt_initial
                    
                    do step_idx = 1, n_steps
                        t_current = t_start + real(step_idx - 1, DP) * dt_step
                        call rk4_integrate(current_state, dt_step, t_current, next_state)
                        current_state = next_state
                    end do
                    
                    ! 处理尾部残余步长
                    t_current = t_start + real(n_steps, DP) * dt_step
                    dt_step = t_end - t_current
                    if (dt_step > 1.0e-6_DP) then
                        call rk4_integrate(current_state, dt_step, t_current, next_state)
                        current_state = next_state
                    end if
            end select
            
            ! 存储该粒子的最终结果
            output_state%samples(:, i) = current_state
            
            ! 极其重要：在 OMP 循环末尾手动释放当前线程的内存，防止 10000 个粒子导致内存溢出
            if (allocated(temp_times)) deallocate(temp_times)
            if (allocated(temp_states)) deallocate(temp_states)
            
        end do
        !$omp end parallel do

        if (this%verbose) write(*,*) '[MC Propagator] 蒙特卡洛传播计算完毕。'

    end subroutine mc_propagate
end module pod_uq_mc_module