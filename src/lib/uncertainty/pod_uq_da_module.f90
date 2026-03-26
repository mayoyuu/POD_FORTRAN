!> @file pod_uq_da_module.f90
!> @brief 基于微分代数的轨道不确定性传播器实现
!> @author Song Yu
!> @date 2026-03-26
!> 这个模块实现了一个基于微分代数 (DA) 的轨道不确定性传播器，继承自 uq_propagator_base。
!> 核心思想是：通过一次 DA 积分得到中心轨迹的 DA 多项式展开，然后利用这个多项式快速映射所有粒子样本，大幅节省计算资源。适用于高维状态空间和大量粒子样本的
!> 不确定性传播场景。
!> 依赖于 pod_uq_base_module 定义的基类接口，以及 pod_uq_state_module 定义的状态数据结构。同时利用了底层的 DA 积分器（如 RKF45、RKF78）进行中心轨迹的 DA 积分。
!> 注意：这个实现假设输入状态的样本是状态的样本，因此需要额外声明中心样本，以便直接作为 DA 多项式的输入。如果输入状态不满足这个假设，可能需要在 propagate 内部进行一次预处理。
module pod_uq_da_module
    use pod_global, only: DP, MAX_STRING_LEN
    ! 引入基类和积分器常量
    use pod_uq_base_module, only: uq_propagator_base, INTEG_RK4, INTEG_RKF45, INTEG_RKF78
    use pod_uq_state_module, only: uq_state_type
    use pod_dace_classes
    ! 引入你的底层 DA 积分器和对应的方法常量
    use pod_da_integrator_module, only: da_adaptive_step_integrate, da_rk4_integrate, &
                                        METHOD_RKF45, METHOD_RKF78
    implicit none
    private
    public :: uq_da_propagator

    type, extends(uq_propagator_base) :: uq_da_propagator
        integer :: da_order = 2
    contains
        ! 绑定传播实现
        procedure, pass :: propagate => da_propagate
        ! 必须实现基类的 get_name_if 契约
        procedure, pass :: get_method_name => da_get_method_name
    end type uq_da_propagator

contains

    !> 获取方法名称
    function da_get_method_name(this) result(name)
        class(uq_da_propagator), intent(in) :: this
        character(len=MAX_STRING_LEN)       :: name
        name = "Differential Algebra (DA) Propagator"
    end function da_get_method_name

    !> 核心传播实现
    subroutine da_propagate(this, t_start, t_end, input_state, output_state)
        class(uq_da_propagator), intent(inout) :: this
        real(DP), intent(in)                   :: t_start, t_end
        type(uq_state_type), intent(in)        :: input_state
        type(uq_state_type), intent(inout)     :: output_state
        
        type(AlgebraicVector) :: state_da_0, state_da_f, state_da_next
        real(DP), allocatable :: times(:)
        type(AlgebraicVector), allocatable :: states(:)
        type(CompiledDA)      :: compiled_state
        integer :: n_steps, i, n_particles, dim
        real(DP) :: eval_inputs(6), eval_results(6)
        real(DP) :: t_current, dt_step

        dim = size(input_state%samples, 1)
        n_particles = size(input_state%samples, 2)
        
        ! 调用修正后的内存分配函数
        call output_state%allocate_memory(dim, n_particles)

        ! 1. 初始化 DA 中心状态 (以输入的均值为展开点)
        call dace_initialize(this%da_order, dim)
        call state_da_0%init(dim)
        do i = 1, dim
            state_da_0%elements(i) = input_state%mean(i) + da_var(i)
        end do

        ! 2. 核心：仅做一次 DA 动力学积分！(大幅节省算力)
        if (this%verbose) write(*,*) '[DA Propagator] 开始中心轨迹 DA 积分传播...'
        
        select case (this%integrator_type)
            case (INTEG_RKF45)
                call da_adaptive_step_integrate(state_da_0, t_start, t_end, 500000, &
                                                this%abs_error, METHOD_RKF45, times, states, n_steps)
                state_da_f = states(n_steps)
                
            case (INTEG_RKF78)
                call da_adaptive_step_integrate(state_da_0, t_start, t_end, 500000, &
                                                this%abs_error, METHOD_RKF78, times, states, n_steps)
                state_da_f = states(n_steps)
                
            case (INTEG_RK4)
                ! RK4 通常为定步长，这里加入一个简单的定步长循环到达 t_end
                state_da_f = state_da_0
                t_current = t_start
                do while (t_current < t_end)
                    ! 防止最后一步越界
                    dt_step = min(this%dt_initial, t_end - t_current) 
                    call da_rk4_integrate(state_da_f, dt_step, t_current, state_da_next)
                    state_da_f = state_da_next
                    t_current = t_current + dt_step
                end do
                
            case default
                write(*,*) '[ERROR] DA Propagator: 未知的积分器类型！'
                return
        end select

        ! 3. OpenMP 多线程粒子映射求值 (替代 TBB)
        ! ========================================================
        ! 3. DA 多项式预编译与多线程快速求值
        ! ========================================================
        if (this%verbose) write(*,*) '[DA Propagator] 积分完成，开始预编译多项式...'
        
        ! 【核心修复】：在 OpenMP 循环之外，单线程下编译，生成只读的求值计划
        compiled_state = state_da_f%compile()
        
        if (this%verbose) write(*,*) '[DA Propagator] 编译完成，开始多线程极速求值...'
        
        !$omp parallel do default(none) &
        !$omp private(i, eval_inputs, eval_results) &
        !$omp shared(n_particles, dim, input_state, output_state, compiled_state) ! <-- 注意共享 compiled_state
        do i = 1, n_particles
            ! 获取当前粒子相对中心点的偏差向量
            eval_inputs(:) = input_state%samples(:, i) - input_state%mean(:)
            
            ! 【核心修复】：直接调用编译好对象的 eval 方法！
            ! 它是纯粹的数学代入，绝对线程安全，且速度是之前的几百倍
            eval_results = compiled_state%eval(eval_inputs)
            
            ! 存储结果
            output_state%samples(:, i) = eval_results
        end do
        !$omp end parallel do
        
        ! 【核心修复】：循环结束后，统一释放一次内存
        call compiled_state%destroy()

        if (this%verbose) write(*,*) '[DA Propagator] 误差传播计算完毕。'

    end subroutine da_propagate
end module pod_uq_da_module