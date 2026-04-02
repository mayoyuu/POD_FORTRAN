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
    use pod_uq_base_module, only: uq_propagator_base
    use pod_uq_state_module, only: uq_state_type
    use pod_dace_classes
    use pod_da_force_model_module, only: set_propagation_epoch
    use pod_config, only: config
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
        if (this%integrator_type == -1) continue
        name = "Differential Algebra (DA) Propagator"
    end function da_get_method_name

    !> 核心传播实现

    !> 核心传播实现
    subroutine da_propagate(this, t_start, t_end, input_state, output_state)
        class(uq_da_propagator), intent(inout) :: this
        real(DP), intent(in)                   :: t_start, t_end
        type(uq_state_type), intent(in)        :: input_state
        type(uq_state_type), intent(inout)     :: output_state
        
        type(AlgebraicVector) :: state_da_0, state_da_f
        real(DP), allocatable :: times(:)
        type(AlgebraicVector), allocatable :: states(:)
        type(CompiledDA)      :: compiled_state
        integer :: n_steps, i, n_particles, dim
        real(DP) :: eval_inputs(6), eval_results(6)
        real(DP) :: t_start_nondim, t_end_nondim
    

        dim = size(input_state%samples, 1)
        n_particles = size(input_state%samples, 2)
        
        ! 1. 为输出分布分配内存
        call output_state%allocate_memory(dim, n_particles)

        ! ========================================================
        ! 核心修改 1：从 this%epoch0 获取基准历元，注入给 DA 物理引擎
        ! ========================================================
        call set_propagation_epoch(this%epoch0)

        ! 2. 初始化 DA 中心状态并注入独立方差
        call dace_initialize(this%da_order, dim)
        call state_da_0%init(dim)
        do i = 1, 3
            ! 物理均值 + 物理摄动量，随后直接除以特征量度进行无量纲化
            state_da_0%elements(i)   = (input_state%mean(i) + da_var(i)) / config%LU
            state_da_0%elements(i+3) = (input_state%mean(i+3) + da_var(i+3)) / config%VU
        end do

        ! 将相对传播时间无量纲化
        t_start_nondim = t_start / config%TU
        t_end_nondim   = t_end / config%TU

        ! 3. 核心：仅做一次 DA 动力学积分
        if (this%verbose) write(*,*) '[DA Propagator] 开始中心轨迹 DA 积分传播...'
        select case (this%integrator_type)
            case (METHOD_RKF78)
                call da_adaptive_step_integrate(state_da_0, t_start_nondim, t_end_nondim, &
                                                METHOD_RKF78, times, states, n_steps)
                state_da_f = states(n_steps)
            case (METHOD_RKF45)
                call da_adaptive_step_integrate(state_da_0, t_start_nondim, t_end_nondim, &
                                                METHOD_RKF45, times, states, n_steps)
                state_da_f = states(n_steps)
            case default
                write(*,*) '[ERROR] DA Propagator: 未知的积分器类型！'
                return
            ! ... [保留其他的 case，例如 RKF45] ...
        end select

        ! 4. 将积分得到的多项式还原回物理量纲
        do i = 1, 3
            state_da_f%elements(i)   = state_da_f%elements(i) * config%LU
            state_da_f%elements(i+3) = state_da_f%elements(i+3) * config%VU
        end do

        ! 5. 编译与多线程求值 (天然免疫量纲与历元，纯数学操作)
        compiled_state = state_da_f%compile()
        
        !$omp parallel do default(none) &
        !$omp private(i, eval_inputs, eval_results) &
        !$omp shared(n_particles, dim, input_state, output_state, compiled_state)
        do i = 1, n_particles
            ! 偏差向量：当前粒子 - 均值
            eval_inputs(:) = input_state%samples(:, i) - input_state%mean(:)
            eval_results = compiled_state%eval(eval_inputs)
            output_state%samples(:, i) = eval_results
        end do
        !$omp end parallel do
        
        call compiled_state%destroy()

        ! ========================================================
        ! 核心修改 2：调用封装好的 compute_moments 计算后验均值与协方差
        ! ========================================================
        call output_state%compute_moments()

        if (this%verbose) write(*,*) '[DA Propagator] 误差传播计算完毕。'
    end subroutine da_propagate
end module pod_uq_da_module