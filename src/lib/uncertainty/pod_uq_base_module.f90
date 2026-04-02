!> 不确定性量化/传播的抽象基类
module pod_uq_base_module
    use pod_global, only: DP, MAX_STRING_LEN
    use pod_uq_state_module, only: uq_state_type
    ! 引入底层积分器常量
    use pod_integrator_module, only: METHOD_RKF45, METHOD_RKF78
    
    implicit none
    private
    public :: uq_propagator_base

    type, abstract :: uq_propagator_base
        ! --- 基础配置参数 ---
        character(len=MAX_STRING_LEN) :: method_name = "Unknown" 
        logical  :: verbose = .true.          
        real(DP) :: epoch0 = 0.0_DP  ! 物理历元基准
        
        ! --- 积分器配置参数 ---
        ! 保留积分器类型，用于子类中的路由分发
        integer  :: integrator_type = METHOD_RKF78 
        
    contains
        ! 1. 核心的多态接口：强制子类必须实现
        procedure(propagate_if), deferred, pass :: propagate
        
        ! 2. 强制子类必须实现返回自己名字的方法
        procedure(get_name_if), deferred, pass :: get_method_name
        
        ! 3. 设置器方法
        procedure, pass :: set_integrator => base_set_integrator
        procedure, pass :: set_verbosity => base_set_verbosity
    end type uq_propagator_base

    ! ==========================================
    ! 接口定义区 (Deferred Interfaces)
    ! 【注意：就是这里！千万不能删】
    ! ==========================================
    abstract interface
        subroutine propagate_if(this, t_start, t_end, input_state, output_state)
            import :: uq_propagator_base, DP, uq_state_type
            class(uq_propagator_base), intent(inout) :: this
            real(DP), intent(in)                     :: t_start, t_end
            type(uq_state_type), intent(in)          :: input_state
            type(uq_state_type), intent(inout)       :: output_state
        end subroutine propagate_if
        
        function get_name_if(this) result(name)
            import :: uq_propagator_base, MAX_STRING_LEN
            class(uq_propagator_base), intent(in) :: this
            character(len=MAX_STRING_LEN)         :: name
        end function get_name_if
    end interface

contains

    ! ==========================================
    ! 基类具体实现区 (Base Implementations)
    ! ==========================================
    
    !> 极简的积分器设置
    subroutine base_set_integrator(this, integ_type)
        class(uq_propagator_base), intent(inout) :: this
        integer, intent(in) :: integ_type
        
        this%integrator_type = integ_type
    end subroutine base_set_integrator
    
    !> 设置是否静默输出
    subroutine base_set_verbosity(this, is_verbose)
        class(uq_propagator_base), intent(inout) :: this
        logical, intent(in) :: is_verbose
        this%verbose = is_verbose
    end subroutine base_set_verbosity

end module pod_uq_base_module