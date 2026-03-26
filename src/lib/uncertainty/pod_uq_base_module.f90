!> 不确定性量化/传播的抽象基类
module pod_uq_base_module
    use pod_global, only: DP, MAX_STRING_LEN
    use pod_uq_state_module, only: uq_state_type
    implicit none
    private
    public :: uq_propagator_base

    ! 定义积分器类型的常量（可以移到全局配置模块中）
    integer, parameter, public :: INTEG_RK4   = 1
    integer, parameter, public :: INTEG_RKF45 = 2
    integer, parameter, public :: INTEG_RKF78 = 3

    type, abstract :: uq_propagator_base
        ! --- 基础配置参数 ---
        character(len=MAX_STRING_LEN) :: method_name = "Unknown" ! 方法名称 (DA, MC, UT)
        logical  :: verbose = .true.          ! 是否在控制台打印进度
        
        ! --- 积分器配置参数 ---
        integer  :: integrator_type = INTEG_RKF78 ! 默认使用 RKF78
        real(DP) :: dt_initial = 60.0_DP      ! 初始步长
        real(DP) :: abs_error = 1.0D-12       ! 绝对容差
        real(DP) :: rel_error = 1.0D-12       ! 相对容差
        logical  :: use_adaptive_ctrl = .true.! 是否变步长
        
    contains
        ! 1. 核心的多态接口：强制子类必须实现
        procedure(propagate_if), deferred, pass :: propagate
        
        ! 2. 强制子类必须实现返回自己名字的方法
        procedure(get_name_if), deferred, pass :: get_method_name
        
        ! 3. 通用的 Setter 方法（基类实现，所有子类继承可用）
        procedure, pass :: set_integration_params => base_set_integration_params
        procedure, pass :: set_verbosity => base_set_verbosity
    end type uq_propagator_base

    ! ==========================================
    ! 接口定义区 (Deferred Interfaces)
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
    
    !> 统一设置积分参数，包含简单的安全性校验
    subroutine base_set_integration_params(this, integ_type, dt_init, abs_tol, rel_tol, adaptive)
        class(uq_propagator_base), intent(inout) :: this
        integer, intent(in), optional  :: integ_type
        real(DP), intent(in), optional :: dt_init, abs_tol, rel_tol
        logical, intent(in), optional  :: adaptive
        
        if (present(integ_type)) this%integrator_type = integ_type
        if (present(dt_init))    this%dt_initial = dt_init
        if (present(adaptive))   this%use_adaptive_ctrl = adaptive
        
        if (present(abs_tol)) then
            if (abs_tol <= 0.0_DP) write(*,*) '[Warn] UP Base: abs_tol must be > 0. Ignored.'
            if (abs_tol > 0.0_DP) this%abs_error = abs_tol
        end if
        
        if (present(rel_tol)) then
            if (rel_tol <= 0.0_DP) write(*,*) '[Warn] UP Base: rel_tol must be > 0. Ignored.'
            if (rel_tol > 0.0_DP) this%rel_error = rel_tol
        end if
    end subroutine base_set_integration_params
    
    !> 设置是否静默输出
    subroutine base_set_verbosity(this, is_verbose)
        class(uq_propagator_base), intent(inout) :: this
        logical, intent(in) :: is_verbose
        this%verbose = is_verbose
    end subroutine base_set_verbosity

end module pod_uq_base_module