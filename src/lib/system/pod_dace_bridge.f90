!> # CAT DACE Bridge Module
!> 负责将 C++ 的 DACE 微分代数库通过 iso_c_binding 引入 Fortran
module pod_dace_bridge
    use, intrinsic :: iso_c_binding
    implicit none
    
    ! 必须设为 private，只暴露安全的 Fortran 封装接口给外部
    private
    public :: dace_initialize
    
    ! ---------------------------------------------------------
    ! 1. 裸露的 C 接口绑定 (严格对应 DACE C-API)
    ! ---------------------------------------------------------
    interface
        ! 对应 C 函数: void daceInit(int maxOrder, int maxVars);
        subroutine c_daceInit(max_order, max_vars) bind(C, name="daceInitialize")
            import :: c_int
            ! C++ 默认是传值传递 (pass by value)，Fortran 必须加 value 关键字！
            integer(c_int), value :: max_order
            integer(c_int), value :: max_vars
        end subroutine c_daceInit
    end interface
    
contains

    ! ---------------------------------------------------------
    ! 2. 高级 Fortran 安全封装层
    ! ---------------------------------------------------------
    !> 初始化 DACE 微分代数环境
    !> @param order DA 多项式的最高阶数 (例如 5 阶)
    !> @param vars  系统的独立变量个数 (例如 6 个状态量 = 6)
    subroutine dace_initialize(order, vars)
        integer, intent(in) :: order
        integer, intent(in) :: vars
        
        write(*, *) '[DACE Bridge] 正在向 C++ 核心发送初始化指令...'
        write(*, *) '[DACE Bridge] 设置最大阶数 (Order) = ', order
        write(*, *) '[DACE Bridge] 设置变量数 (Vars)  = ', vars
        
        ! 将 Fortran 的默认 integer 转换为 C 的 int 并调用
        call c_daceInit(int(order, c_int), int(vars, c_int))
        
        write(*, *) '[DACE Bridge] DACE 引擎初始化成功！'
    end subroutine dace_initialize

end module pod_dace_bridge