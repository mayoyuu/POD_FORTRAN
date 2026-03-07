!> # CAT DACE Core Module
module pod_dace_core
    use, intrinsic :: iso_c_binding
    implicit none
    private

    public :: DA

    ! 1. 绑定 C 层的基础 DA 操作
    interface
        subroutine c_daceAllocateDA(handle) bind(C, name="daceAllocateDA")
            import :: c_int
            integer(c_int), intent(out) :: handle
        end subroutine c_daceAllocateDA
        
        subroutine c_daceAdd(h1, h2, h_out) bind(C, name="daceAdd")
            import :: c_int
            integer(c_int), value :: h1, h2
            integer(c_int), value :: h_out
        end subroutine c_daceAdd
    end interface

    ! 2. 定义 Fortran 端的 DA 对象
    type :: DA
        integer(c_int) :: handle = -1 ! 默认无效句柄
    contains
        procedure :: init => da_init
        ! 可以添加析构函数 (Finalizer) 用于自动释放内存
        ! final :: da_destroy 
    end type DA

    ! 3. 重载操作符
    interface operator(+)
        module procedure da_add_da
    end interface
    public :: operator(+)

contains

    subroutine da_init(this)
        class(DA), intent(inout) :: this
        call c_daceAllocateDA(this%handle)
    end subroutine da_init

    ! 实现 DA + DA
    type(DA) function da_add_da(da1, da2) result(res)
        class(DA), intent(in) :: da1, da2
        call res%init() ! 分配结果的内存/句柄
        call c_daceAdd(da1%handle, da2%handle, res%handle)
    end function da_add_da

end module pod_dace_core