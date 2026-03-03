module dace_bridge_module
    use, intrinsic :: iso_c_binding
    implicit none
    
    ! 声明与 C 语言指针等价的 Fortran 类型
    type, bind(c) :: c_da_ptr
        type(c_ptr) :: ptr
    end type c_da_ptr

    ! 绑定 C 函数 (严格对应 Step 1 中的 C 接口)
    interface
        subroutine c_daceInit(max_order, max_vars) bind(c, name="daceInit")
            import :: c_int
            integer(c_int), value :: max_order
            integer(c_int), value :: max_vars
        end subroutine c_daceInit

        function c_da_alloc() bind(c, name="da_alloc") result(res_ptr)
            import :: c_ptr
            type(c_ptr) :: res_ptr
        end function c_da_alloc

        subroutine c_da_free(da_ptr) bind(c, name="da_free")
            import :: c_ptr
            type(c_ptr), value :: da_ptr
        end subroutine c_da_free

        subroutine c_da_add(res, a, b) bind(c, name="da_add")
            import :: c_ptr
            type(c_ptr), value :: res, a, b
        end subroutine c_da_add
    end interface
end module dace_bridge_module