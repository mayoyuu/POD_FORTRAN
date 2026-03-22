!> # CAT DACE Classes Module
module pod_dace_classes
    use, intrinsic :: iso_c_binding
    implicit none
    private
    
    public :: dace_initilize
    public :: DA, AlgebraicVector
    public :: operator(+), operator(-), operator(*), operator(/)
    public :: assignment(=)
    public :: sin, cos, exp, sqrt

    ! =========================================================
    ! 1. C 接口绑定
    ! =========================================================
    interface
        ! [新增] 引擎全局初始化
        subroutine c_daceInit(max_order, max_vars) bind(C, name="fdace_core_init")
            import :: c_int
            integer(c_int), value :: max_order
            integer(c_int), value :: max_vars
        end subroutine c_daceInit
        ! 基础接口：分配、释放、复制 DA 句柄
        subroutine c_fdace_allocate(h) bind(C, name="fdace_allocate")
            import :: c_int; integer(c_int), intent(out) :: h
        end subroutine c_fdace_allocate

        subroutine c_fdace_allocate_var(h, var_idx) bind(C, name="fdace_allocate_var")
            import :: c_int; integer(c_int), intent(out) :: h; integer(c_int), value :: var_idx
        end subroutine c_fdace_allocate_var

        subroutine c_fdace_free(h) bind(C, name="fdace_free")
            import :: c_int; integer(c_int), value :: h
        end subroutine c_fdace_free

        subroutine c_fdace_copy(h_src, h_dest) bind(C, name="fdace_copy")
            import :: c_int; integer(c_int), value :: h_src, h_dest
        end subroutine c_fdace_copy

        ! 编译接口：接受多个 DA 句柄，返回一个新的编译后 DA 句柄 泛化的需求
        subroutine c_fdace_compile_vector(da_handles, size, cda_handle) bind(C, name="fdace_compile_vector")
            import :: c_int
            integer(c_int), intent(in) :: da_handles(*)
            integer(c_int), value :: size
            integer(c_int), intent(out) :: cda_handle
        end subroutine c_fdace_compile_vector

        subroutine c_fdace_compiled_free(h) bind(C, name="fdace_compiled_free")
            import :: c_int
            integer(c_int), value :: h
        end subroutine c_fdace_compiled_free

        subroutine c_fdace_compiled_eval_double(cda_handle, in_args, n_args, out_res, n_res) bind(C, name="fdace_compiled_eval_double")
            import :: c_int, c_double
            integer(c_int), value :: cda_handle, n_args, n_res
            real(c_double), intent(in) :: in_args(*)
            real(c_double), intent(out) :: out_res(*)
        end subroutine c_fdace_compiled_eval_double

        ! 基础代数运算
        subroutine c_fdace_add(h1, h2, ho) bind(C, name="fdace_add")
            import :: c_int; integer(c_int), value :: h1, h2, ho
        end subroutine c_fdace_add

        subroutine c_fdace_sub(h1, h2, ho) bind(C, name="fdace_sub")
            import :: c_int; integer(c_int), value :: h1, h2, ho
        end subroutine c_fdace_sub

        subroutine c_fdace_mul(h1, h2, ho) bind(C, name="fdace_mul")
            import :: c_int; integer(c_int), value :: h1, h2, ho
        end subroutine c_fdace_mul

        subroutine c_fdace_div(h1, h2, ho) bind(C, name="fdace_div")
            import :: c_int; integer(c_int), value :: h1, h2, ho
        end subroutine c_fdace_div

        subroutine c_fdace_add_double(h1, val, ho) bind(C, name="fdace_add_double")
            import :: c_int, c_double; integer(c_int), value :: h1, ho; real(c_double), value :: val
        end subroutine c_fdace_add_double

        subroutine c_fdace_mul_double(h1, val, ho) bind(C, name="fdace_mul_double")
            import :: c_int, c_double; integer(c_int), value :: h1, ho; real(c_double), value :: val
        end subroutine c_fdace_mul_double

        subroutine c_fdace_negate(hi, ho) bind(C, name="fdace_negate")
            import :: c_int; integer(c_int), value :: hi, ho
        end subroutine c_fdace_negate

        subroutine c_fdace_sub_double(h1, val, ho) bind(C, name="fdace_sub_double")
            import :: c_int, c_double; integer(c_int), value :: h1, ho; real(c_double), value :: val
        end subroutine c_fdace_sub_double

        subroutine c_fdace_double_sub(val, h2, ho) bind(C, name="fdace_double_sub")
            import :: c_int, c_double; integer(c_int), value :: h2, ho; real(c_double), value :: val
        end subroutine c_fdace_double_sub

        subroutine c_fdace_double_div(val, h2, ho) bind(C, name="fdace_double_div")
            import :: c_int, c_double; integer(c_int), value :: h2, ho; real(c_double), value :: val
        end subroutine c_fdace_double_div

        subroutine c_fdace_sin(hi, ho) bind(C, name="fdace_sin")
            import :: c_int; integer(c_int), value :: hi, ho
        end subroutine c_fdace_sin

        ! [新增] cos, exp, sqrt 的 C 接口
        subroutine c_fdace_cos(hi, ho) bind(C, name="fdace_cos")
            import :: c_int; integer(c_int), value :: hi, ho
        end subroutine c_fdace_cos

        subroutine c_fdace_exp(hi, ho) bind(C, name="fdace_exp")
            import :: c_int; integer(c_int), value :: hi, ho
        end subroutine c_fdace_exp

        subroutine c_fdace_sqrt(hi, ho) bind(C, name="fdace_sqrt")
            import :: c_int; integer(c_int), value :: hi, ho
        end subroutine c_fdace_sqrt

        subroutine c_fdace_deriv(hi, var_idx, ho) bind(C, name="fdace_deriv")
            import :: c_int; integer(c_int), value :: hi, var_idx, ho
        end subroutine c_fdace_deriv

        subroutine c_fdace_get_cons(h, val) bind(C, name="fdace_get_cons")
            import :: c_int, c_double; integer(c_int), value :: h; real(c_double), intent(out) :: val
        end subroutine c_fdace_get_cons

        subroutine c_fdace_print(h) bind(C, name="fdace_print")
            import :: c_int; integer(c_int), value :: h
        end subroutine c_fdace_print

        subroutine c_fdace_eval_var(hi, var_idx, val, ho) bind(C, name="fdace_eval_var")
            import :: c_int, c_double
            integer(c_int), value :: hi, var_idx, ho
            real(c_double), value :: val
        end subroutine c_fdace_eval_var

        subroutine c_fdace_eval_all(hi, vals, n, res) bind(C, name="fdace_eval_all")
            import :: c_int, c_double
            integer(c_int), value :: hi, n
            real(c_double), intent(in) :: vals(*)
            real(c_double), intent(out) :: res
        end subroutine c_fdace_eval_all
    end interface

    ! 在模块顶部类型定义区补充
    type :: CompiledDA
        integer(c_int) :: handle = -1
        integer :: dim = 0  ! 记录编译了几个多项式 (比如位置速度是 6)
    contains
        procedure :: destroy => compiled_destroy
        procedure :: eval => compiled_eval
    end type CompiledDA

    ! =========================================================
    ! 2. 面向对象 DA 类
    ! =========================================================
    type :: DA
        integer(c_int) :: handle = -1
    contains
        procedure :: init => da_init
        procedure :: init_var => da_init_var
        procedure :: destroy => da_destroy
        procedure :: print => da_print
        procedure :: cons => da_get_cons
        procedure :: deriv => da_deriv
        procedure, private :: da_eval_var
        procedure, private :: da_eval_all
        generic :: eval => da_eval_var, da_eval_all
    end type DA

    ! ---------------------------------------------------------
    ! 算符重载：标量 DA 与 实数
    ! ---------------------------------------------------------
    interface operator(+)
        module procedure da_add_da, da_add_real, real_add_da
    end interface
    
    interface operator(-)
        ! 包含了一元负号 (例如: -x)
        module procedure da_sub_da, da_sub_real, real_sub_da, unary_minus_da
    end interface
    
    interface operator(*)
        module procedure da_mul_da, da_mul_real, real_mul_da
    end interface
    
    interface operator(/)
        module procedure da_div_da, da_div_real, real_div_da
    end interface
    
    interface assignment(=)
        module procedure da_assign_da, da_assign_real
    end interface

    
    ! 科学函数重载
    interface sin
        module procedure da_sin
    end interface
    interface cos
        module procedure da_cos
    end interface
    interface exp
        module procedure da_exp
    end interface
    interface sqrt
        module procedure da_sqrt
    end interface

    ! =========================================================
    ! 3. 向量类 (AlgebraicVector)
    ! =========================================================
    type :: AlgebraicVector
        type(DA), allocatable :: elements(:)
        integer :: size = 0
    contains
        procedure :: init => vector_init
        procedure :: destroy => vector_destroy
        procedure :: print => vector_print
        procedure :: get => vector_get
        procedure :: set => vector_set
        procedure :: compile => vector_compile
        procedure, private :: vector_eval_var
        procedure, private :: vector_eval_all
        generic :: eval => vector_eval_var, vector_eval_all

    end type AlgebraicVector

    ! ---------------------------------------------------------
    ! 算符重载：向量与 标量(DA/实数)
    ! ---------------------------------------------------------
    interface operator(+)
        module procedure vector_add_vector
    end interface
    
    interface operator(-)
        ! 包含了一元负号 (例如: -v)
        module procedure vector_sub_vector, unary_minus_vector
    end interface
    
    interface operator(*)
        ! 向量点乘
        module procedure vector_dot_vector
        ! 标量 * 向量
        module procedure real_mul_vector, da_mul_vector
        ! 向量 * 标量 (交换律)
        module procedure vector_mul_real, vector_mul_da
    end interface
    
    interface operator(/)
        ! 向量 / 标量
        module procedure vector_div_real, vector_div_da
    end interface
    
    interface assignment(=)
        module procedure vector_assign_vector
    end interface

contains
    ! [新增] 全局初始化子程序
    subroutine dace_initialize(order, vars)
        integer, intent(in) :: order, vars
        write(*, *) '[DACE Engine] Initializing... Max Order = ', order, ', Variables = ', vars
        call c_daceInit(int(order, c_int), int(vars, c_int))
        write(*, *) '[DACE Engine] Initialization Complete.'
    end subroutine dace_initialize

    ! --- DA 生命周期管理 ---
    subroutine da_init(this)
        class(DA), intent(inout) :: this
        if (this%handle == -1) call c_fdace_allocate(this%handle)
    end subroutine da_init

    subroutine da_init_var(this, var_idx)
        class(DA), intent(inout) :: this
        integer, intent(in) :: var_idx
        if (this%handle == -1) call c_fdace_allocate_var(this%handle, int(var_idx, c_int))
    end subroutine da_init_var

    subroutine da_destroy(this)
        class(DA), intent(inout) :: this
        if (this%handle /= -1) then
            call c_fdace_free(this%handle)
            this%handle = -1
        end if
    end subroutine da_destroy

    ! === CompiledDA 的生命周期与求值 ===
    subroutine compiled_destroy(this)
        class(CompiledDA), intent(inout) :: this
        if (this%handle /= -1) then
            call c_fdace_compiled_free(this%handle)
            this%handle = -1
        end if
    end subroutine compiled_destroy

    function compiled_eval(this, args) result(res)
        class(CompiledDA), intent(in) :: this
        real(8), intent(in) :: args(:)
        real(8), allocatable :: res(:)
        
        allocate(res(this%dim))
        ! 调用 C++ 端的极致优化模板
        call c_fdace_compiled_eval_double(this%handle, args, size(args), res, this%dim)
    end function compiled_eval

    ! === AlgebraicVector 的编译方法 ===
    function vector_compile(this) result(cda)
        class(AlgebraicVector), intent(in) :: this
        type(CompiledDA) :: cda
        integer(c_int), allocatable :: handles(:)
        integer :: i
        
        ! 提取向量中所有 DA 的底层句柄
        allocate(handles(this%size))
        do i = 1, this%size
            handles(i) = this%elements(i)%handle
        end do
        
        ! 批量编译
        call c_fdace_compile_vector(handles, this%size, cda%handle)
        cda%dim = this%size
    end function vector_compile

    ! --- DA 赋值重载 (=) ---
    subroutine da_assign_da(lhs, rhs)
        class(DA), intent(inout) :: lhs
        class(DA), intent(in) :: rhs
        if (lhs%handle == -1) call lhs%init()
        call c_fdace_copy(rhs%handle, lhs%handle)
    end subroutine da_assign_da

    subroutine da_assign_real(lhs, val)
        class(DA), intent(inout) :: lhs
        real(8), intent(in) :: val
        if (lhs%handle == -1) call lhs%init()
        call c_fdace_add_double(lhs%handle, val, lhs%handle)
    end subroutine da_assign_real

    ! --- DA 数学算符 ---
    type(DA) function da_add_da(da1, da2) result(res)
        class(DA), intent(in) :: da1, da2
        call res%init()
        call c_fdace_add(da1%handle, da2%handle, res%handle)
    end function da_add_da

    type(DA) function da_sub_da(da1, da2) result(res)
        class(DA), intent(in) :: da1, da2
        call res%init()
        call c_fdace_sub(da1%handle, da2%handle, res%handle)
    end function da_sub_da

    type(DA) function da_mul_da(da1, da2) result(res)
        class(DA), intent(in) :: da1, da2
        call res%init()
        call c_fdace_mul(da1%handle, da2%handle, res%handle)
    end function da_mul_da

    type(DA) function da_div_da(da1, da2) result(res)
        class(DA), intent(in) :: da1, da2
        call res%init()
        call c_fdace_div(da1%handle, da2%handle, res%handle)
    end function da_div_da

    type(DA) function da_add_real(da1, val) result(res)
        class(DA), intent(in) :: da1
        real(8), intent(in) :: val
        call res%init()
        call c_fdace_add_double(da1%handle, val, res%handle)
    end function da_add_real

    type(DA) function real_add_da(val, da1) result(res)
        real(8), intent(in) :: val
        class(DA), intent(in) :: da1
        res = da_add_real(da1, val)
    end function real_add_da

    type(DA) function da_mul_real(da1, val) result(res)
        class(DA), intent(in) :: da1
        real(8), intent(in) :: val
        call res%init()
        call c_fdace_mul_double(da1%handle, val, res%handle)
    end function da_mul_real

    type(DA) function real_mul_da(val, da1) result(res)
        real(8), intent(in) :: val
        class(DA), intent(in) :: da1
        res = da_mul_real(da1, val)
    end function real_mul_da

    ! --- 一元负号: -DA ---
    type(DA) function unary_minus_da(da1) result(res)
        class(DA), intent(in) :: da1
        call res%init()
        call c_fdace_negate(da1%handle, res%handle)
    end function unary_minus_da

    ! --- DA - 实数 ---
    type(DA) function da_sub_real(da1, val) result(res)
        class(DA), intent(in) :: da1
        real(8), intent(in) :: val
        call res%init()
        call c_fdace_sub_double(da1%handle, val, res%handle)
    end function da_sub_real

    ! --- 实数 - DA ---
    type(DA) function real_sub_da(val, da1) result(res)
        real(8), intent(in) :: val
        class(DA), intent(in) :: da1
        call res%init()
        call c_fdace_double_sub(val, da1%handle, res%handle)
    end function real_sub_da

    ! --- 实数 / DA ---
    type(DA) function real_div_da(val, da2) result(res)
        real(8), intent(in) :: val
        class(DA), intent(in) :: da2
        call res%init()
        call c_fdace_double_div(val, da2%handle, res%handle)
    end function real_div_da

    ! --- DA 除以 实数 ---
    type(DA) function da_div_real(da1, val) result(res)
        class(DA), intent(in) :: da1
        real(8), intent(in) :: val
        res = da1 * (1.0d0 / val)
    end function da_div_real

    ! --- 科学函数与微积分 ---
    type(DA) function da_sin(da1) result(res)
        class(DA), intent(in) :: da1
        call res%init()
        call c_fdace_sin(da1%handle, res%handle)
    end function da_sin

    type(DA) function da_cos(da1) result(res)
        class(DA), intent(in) :: da1
        call res%init()
        call c_fdace_cos(da1%handle, res%handle)
    end function da_cos

    type(DA) function da_exp(da1) result(res)
        class(DA), intent(in) :: da1
        call res%init()
        call c_fdace_exp(da1%handle, res%handle)
    end function da_exp

    type(DA) function da_sqrt(da1) result(res)
        class(DA), intent(in) :: da1
        call res%init()
        call c_fdace_sqrt(da1%handle, res%handle)
    end function da_sqrt

    type(DA) function da_deriv(this, var_idx) result(res)
        class(DA), intent(in) :: this
        integer, intent(in) :: var_idx
        call res%init()
        call c_fdace_deriv(this%handle, int(var_idx, c_int), res%handle)
    end function da_deriv

    ! --- 工具函数 ---
    real(8) function da_get_cons(this)
        class(DA), intent(in) :: this
        call c_fdace_get_cons(this%handle, da_get_cons)
    end function da_get_cons

    subroutine da_print(this)
        class(DA), intent(in) :: this
        call c_fdace_print(this%handle)
    end subroutine da_print

    ! =========================================================
    ! AlgebraicVector 实现
    ! =========================================================
    subroutine vector_init(this, n)
        class(AlgebraicVector), intent(inout) :: this
        integer, intent(in) :: n
        integer :: i
        this%size = n
        allocate(this%elements(n))
        do i = 1, n
            call this%elements(i)%init()
        end do
    end subroutine vector_init

    subroutine vector_destroy(this)
        class(AlgebraicVector), intent(inout) :: this
        integer :: i
        if (allocated(this%elements)) then
            do i = 1, this%size
                call this%elements(i)%destroy()
            end do
            deallocate(this%elements)
            this%size = 0
        end if
    end subroutine vector_destroy

    subroutine vector_print(this)
        class(AlgebraicVector), intent(in) :: this
        integer :: i
        do i = 1, this%size
            write(*,"(A,I0,A)", advance="no") "Vector[", i, "] = "
            call this%elements(i)%print()
        end do
    end subroutine vector_print

    type(DA) function vector_get(this, i)
        class(AlgebraicVector), intent(in) :: this
        integer, intent(in) :: i
        vector_get = this%elements(i)
    end function vector_get

    subroutine vector_set(this, i, val)
        class(AlgebraicVector), intent(inout) :: this
        integer, intent(in) :: i
        type(DA), intent(in) :: val
        this%elements(i) = val
    end subroutine vector_set

    subroutine vector_assign_vector(lhs, rhs)
        class(AlgebraicVector), intent(inout) :: lhs
        class(AlgebraicVector), intent(in) :: rhs
        integer :: i
        if (lhs%size /= rhs%size) then
            call lhs%destroy()
            call lhs%init(rhs%size)
        end if
        do i = 1, rhs%size
            lhs%elements(i) = rhs%elements(i)
        end do
    end subroutine vector_assign_vector

    type(AlgebraicVector) function vector_sub_vector(v1, v2) result(res)
        class(AlgebraicVector), intent(in) :: v1, v2
        integer :: i
        if (v1%size /= v2%size) stop "ERROR: Vector Subtraction Dimension Mismatch!"
        call res%init(v1%size)
        do i = 1, v1%size
            res%elements(i) = v1%elements(i) - v2%elements(i)
        end do
    end function vector_sub_vector

    type(AlgebraicVector) function vector_add_vector(v1, v2) result(res)
        class(AlgebraicVector), intent(in) :: v1, v2
        integer :: i
        call res%init(v1%size)
        do i = 1, v1%size
            res%elements(i) = v1%elements(i) + v2%elements(i)
        end do
    end function vector_add_vector

    type(AlgebraicVector) function real_mul_vector(val, vec) result(res)
        real(8), intent(in) :: val
        class(AlgebraicVector), intent(in) :: vec
        integer :: i
        call res%init(vec%size)
        do i = 1, vec%size
            res%elements(i) = val * vec%elements(i)
        end do
    end function real_mul_vector

    type(AlgebraicVector) function da_mul_vector(val, vec) result(res)
        class(DA), intent(in) :: val
        class(AlgebraicVector), intent(in) :: vec
        integer :: i
        call res%init(vec%size)
        do i = 1, vec%size
            res%elements(i) = val * vec%elements(i)
        end do
    end function da_mul_vector

    type(DA) function vector_dot_vector(v1, v2) result(res)
        class(AlgebraicVector), intent(in) :: v1, v2
        integer :: i
        if (v1%size /= v2%size) stop "ERROR: Vector Dot Product Dimension Mismatch!"
        call res%init() 
        res = 0.0d0 
        do i = 1, v1%size
            res = res + (v1%elements(i) * v2%elements(i))
        end do
    end function vector_dot_vector

    ! 一元负号: -Vector
    type(AlgebraicVector) function unary_minus_vector(v1) result(res)
        class(AlgebraicVector), intent(in) :: v1
        integer :: i
        call res%init(v1%size)
        do i = 1, v1%size
            res%elements(i) = -v1%elements(i) ! 自动调用 DA 的一元负号
        end do
    end function unary_minus_vector

    ! 向量 * 实数 (补齐交换律)
    type(AlgebraicVector) function vector_mul_real(vec, val) result(res)
        class(AlgebraicVector), intent(in) :: vec
        real(8), intent(in) :: val
        res = val * vec  ! 直接复用已有的 real_mul_vector
    end function vector_mul_real

    ! 向量 * DA (补齐交换律)
    type(AlgebraicVector) function vector_mul_da(vec, val) result(res)
        class(AlgebraicVector), intent(in) :: vec
        class(DA), intent(in) :: val
        res = val * vec  ! 直接复用已有的 da_mul_vector
    end function vector_mul_da

    ! 向量 / 实数
    type(AlgebraicVector) function vector_div_real(vec, val) result(res)
        class(AlgebraicVector), intent(in) :: vec
        real(8), intent(in) :: val
        res = vec * (1.0d0 / val)
    end function vector_div_real

    ! 向量 / DA
    type(AlgebraicVector) function vector_div_da(vec, val) result(res)
        class(AlgebraicVector), intent(in) :: vec
        class(DA), intent(in) :: val
        integer :: i
        call res%init(vec%size)
        do i = 1, vec%size
            res%elements(i) = vec%elements(i) / val ! 自动调用 DA / DA
        end do
    end function vector_div_da

    ! --- DA 的 eval 实现 ---
    type(DA) function da_eval_var(this, var_idx, val) result(res)
        class(DA), intent(in) :: this
        integer, intent(in) :: var_idx
        real(8), intent(in) :: val
        
        call res%init()
        call c_fdace_eval_var(this%handle, int(var_idx, c_int), val, res%handle)
    end function da_eval_var

    ! --- Vector 的 eval 实现 (天然的数组循环优势) ---
    type(AlgebraicVector) function vector_eval_var(this, var_idx, val) result(res)
        class(AlgebraicVector), intent(in) :: this
        integer, intent(in) :: var_idx
        real(8), intent(in) :: val
        integer :: i
        
        call res%init(this%size)
        do i = 1, this%size
            ! 对向量里的每一个 DA 元素执行 eval
            res%elements(i) = this%elements(i)%eval(var_idx, val)
        end do
    end function vector_eval_var

    ! --- DA 的全代入实现 (返回实数标量) ---
    real(8) function da_eval_all(this, vals)
        class(DA), intent(in) :: this
        real(8), intent(in) :: vals(:) ! 传入形如 [0.5d0, 1.2d0, ...] 的数组
        
        call c_fdace_eval_all(this%handle, vals, size(vals), da_eval_all)
    end function da_eval_all

    ! --- Vector 的全代入实现 (返回实数数组) ---
    ! --- Vector 的全代入实现 (单次极速版，替换原有的 do 循环) ---
    function vector_eval_all(this, vals) result(res)
        ! 只适用于单次求值场景
        !! 如果要进行大规模MC求值，请调用compile方法预编译后再用compile eval
        class(AlgebraicVector), intent(in) :: this
        real(8), intent(in) :: vals(:)
        real(8), allocatable :: res(:)
        
        integer(c_int) :: cda_handle
        integer(c_int), allocatable :: handles(:)
        integer :: i
        
        ! 1. 提取向量中所有 DA 的底层句柄
        allocate(handles(this%size))
        do i = 1, this%size
            handles(i) = this%elements(i)%handle
        end do
        
        ! 2. 整体编译：只需调用 1 次 C 接口，将 6 个多项式一次性打包展平
        call c_fdace_compile_vector(handles, this%size, cda_handle)
        
        ! 3. 极速求值：调用 C++ 底层的 double 批量求值模板
        allocate(res(this%size))
        call c_fdace_compiled_eval_double(cda_handle, vals, size(vals), res, this%size)
        
        ! 4. 立即释放临时计划
        call c_fdace_compiled_free(cda_handle)
    end function vector_eval_all


end module pod_dace_classes