!> # CAT DACE Classes Module
module pod_dace_classes
    use, intrinsic :: iso_c_binding
    implicit none
    private
    
    public :: DA, AlgebraicVector
    public :: operator(+), operator(-), operator(*), operator(/)
    public :: assignment(=)
    public :: sin, cos, exp, sqrt

    ! =========================================================
    ! 1. C 接口绑定
    ! =========================================================
    interface
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
    end interface

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
    end type DA

    ! 算符重载
    interface operator(+)
        module procedure da_add_da, da_add_real, real_add_da
    end interface
    interface operator(-)
        module procedure da_sub_da
    end interface
    interface operator(*)
        module procedure da_mul_da, da_mul_real, real_mul_da
    end interface
    interface operator(/)
        module procedure da_div_da
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
    end type AlgebraicVector

    interface operator(+)
        module procedure vector_add_vector
    end interface
    interface operator(-)
        module procedure vector_sub_vector
    end interface
    interface operator(*)
        module procedure real_mul_vector, da_mul_vector, vector_dot_vector
    end interface
    interface assignment(=)
        module procedure vector_assign_vector
    end interface

contains

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

end module pod_dace_classes