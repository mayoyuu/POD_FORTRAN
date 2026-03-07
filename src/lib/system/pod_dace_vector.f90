!> # CAT DACE Vector Module
module pod_dace_vector
    use pod_dace_core
    implicit none
    private

    public :: AlgebraicVector

    !> 微分代数向量类
    type :: AlgebraicVector
        type(DA), allocatable :: elements(:)
        integer :: size = 0
    contains
        procedure :: init => vector_init
        procedure :: get => vector_get
        procedure :: set => vector_set
    end type AlgebraicVector

    ! 重载向量加法
    interface operator(+)
        module procedure vector_add_vector
    end interface
    public :: operator(+)

contains

    ! 初始化向量长度
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

    function vector_get(this, i) result(res)
        class(AlgebraicVector), intent(in) :: this
        integer, intent(in) :: i
        type(DA) :: res
        res = this%elements(i)
    end function vector_get

    subroutine vector_set(this, i, val)
        class(AlgebraicVector), intent(inout) :: this
        integer, intent(in) :: i
        type(DA), intent(in) :: val
        this%elements(i) = val
    end subroutine vector_set

    ! 实现向量相加 (对应精密定轨中的状态向量更新)
    type(AlgebraicVector) function vector_add_vector(v1, v2) result(res)
        class(AlgebraicVector), intent(in) :: v1, v2
        integer :: i
        
        if (v1%size /= v2%size) then
            write(*,*) "Error: Vector sizes do not match!"
            stop
        end if
        
        call res%init(v1%size)
        do i = 1, v1%size
            ! 这里会自动调用 pod_dace_core 中重载的 DA + DA
            res%elements(i) = v1%elements(i) + v2%elements(i) 
        end do
    end function vector_add_vector

end module pod_dace_vector