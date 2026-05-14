! ==============================================================================
! test_vector_leak.f90
! 测试目标：证明 DA 本身无泄露，泄露源于 AlgebraicVector 算符重载产生的隐式临时对象。
! ==============================================================================

program test_vector_leak
    use pod_global, only: DP
    use pod_dace_classes
    implicit none

    integer, parameter :: n_iter = 1000
    integer :: i
    type(AlgebraicVector) :: v1, v2, v3, temp_vec
    integer :: count_start, count_end

    call dace_initialize(2, 3)

    ! 初始化 3 维测试向量
    call v1%init(3)
    call v2%init(3)
    call v3%init(3)
    call temp_vec%init(3)

    print *, "============================================================"
    print *, "   DA Vector Operator Leak Verification Test"
    print *, "============================================================"

    ! ------------------------------------------------------------------
    ! 测试 1: 重现泄露 (算符重载作为子程序参数)
    ! ------------------------------------------------------------------
    count_start = active_da_count()
    do i = 1, n_iter
        ! 这行代码完美复刻了你在引力场 f_zonal_da 中的写法：
        ! `0.0_DP * v1` 会触发 real_mul_vector 生成一个隐式临时 AlgebraicVector 对象。
        ! 放入 call 语句后，Fortran 编译器极大概率会“遗忘”销毁这个临时对象。
        call vec_add(v2, 0.0_DP * v1, v3)
    end do
    count_end = active_da_count()

    print *, ">>> Test 1: Using overloaded operator (0.0_DP * v1) in args"
    print *, "    Active handles before: ", count_start
    print *, "    Active handles after : ", count_end
    print *, "    *** Leaked handles *** : ", count_end - count_start
    print *, "    (预期: 泄露 ", n_iter * 3, " 个句柄)"
    print *, ""

    ! ------------------------------------------------------------------
    ! 测试 2: 纯 Subroutine 安全写法
    ! ------------------------------------------------------------------
    count_start = active_da_count()
    do i = 1, n_iter
        ! 彻底抛弃算符，将隐式的临时变量“显式化”并使用无临时变量的接口
        ! 等价于: temp_vec = 0.0_DP * v1
        call vec_mul(0.0_DP, v1, temp_vec)  
        ! 等价于: v3 = v2 + temp_vec
        call vec_add(v2, temp_vec, v3)      
    end do
    count_end = active_da_count()

    print *, ">>> Test 2: Using pure _sub routines (vec_mul + vec_add)"
    print *, "    Active handles before: ", count_start
    print *, "    Active handles after : ", count_end
    print *, "    ✓ Leaked handles     : ", count_end - count_start
    print *, "    (预期: 0 泄露)"

    ! 清理
    call v1%destroy()
    call v2%destroy()
    call v3%destroy()
    call temp_vec%destroy()

end program test_vector_leak