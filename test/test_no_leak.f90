program test_no_leak
    use pod_global, only: DP
    use pod_dace_classes
    implicit none

    integer, parameter :: n_iter = 200000   ! 20 万次迭代，避免太慢
    integer :: i
    type(DA) :: a, b, c
    integer :: size_start, size_end, active_start, active_end
    integer :: size_start2, size_end2, active_start2, active_end2

    ! 初始化 DACE (4 阶，2 个变量)
    call dace_initialize(4, 2)

    a = da_var(1)
    b = da_var(2)

    ! ==================================================
    ! 测试1：原始方式（隐式临时 + 无手动销毁）
    ! ==================================================
    print *, "=== Test 1: Implicit temporary, no manual destroy ==="
    size_start = da_registry_size()
    active_start = active_da_count()
    print *, "Before loop: registry = ", size_start, ", active = ", active_start

    do i = 1, n_iter
        c = a + b          ! 每步产生一个临时 DA 对象
        ! 注意：c 的旧值会在赋值时被释放（因为已经修改了赋值操作符）
        ! 但临时对象 a+b 不会被 Fortran 自动销毁
    end do

    size_end = da_registry_size()
    active_end = active_da_count()
    print *, "After loop:  registry = ", size_end, ", active = ", active_end
    print *, "Registry growth = ", size_end - size_start
    print *, "Active growth   = ", active_end - active_start

    ! 释放 c 当前持有的资源（避免干扰后续测试）
    call c%destroy()
    ! 注意：a 和 b 还保留，但 registry 不会缩小，因为 free_handles 会回收

    ! ==================================================
    ! 测试2：显式临时变量 + 手动销毁
    ! ==================================================
    print *, "=== Test 2: Explicit temporary + manual destroy ==="
    size_start2 = da_registry_size()
    active_start2 = active_da_count()
    print *, "Before loop: registry = ", size_start2, ", active = ", active_start2

    do i = 1, n_iter
        block
            type(DA) :: temp
            temp = a + b     ! 临时变量接管结果
            c = temp         ! 赋值给 c
            call temp%destroy()   ! 显式销毁临时对象，释放底层句柄
        end block
    end do

    size_end2 = da_registry_size()
    active_end2 = active_da_count()
    print *, "After loop:  registry = ", size_end2, ", active = ", active_end2
    print *, "Registry growth = ", size_end2 - size_start2
    print *, "Active growth   = ", active_end2 - active_start2

    ! 最终清理
    call a%destroy()
    call b%destroy()
    call c%destroy()

end program test_no_leak