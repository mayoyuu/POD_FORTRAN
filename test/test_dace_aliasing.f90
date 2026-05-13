program test_dace_final
    use :: pod_global, only: DP
    use :: pod_dace_classes
    implicit none

    type(AlgebraicVector) :: v_state, v_k, v_ref
    type(DA) :: tmp_var
    real(DP) :: dt_beta
    integer :: i, n, count_before, count_after
    real(DP), allocatable :: res_val(:), ref_val(:)

    ! 1. 初始化 DACE 引擎
    call dace_initialize(2, 3)
    n = 3

    print *, "========================================="
    print *, " 验证 1: 内存泄露测试 (确认 da_var 的危害)"
    print *, "========================================="
    
    ! --- 情况 A: 演示泄露 ---
    count_before = active_da_count()
    do i = 1, 100
        ! 重点：da_var(1) 返回的对象包含一个新句柄，但赋值后该对象消失，句柄无法 destroy
        call da_add(da_var(1), 1.0_DP, tmp_var) 
    end do
    count_after = active_da_count()
    print *, "直接调用 da_var 100 次，句柄增加: ", count_after - count_before

    ! --- 情况 B: 演示 tmp_var 的手动控制 ---
    count_before = active_da_count()
    do i = 1, 100
        tmp_var = da_var(1)                 ! 将返回的对象接住（拷贝句柄）
        call da_add(tmp_var, 1.0_DP, tmp_var)
        call tmp_var%destroy()              ! 手动销毁，归还句柄
    end do
    count_after = active_da_count()
    print *, "使用 tmp_var 并销毁 100 次，句柄增加: ", count_after - count_before

    print *, ""
    print *, "========================================="
    print *, " 验证 2: vec_add_scaled_inplace 正确性"
    print *, "========================================="

    call v_state%init(n)
    call v_k%init(n)
    call v_ref%init(n)

    ! 初始化 v_state = [1.0+x1, 2.0+x2, 3.0+x3]
    do i = 1, n
        tmp_var = da_var(i)
        call da_add(tmp_var, real(i, DP), v_state%elements(i))
        call tmp_var%destroy() ! 必须销毁，否则初始化就泄露了
    end do

    ! 初始化 v_k = [0.5, 0.5, 0.5]
    do i = 1, n
        ! 直接设置常数项
       v_k%elements(i) = 0.5_DP
    end do

    dt_beta = 0.1_DP

    ! 3. 计算参考值 v_ref = v_state + dt_beta * v_k
    do i = 1, n
        call da_mul(v_k%elements(i), dt_beta, v_ref%elements(i))
        call da_add(v_ref%elements(i), v_state%elements(i), v_ref%elements(i))
    end do

    ! 4. 执行原位测试 (Aliasing 测试)
    ! v_state = v_state + dt_beta * v_k
    call vec_add_scaled_inplace(v_state, dt_beta, v_k, v_state)

    ! 5. 结果校验 (使用 .cons() 提取)
    allocate(res_val(n), ref_val(n))
    do i = 1, n
        res_val(i) = v_state%elements(i)%cons()
        ref_val(i) = v_ref%elements(i)%cons()
        print "(A,I1,A,F6.3,A,F6.3)", "Element ", i, ": Result = ", res_val(i), " | Expected = ", ref_val(i)
    end do

    if (all(abs(res_val - ref_val) < 1e-12)) then
        print *, "SUCCESS: vec_add_scaled_inplace 逻辑正确且地址安全!"
    else
        print *, "FAILURE: 计算结果不一致!"
    end if

    ! 6. 最终清理
    call v_state%destroy()
    call v_k%destroy()
    call v_ref%destroy()
    
    print *, "测试结束，当前活动句柄数: ", active_da_count()

end program test_dace_final