program test_dace_link
    use pod_dace_bridge, only: dace_initialize
    use pod_dace_classes
    implicit none
    ! A function to calculate square root
    
    ! 声明标量 DA 对象
    type(DA) :: x, y, z, f, df_dx, dot_res, ff, f_eval
    ! 声明向量对象
    type(AlgebraicVector) :: v1, v2, v3, v_eval

    real(8) :: f_final
    real(8), allocatable :: final_vec_vals(:)
    
    write(*, *) '=========================================='
    write(*, *) '      DACE 终极运算能力与内存安全测试'
    write(*, *) '=========================================='
    write(*, *) ''
    
    ! ==========================================
    ! 阶段 1: 引擎初始化
    ! ==========================================
    ! 我们设最大阶数为 3 阶，包含 3 个独立变量
    call dace_initialize(3, 3)
    write(*, *) '>>> [阶段 1] 引擎初始化成功 (Order=3, Vars=3)'
    write(*, *) ''
    
    ! ==========================================
    ! 阶段 2: 标量微积分测试
    ! ==========================================
    write(*, *) '>>> [阶段 2] 标量 DA 微积分运算测试'
    
    ! 把 x 设为第 1 个独立变量: DA(1)
    call x%init_var(1) 
    ! 把 y 设为第 2 个独立变量: DA(2)
    call y%init_var(2) 

    call z%init_var(3)

    ! 极度优雅的数学表达式！
    ! 构造函数: f(x,y) = sin(x) * y + 2.5
    f = sin(x) * y + 2.5d0
    
    write(*, *) '构造函数 f(x,y) = sin(x)*y + 2.5 :'
    call f%print()
    
    ! 对 x 求偏导数: df/dx = cos(x) * y
    df_dx = f%deriv(1)
    write(*, *) '对 x 求导 df/dx (期望包含 cos(x)*y 的泰勒展开) :'
    call df_dx%print()
    write(*, *) ''

    ! ==========================================
    ! 阶段 3: 向量代数测试
    ! ==========================================
    write(*, *) '>>> [阶段 3] 向量 AlgebraicVector 运算测试'
    
    ! 初始化 3 维向量
    call v1%init(3)
    call v2%init(3)
    
    ! 给 v1 赋值: [x, y, 1.0]
    v1%elements(1) = x
    v1%elements(2) = y
    v1%elements(3) = 1.0d0  ! 触发 da_assign_real 重载
    
    ! 给 v2 赋值: [2.0, 3.0, x]
    v2%elements(1) = 2.0d0
    v2%elements(2) = 3.0d0
    v2%elements(3) = x
    
    write(*, *) '--- 向量 v1 ---'
    call v1%print()
    
    write(*, *) '--- 测试: 标量乘向量 (v3 = 2.0 * v1) ---'
    v3 = 2.0d0 * v1
    call v3%print()
    
    write(*, *) '--- 测试: DA变量乘向量 (v3 = y * v2) ---'
    v3 = y * v2
    call v3%print()
    
    write(*, *) '--- 测试: 向量加法 (v3 = v1 + v2) ---'
    v3 = v1 + v2
    call v3%print()
    
    write(*, *) '--- 测试: 向量点乘 (dot_res = v1 * v2) ---'
    ! 预期结果: x*2.0 + y*3.0 + 1.0*x = 3.0*x + 3.0*y
    dot_res = v1 * v2
    call dot_res%print()
    write(*, *) ''

    ! ==========================================
    ! 阶段 3.5: Eval 代入求值测试
    ! ==========================================
    write(*, *) '>>> [阶段 3.5] 代入求值 (Eval/Plug) 测试'
    
    ! 此时 x 和 y 还是存活状态，可以安全使用
    ff = 3.0d0 * x + 5.0d0 * y  
    write(*, *) '原始多项式 ff = 3.0*x + 5.0*y :'
    call ff%print()
    
    ! 把 delta_x (1号变量) 代入为 0.5
    ! f_eval 应该变成: 1.5 + 5.0 * y
    f_eval = ff%eval(1, 0.5d0) 
    write(*, *) '将 x (变量1) 代入为 0.5 后 :'
    call f_eval%print()

    f_final = f_eval%eval([0.5d0])
    write(*, *) '将 f_eval 的 y (变量1) 代入为 0.5 后 :'
    write(*, *) '将 f_eval 的所有剩余变量代入后得到的纯实数 :'
    write(*, *) f_final  ! 直接用 write 打印，不要用 %print()

    ! 向量 eval 测试
    ! 我们直接借用阶段 3 已经赋好值的 v1: [x, y, 1.0]
    write(*, *) '原始向量 v1 :'
    call v1%print()
    
    ! 一键将向量里所有的 x 替换为 0.5
    v_eval = v1%eval(1, 0.5d0)
    write(*, *) '将 v1 中的 x 代入为 0.5 后 :'
    call v_eval%print()
    write(*, *) ''

    final_vec_vals = v_eval%eval([0.0d0,0.0d0])
    write(*, *) '将 v_eval 中的 x 代入为 0.5 后 :'
    write(*, *) final_vec_vals
    

    ! ==========================================
    ! 阶段 4: 内存生命周期测试
    ! ==========================================
    write(*, *) '>>> [阶段 4] 内存自动释放与回收池测试...'
    ! 在复杂的定轨循环中，每步积分结束都需要释放临时对象
    call x%destroy()
    call y%destroy()
    call f%destroy()
    call df_dx%destroy()
    call dot_res%destroy()
    call v1%destroy()
    call v2%destroy()
    call v3%destroy()
    
    write(*, *) ''
    write(*, *) '🎉 恭喜！DA 算符与 Vector 引擎全部正常工作，内存安全退出！'
    write(*, *) '=========================================='



end program test_dace_link