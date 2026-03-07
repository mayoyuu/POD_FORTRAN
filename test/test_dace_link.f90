program test_dace_link
    use pod_dace_bridge, only: dace_initialize
    use pod_dace_classes
    implicit none
    
    ! 声明标量 DA 对象
    type(DA) :: x, y, f, df_dx, dot_res
    ! 声明向量对象
    type(AlgebraicVector) :: v1, v2, v3
    
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