program test_orbit_propagation
    use pod_engine_module, only: pod_engine_init
    use pod_global, only: DP
    use pod_spice, only: str2et
    use pod_force_model_module, only: init_gravity_network
    
    ! 仅仅引入顶层传播模块
    use pod_orbit_propagation, only: orbit_state, propagation_result, &
                                     propagate_orbit, display_propagation_results, &
                                     save_propagation_results, cleanup_propagation_result
    use pod_da_integrator_module, only: METHOD_RKF45, METHOD_RKF78
    
    implicit none
    
    type(orbit_state) :: initial_state
    type(propagation_result) :: result
    
    write(*,*) "=================================================="
    write(*,*) "       POD System - 轨道传播顶层集成测试        "
    write(*,*) "=================================================="
    
    ! =========================================================
    ! 1. 系统底层初始化 
    ! =========================================================
    call pod_engine_init('dummy_test_config.txt')
    
    ! =========================================================
    ! 2. 用户调用层 (面向物理真实量)
    ! =========================================================
    
    ! A. 准备初始状态 ( km, km/s, UTC)
    write(*,*) ">>> 正在装载轨道初始状态..."
    call str2et('2024-03-09T12:00:00', initial_state%epoch)
    initial_state%state = [100000.0_DP, 50000.0_DP, 20000.0_DP, &  ! 位置 (km)
                                1.5_DP,      2.5_DP,      0.5_DP]  ! 速度 (km/s)
    
    ! B. 传播
    ! 参数说明: 
    ! - initial_state: 刚才设置的物理状态
    ! - 86400.0_DP:    总共传播 1 天 (86400秒)
    ! - 2:             选择 RKF78 (1代表RKF45，2代表RKF78)
    ! - result:        用于接收输出结果的结构体
    write(*,*) ">>> 开始轨道积分..."
    call propagate_orbit(initial_state, 86400.0_DP, METHOD_RKF78, result)
    
    ! =========================================================
    ! 3. 结果处理
    ! =========================================================
    
    ! 在屏幕上打印摘要
    call display_propagation_results(result)
    
    ! 一键导出 CSV 文件到 results 目录
    call save_propagation_results(result)
    
    ! 释放内部产生的巨量数组内存
    call cleanup_propagation_result(result)
    
    write(*,*) "测试顺利结束！"

end program test_orbit_propagation