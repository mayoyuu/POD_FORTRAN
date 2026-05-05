program test_da_orbit_propagation
    ! 1. 引入引擎点火 API
    use pod_engine_module, only: pod_engine_init
    use pod_global, only: DP
    use pod_spice, only: str2et
    
    ! 仅仅引入顶层 DA 传播模块
    use pod_da_orbit_propagation, only: da_orbit_state, da_propagation_result, &
                                        propagate_da_orbit, display_da_propagation_results, &
                                        save_da_propagation_results, cleanup_da_propagation_result
    use pod_da_integrator_module, only: METHOD_RKF45, METHOD_RKF78
    
    implicit none
    
    type(da_orbit_state) :: initial_state
    type(da_propagation_result) :: result
    real(DP) :: propagation_time, epoch_end
    
    write(*,*) "=================================================="
    write(*,*) "      POD System - DA 轨道传播顶层集成测试       "
    write(*,*) "=================================================="
    
    ! =========================================================
    ! 1. 系统底层初始化 
    ! =========================================================
    call pod_engine_init('dummy_test_config.txt')
    
    ! =========================================================
    ! 2. 用户调用层 
    ! =========================================================
    
    ! A. 准备初始状态 (标称值 km, km/s, UTC, 以及极度关键的 DA 阶数)
    ! "X": -321297.69007257116,
    ! "X_DOT": -0.8302491438122398,
    ! "Y": 133068.22288313005,
    ! "Y_DOT": -0.6171754616163115,
    ! "Z": 73161.5765643454,
    ! "Z_DOT": -0.3597639886185437
    write(*,*) ">>> 正在装载 DA 轨道初始状态..."
    call str2et('2025-12-15T00:00:01.000000', initial_state%epoch)
    call str2et('2026-01-28T17:22:00.999000 ', epoch_end)
    ! epoch_end = 818669185.98332703 

    ! DRO UTC 2025 12 15 00 00   1.000000  
    ! -402779.291910,-181111.455576,-109249.167033,0.291707,-0.504100,-0.263272

    ! DRO UTC 2026 01 28 17 22   0.999000   
    ! 136191.971183  343565.836997  189807.408288      -0.706857       0.168891       0.065850


    write(*,*) ">>> 初始状态历元: ", initial_state%epoch
    write(*,*) ">>> 传播结束历元: ", epoch_end
    
    initial_state%nominal_state = [-402779.291910_DP,-181111.455576_DP,-109249.167033_DP, &  ! 位置标称值 (km)
                                        0.291707_DP,-0.504100_DP,-0.263272_DP]  ! 速度标称值 (km/s)
                                        
    ! ✨ 唯一比普通传播多出来的一行配置：告诉系统你想算到几阶泰勒展开！
    ! 设为 1 阶通常用于提取 STM；设为 2~5 阶用于高阶不确定性/蒙特卡洛替代模型
    initial_state%da_order = 3
    


    propagation_time = epoch_end - initial_state%epoch


    ! B. 传播
    ! 参数说明: 
    ! - initial_state: 刚才设置的物理状态 + DA 阶数
    ! - propagation_time:    总共传播时间
    ! - 2:             选择 RKF78 (1代表RKF45，2代表RKF78)
    ! - result:        用于接收输出结果的结构体 (里面装的是 DA 向量)
    write(*,*) ">>> 开始微分代数 (DA) 轨道积分..."
    call propagate_da_orbit(initial_state, propagation_time, METHOD_RKF78, result)
    
    ! =========================================================
    ! 3. 结果处理
    ! =========================================================
    
    ! 在屏幕上打印摘要 (会自动打印出 0 阶标称轨迹，以及左上角 3x3 物理 STM)
    call display_da_propagation_results(result)
    
    ! 一键导出标称轨迹 CSV 文件到 output/ 目录
    call save_da_propagation_results(result)
    
    ! 释放内部产生的巨量 DA 代数多项式内存与 C++ 句柄
    call cleanup_da_propagation_result(result)
    
    write(*,*) "DA 顶层测试顺利结束！"

end program test_da_orbit_propagation