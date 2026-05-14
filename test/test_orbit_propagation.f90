program test_orbit_propagation
    use pod_engine_module, only: pod_engine_init
    use pod_global, only: DP
    use pod_spice, only: str2et
    
    ! 仅仅引入顶层传播模块
    use pod_orbit_propagation, only: orbit_state, propagation_result, &
                                     propagate_orbit, display_propagation_results, &
                                     save_propagation_results, cleanup_propagation_result
    use pod_da_integrator_module, only: METHOD_RKF45, METHOD_RKF78
    
    implicit none
    
    type(orbit_state) :: initial_state
    type(propagation_result) :: result
    real(DP) :: propagation_time, epoch_end
    
    write(*,*) "=================================================="
    write(*,*) "       POD System - 轨道传播顶层集成测试        "
    write(*,*) "=================================================="
    
    ! =========================================================
    ! 1. 系统底层初始化 
    ! =========================================================
    call pod_engine_init('config/dummy_test_config.txt')
    
    ! =========================================================
    ! 2. 用户调用层 (面向物理真实量)
    ! =========================================================
    
    ! A. 准备初始状态 ( km, km/s, UTC)
    write(*,*) ">>> 正在装载轨道初始状态..."
    call str2et('2025-12-15T00:00:01.000000', initial_state%epoch)
    call str2et('2026-01-01T17:22:00.999000 ', epoch_end)
   
    
    initial_state%state = [-402779.291910_DP,-181111.455576_DP,-109249.167033_DP, &  ! 位置标称值 (km)
                                        0.291707_DP,-0.504100_DP,-0.263272_DP]  ! 速度标称值 (km/s)
                    
    propagation_time = epoch_end - initial_state%epoch
    ! B. 传播
    ! 参数说明: 
    ! - initial_state: 刚才设置的物理状态
    ! - 86400.0_DP:    总共传播 1 天 (86400秒)
    ! - 2:             选择 RKF78 (1代表RKF45，2代表RKF78)
    ! - result:        用于接收输出结果的结构体
    write(*,*) ">>> 开始轨道积分..."
    call propagate_orbit(initial_state, propagation_time, METHOD_RKF78, result)
    
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