! ==============================================================================
! test_uq_api.f90
!
! 测试目标：
!   1. 验证 DA 误差传播算法（METHOD_DA）与传统蒙特卡洛（METHOD_MC）的
!      均值和协方差结果一致性（验证 DA 泰勒展开映射的正确性）。
!   2. 验证 DA 传播算法在多次循环调用下无内存泄露。
!
! 测试流程：
!   Phase 0: 系统初始化（开启太阳、月球、地球高阶引力场与 DACE）
!   Phase 1: 设置测试轨道与初始协方差矩阵
!   Phase 2: DA vs MC 一致性测试 (传播 1 天)
!   Phase 3: DA 内存泄露测试（高频调用 run_uq_propagation）
!   Phase 4: 清理与结果汇总
! ==============================================================================

program test_uq_api
    use pod_global, only: DP
    use pod_engine_module, only: pod_engine_init
    use pod_config, only: config
    use pod_spice, only: str2et

    use pod_da_force_model_module, only: &
        init_gravity_network
        
    use pod_uq_propagation, only: &
        run_uq_propagation, METHOD_MC, METHOD_DA, METHOD_UT
    use pod_uq_state_module, only: uq_state_type
    use pod_dace_classes

    implicit none

    ! ---- 测试参数 ----
    character(len=*), parameter :: CONFIG_FILE = 'config/dummy_test_config.txt'
    character(len=*), parameter :: TEST_EPOCH  = '2024-03-09T12:00:00'
    integer,          parameter :: DA_ORDER    = 4 ! 为了匹配 MC 的非线性，可以用 3 阶
    integer,          parameter :: DA_NVARS    = 6

    ! ---- 轨道与分布状态 ----
    real(DP) :: nominal_orbit(6)
    real(DP) :: initial_covariance(6,6)
    real(DP) :: t0, tf, t_leak, epoch_start
    
    ! ---- UQ 状态对象 ----
    type(uq_state_type) :: initial_dist_da, final_dist_da
    type(uq_state_type) :: initial_dist_mc, final_dist_mc
    type(uq_state_type) :: temp_initial, temp_final

    ! ---- 杂项变量 ----
    integer :: i, count_before, count_after, count_final
    integer :: n_pass, n_fail
    real(DP) :: max_mean_diff, max_cov_diff
    
    ! MC vs DA 统计比较容差（由于 MC 的随机性，容差不能设得极小，仅作一致性范围判断）
    real(DP), parameter :: TOL_MEAN = 5.0_DP     ! 均值偏差容差 (km 或 km/s)
    real(DP), parameter :: TOL_COV  = 50.0_DP    ! 协方差偏差容差

    n_pass = 0
    n_fail = 0

    write(*,*) '============================================================'
    write(*,*) '  POD UQ API — DA vs MC 一致性 & 内存泄露测试'
    write(*,*) '============================================================'

    ! ============================================================
    ! Phase 0: 系统初始化
    ! ============================================================
    write(*,*) ''
    write(*,*) '>>> [Phase 0] 系统初始化...'

    call pod_engine_init(CONFIG_FILE)
    
    ! 开启引力场
    config%use_earth_nspheric = .true.
    config%earth_degree       = 10
    config%use_moon_nspheric  = .true.
    config%moon_degree        = 10
    call init_gravity_network()

    ! 初始化 DACE
    call dace_initialize(DA_ORDER, DA_NVARS)
    
    call str2et(TEST_EPOCH, epoch_start)
    write(*,*) '  测试历元 (TDB): ', epoch_start, ' 秒'

    ! ============================================================
    ! Phase 1: 设置测试轨道与协方差
    ! ============================================================
    write(*,*) ''
    write(*,*) '>>> [Phase 1] 设置测试状态与初始分布...'

    ! 初始轨道：典型 LTO
    nominal_orbit = [100000.0_DP, 50000.0_DP, 20000.0_DP, 1.5_DP, 2.5_DP, 0.5_DP]
    
    ! 初始协方差矩阵 (位置误差 1km，速度误差 1m/s)
    initial_covariance = 0.0_DP
    initial_covariance(1,1) = 1.0_DP**2
    initial_covariance(2,2) = 1.0_DP**2
    initial_covariance(3,3) = 1.0_DP**2
    initial_covariance(4,4) = 1.0D-3**2
    initial_covariance(5,5) = 1.0D-3**2
    initial_covariance(6,6) = 1.0D-3**2

    t0 = 0.0_DP
    tf = 86400.0_DP*4.0_DP ! 积分 4 天

    write(*,*) '  初始状态 (km, km/s): ', nominal_orbit
    write(*,*) '  粒子数量: 10000, 积分时长: 4 天'

    ! ============================================================
    ! Phase 2: DA vs MC 一致性测试
    ! ============================================================
    write(*,*) ''
    write(*,*) '>>> [Phase 2] DA vs MC 一致性测试...'

    ! ---- 运行 DA 传播 ----
    write(*,*) '  [执行 DA 算法...]'
    call run_uq_propagation( &
        nominal_state = nominal_orbit, &
        initial_cov   = initial_covariance, &
        epoch0        = epoch_start, &
        t_start       = t0, t_end = tf, &
        method_switch = METHOD_DA, &
        n_particles   = 10000, &
        save_results_to_file = .false., &
        initial_state_out = initial_dist_da, &
        da_order      = DA_ORDER, &
        final_state_out   = final_dist_da)

    ! ---- 运行 MC 传播 ----
    write(*,*) '  [执行 MC 算法...]'
    call run_uq_propagation( &
        nominal_state = nominal_orbit, &
        initial_cov   = initial_covariance, &
        epoch0        = epoch_start, &
        t_start       = t0, t_end = tf, &
        method_switch = METHOD_MC, &
        n_particles   = 10000, &
        save_results_to_file = .false., &
        initial_state_out = initial_dist_mc, &
        final_state_out   = final_dist_mc)

    ! ---- 比对统计结果 ----
    max_mean_diff = maxval(abs(final_dist_da%mean - final_dist_mc%mean))
    max_cov_diff  = maxval(abs(final_dist_da%cov - final_dist_mc%cov))

    write(*,*) '  均值最大偏差    : ', max_mean_diff
    write(*,*) '  协方差最大偏差  : ', max_cov_diff

    ! 检查均值
    if (max_mean_diff < TOL_MEAN) then
        write(*,*) '  ✓ 均值一致 (在随机散布容差内)'
        n_pass = n_pass + 1
    else
        write(*,*) '  *** 均值偏差过大!'
        n_fail = n_fail + 1
    end if

    ! 检查协方差
    if (max_cov_diff < TOL_COV) then
        write(*,*) '  ✓ 协方差一致 (在随机散布容差内)'
        n_pass = n_pass + 1
    else
        write(*,*) '  *** 协方差偏差过大!'
        n_fail = n_fail + 1
    end if

    ! 释放这部分状态数组
    if (allocated(initial_dist_da%samples)) deallocate(initial_dist_da%samples)
    if (allocated(final_dist_da%samples)) deallocate(final_dist_da%samples)
    if (allocated(initial_dist_mc%samples)) deallocate(initial_dist_mc%samples)
    if (allocated(final_dist_mc%samples)) deallocate(final_dist_mc%samples)

    ! ============================================================
    ! Phase 2.5: UT 传播测试
    ! ============================================================
    write(*,*) ''
    write(*,*) '>>> [Phase 2.5] UT 传播测试...'

    call run_uq_propagation( &
        nominal_state = nominal_orbit, &
        initial_cov   = initial_covariance, &
        epoch0        = epoch_start, &
        t_start       = t0, t_end = tf, &
        method_switch = METHOD_UT, &
        n_particles   = 100, &
        save_results_to_file = .false., &
        initial_state_out = temp_initial, &
        final_state_out   = temp_final)

    write(*,*) '  UT 传播完成'
    write(*,*) '  最终均值: ', temp_final%mean
    write(*,*) '  最终标准差: '
    write(*,*) '    Pos: ', sqrt(temp_final%cov(1,1)), sqrt(temp_final%cov(2,2)), sqrt(temp_final%cov(3,3))
    write(*,*) '    Vel: ', sqrt(temp_final%cov(4,4)), sqrt(temp_final%cov(5,5)), sqrt(temp_final%cov(6,6))

    if (any(temp_final%mean /= temp_final%mean)) then
        write(*,*) '  *** UT 均值包含 NaN!'
        n_fail = n_fail + 1
    else if (any([sqrt(temp_final%cov(1,1)), sqrt(temp_final%cov(2,2)), sqrt(temp_final%cov(3,3)), &
                   sqrt(temp_final%cov(4,4)), sqrt(temp_final%cov(5,5)), sqrt(temp_final%cov(6,6))] < 0.0_DP)) then
        write(*,*) '  *** UT 标准差为负!'
        n_fail = n_fail + 1
    else
        write(*,*) '  ✓ UT 传播结果合理'
        n_pass = n_pass + 1
    end if

    if (allocated(temp_initial%samples)) deallocate(temp_initial%samples)
    if (allocated(temp_final%samples)) deallocate(temp_final%samples)

    ! ============================================================
    ! Phase 3: DA 内存泄露测试
    ! ============================================================
    write(*,*) ''
    write(*,*) '>>> [Phase 3] DA UQ 传播内存泄露测试 (循环调用 API)...'

    count_before = active_da_count()
    t_leak = 3600.0_DP  ! 积分 1 小时即可
    
    do i = 1, 10
        call run_uq_propagation( &
            nominal_state = nominal_orbit, &
            initial_cov   = initial_covariance, &
            epoch0        = epoch_start, &
            t_start       = t0, t_end = t_leak, &
            method_switch = METHOD_DA, &
            n_particles   = 100, &          ! 测试泄露时少生成点粒子，加速测试
            save_results_to_file = .false., &
            initial_state_out = temp_initial, &
            da_order      = DA_ORDER, &
            final_state_out   = temp_final)
            
        ! Fortran 对象内存在生命周期中会覆盖或随作用域销毁，手动释放确保干干净净
        if (allocated(temp_initial%samples)) deallocate(temp_initial%samples)
        if (allocated(temp_final%samples)) deallocate(temp_final%samples)
    end do
    
    count_after = active_da_count()
    write(*,*) '  高强度调用 UQ DA API 10 次后句柄变化: ', count_after - count_before

    if (count_after == count_before) then
        write(*,*) '  ✓ 无内存泄露（复杂封装下 DA 句柄完全释放）'
        n_pass = n_pass + 1
    else
        write(*,*) '  *** 内存泄露! 句柄增加了 ', count_after - count_before
        n_fail = n_fail + 1
    end if

    ! ============================================================
    ! Phase 4: 最终清理与汇总
    ! ============================================================
    ! 清理全局引力场挂载的 DA 句柄

    count_final = active_da_count()
    write(*,*) ''
    write(*,*) '  最终活动句柄数: ', count_final
    if (count_final == 0) then
        write(*,*) '  ✓ 所有 DA 句柄已完美释放'
        n_pass = n_pass + 1
    else
        write(*,*) '  *** 仍有 ', count_final, ' 个句柄未释放'
        n_fail = n_fail + 1
    end if

    write(*,*) ''
    write(*,*) '============================================================'
    write(*,*) '  测试结果汇总'
    write(*,*) '============================================================'
    write(*,*) '  通过: ', n_pass, ' / ', n_pass + n_fail
    write(*,*) '  失败: ', n_fail

    if (n_fail == 0) then
        write(*,*) ''
        write(*,*) '  ✓ 全部测试通过!'
        write(*,*) '  ✓ DA 多项式映射的均值/协方差与 MC 分布高度一致！'
        write(*,*) '  ✓ DA 在 UQ 上层复杂结构封装中做到 0 泄露！'
    else
        write(*,*) ''
        write(*,*) '  *** 存在测试失败，请检查上述输出 ***'
        stop 1
    end if

    write(*,*) '============================================================'

end program test_uq_api