! ==============================================================================
! test_da_orbit_propagation.f90
!
! 测试目标：
!   1. 验证 DA 版轨道传播器（propagate_da_orbit）的常数项轨迹
!      与实数版轨道传播器（propagate_orbit）计算结果一致。
!   2. 验证 extract_stm_from_result 能正确提取 6×6 STM。
!   3. 验证 DA 版轨道传播器在多次调用后无内存泄露。
!
! 测试流程：
!   Phase 0: 系统初始化（开启太阳、月球、地球高阶引力场）
!   Phase 1: 设置 DRO 初始状态
!   Phase 2: RKF45 传播一致性测试（积分 4 天）
!   Phase 3: RKF78 传播一致性测试（积分 4 天）
!   Phase 4: STM 提取验证测试
!   Phase 5: DA 内存泄露测试（循环调用 propagate_da_orbit）
!   Phase 6: 清理与结果汇总
! ==============================================================================

program test_da_orbit_propagation
    use pod_global, only: DP
    use pod_engine_module, only: pod_engine_init
    use pod_config, only: config
    use pod_spice, only: str2et

    ! ---- 实数版力模型与传播器 ----
    use pod_force_model_module, only: &
        set_propagation_epoch_real => set_propagation_epoch, &
        current_epoch0_real => current_epoch0
    use pod_integrator_module, only: &
        adaptive_step_integrate, METHOD_RKF45, METHOD_RKF78
    use pod_orbit_propagation, only: &
        orbit_state, propagation_result, &
        propagate_orbit, cleanup_propagation_result

    ! ---- DA 版力模型与传播器 ----
    use pod_da_force_model_module, only: &
        set_propagation_epoch_da => set_propagation_epoch, &
        current_epoch0_da => current_epoch0, &
        init_gravity_network
    use pod_da_integrator_module, only: &
        da_adaptive_step_integrate
    use pod_da_orbit_propagation, only: &
        da_orbit_state, da_propagation_result, &
        propagate_da_orbit, display_da_propagation_results, &
        save_da_propagation_results, cleanup_da_propagation_result, &
        extract_stm_from_result

    use pod_dace_classes
    use pod_basicmath_module, only: inverse_and_determinant

    implicit none

    ! ---- 测试参数 ----
    character(len=*), parameter :: CONFIG_FILE = 'dummy_test_config.txt'
    character(len=*), parameter :: TEST_EPOCH  = '2025-12-15T00:00:01.000000'
    integer,          parameter :: DA_ORDER    = 3
    integer,          parameter :: DA_NVARS    = 6

    ! ---- 物理状态（J2000 惯性系，有量纲） ----
    real(DP), dimension(6) :: state_physical

    ! ---- 实数版传播 ----
    type(orbit_state) :: real_initial_state
    type(propagation_result) :: real_result
    real(DP) :: tdb_epoch, propagation_time

    ! ---- DA 版传播 ----
    type(da_orbit_state) :: da_initial_state
    type(da_propagation_result) :: da_result

    ! ---- 比对缓存 ----
    real(DP), dimension(6) :: real_final_state, da_final_cons
    real(DP), dimension(6,6) :: stm, stm_inv
    real(DP) :: stm_det
    integer :: info

    ! ---- 内存泄露检测 ----
    integer :: count_before, count_after, count_final
    integer :: i, j

    ! ---- 测试统计 ----
    integer :: n_pass, n_fail
    real(DP), parameter :: TOL = 1.0e-10_DP

    n_pass = 0
    n_fail = 0

    write(*,*) '============================================================'
    write(*,*) '  POD Orbit Propagation — DA vs Real 一致性 & 内存泄露测试'
    write(*,*) '============================================================'

    ! ============================================================
    ! Phase 0: 系统初始化
    ! ============================================================
    write(*,*) ''
    write(*,*) '>>> [Phase 0] 系统初始化...'

    ! 0.1 一键初始化
    call pod_engine_init(CONFIG_FILE)

    ! 0.2 强制开启复杂环境（地球 10 阶，月球 10 阶，太阳第三体）
    config%use_earth_nspheric = .true.
    config%earth_degree       = 10
    config%use_moon_nspheric  = .true.
    config%moon_degree        = 10
    config%use_third_body     = .true.
    config%use_planet(10)     = .true. ! 激活月球
    config%use_planet(11)     = .true. ! 激活太阳

    ! 0.3 初始化 DA 版引力网
    call init_gravity_network()

    ! 0.4 初始化 DACE 引擎
    call dace_initialize(DA_ORDER, DA_NVARS)

    ! 0.5 获取并设置测试历元
    call str2et(TEST_EPOCH, tdb_epoch)
    call set_propagation_epoch_real(tdb_epoch)
    call set_propagation_epoch_da(tdb_epoch)
    write(*,*) '  测试历元 (TDB): ', tdb_epoch, ' 秒'

    ! ============================================================
    ! Phase 1: 设置 DRO 初始状态
    ! ============================================================
    write(*,*) ''
    write(*,*) '>>> [Phase 1] 设置测试状态 (DRO 轨道)...'

    ! DRO UTC 2025 12 15 00 00   1.000000  
    state_physical = [-402779.291910_DP,-181111.455576_DP,-109249.167033_DP, &  ! 位置 (km)
                            0.291707_DP,-0.504100_DP,-0.263272_DP]              ! 速度 (km/s)

    ! ---- 实数版初始状态 ----
    real_initial_state%state = state_physical
    real_initial_state%epoch = tdb_epoch

    ! ---- DA 版初始状态 ----
    da_initial_state%nominal_state = state_physical
    da_initial_state%epoch         = tdb_epoch
    da_initial_state%da_order      = DA_ORDER

    propagation_time = 86400.0_DP * 4  ! 4 天

    write(*,*) '  初始状态 (km, km/s): ', state_physical
    write(*,*) '  积分总时长: 4 天 (', propagation_time, ' 秒)'

    ! ============================================================
    ! Phase 2: RKF45 传播一致性测试
    ! ============================================================
    write(*,*) ''
    write(*,*) '>>> [Phase 2] RKF45 传播一致性测试 (4 天)...'

    ! ---- 实数版 RKF45 ----
    call propagate_orbit(real_initial_state, propagation_time, 1, real_result)
    real_final_state = real_result%states(real_result%n_steps, :)
    call cleanup_propagation_result(real_result)

    ! ---- DA 版 RKF45 ----
    call propagate_da_orbit(da_initial_state, propagation_time, 1, da_result)
    da_final_cons = da_result%final_state%cons()

    ! ---- 比对结果 ----
    write(*,*) '  实数版最终状态 (km, km/s): ', real_final_state
    write(*,*) '  DA 版  最终状态 (km, km/s): ', da_final_cons

    do i = 1, 6
        if (abs(real_final_state(i) - da_final_cons(i)) < TOL) then
            n_pass = n_pass + 1
        else
            n_fail = n_fail + 1
            write(*,*) '  *** RKF45 不一致! 分量 ', i, &
                       ' 差异 = ', abs(real_final_state(i) - da_final_cons(i))
        end if
    end do

    ! 清理 DA 结果
    call cleanup_da_propagation_result(da_result)

    ! ============================================================
    ! Phase 3: RKF78 传播一致性测试
    ! ============================================================
    write(*,*) ''
    write(*,*) '>>> [Phase 3] RKF78 传播一致性测试 (4 天)...'

    ! ---- 实数版 RKF78 ----
    call propagate_orbit(real_initial_state, propagation_time, 2, real_result)
    real_final_state = real_result%states(real_result%n_steps, :)
    call cleanup_propagation_result(real_result)

    ! ---- DA 版 RKF78 ----
    call propagate_da_orbit(da_initial_state, propagation_time, 2, da_result)
    da_final_cons = da_result%final_state%cons()

    ! ---- 比对结果 ----
    write(*,*) '  实数版最终状态 (km, km/s): ', real_final_state
    write(*,*) '  DA 版  最终状态 (km, km/s): ', da_final_cons

    do i = 1, 6
        if (abs(real_final_state(i) - da_final_cons(i)) < TOL) then
            n_pass = n_pass + 1
        else
            n_fail = n_fail + 1
            write(*,*) '  *** RKF78 不一致! 分量 ', i, &
                       ' 差异 = ', abs(real_final_state(i) - da_final_cons(i))
        end if
    end do

    ! ============================================================
    ! Phase 4: STM 提取验证测试
    ! ============================================================
    write(*,*) ''
    write(*,*) '>>> [Phase 4] STM 提取验证测试...'

    ! 从 DA 结果中提取 6×6 STM
    call extract_stm_from_result(da_result, stm)

    write(*,*) '  提取的 STM (6x6):'
    do i = 1, 6
        write(*, "(4X, 6(ES15.6, 2X))") (stm(i, j), j = 1, 6)
    end do

    ! 验证 STM 行列式 ≈ 1.0（辛结构）
    call inverse_and_determinant(stm, stm_inv, stm_det, info)
    write(*,*) '  STM 行列式: ', stm_det

    if (info == 0 .and. abs(stm_det - 1.0_DP) < 1.0e-4_DP) then
        write(*,*) '  ✓ STM 行列式 ≈ 1.0，辛结构保持良好'
        n_pass = n_pass + 1
    else
        write(*,*) '  *** STM 行列式偏离 1.0 较大: ', abs(stm_det - 1.0_DP)
        n_fail = n_fail + 1
    end if

    ! 清理 DA 结果（Phase 3 的结果）
    call cleanup_da_propagation_result(da_result)

    ! ============================================================
    ! Phase 5: DA 内存泄露测试
    ! ============================================================
    write(*,*) ''
    write(*,*) '>>> [Phase 5] DA 内存泄露测试 (循环调用 propagate_da_orbit)...'

    count_before = active_da_count()
    write(*,*) '  初始活动句柄数: ', count_before

    ! 循环 10 次 DA 轨道传播（短时间积分）
    do i = 1, 10
        call propagate_da_orbit(da_initial_state, 3600.0_DP, 2, da_result)
        call cleanup_da_propagation_result(da_result)
    end do

    count_after = active_da_count()
    write(*,*) '  循环 10 次后活动句柄数: ', count_after
    write(*,*) '  句柄变化: ', count_after - count_before

    if (count_after == count_before) then
        write(*,*) '  ✓ 无内存泄露（句柄数未增长）'
        n_pass = n_pass + 1
    else
        write(*,*) '  *** 内存泄露! 句柄增加了 ', count_after - count_before
        n_fail = n_fail + 1
    end if

    ! ============================================================
    ! Phase 6: 最终清理与汇总
    ! ============================================================
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
        write(*,*) '  ✓ DA 版轨道传播器与实数版计算结果一致'
        write(*,*) '  ✓ STM 提取功能正常，辛结构保持良好'
        write(*,*) '  ✓ DA 版轨道传播器无内存泄露'
    else
        write(*,*) ''
        write(*,*) '  *** 存在测试失败，请检查上述输出 ***'
        stop 1
    end if

    write(*,*) '============================================================'

end program test_da_orbit_propagation
