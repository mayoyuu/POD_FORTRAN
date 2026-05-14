! ==============================================================================
! test_da_integrator.f90
!
! 测试目标：
!   1. 验证 DA 版变步长积分器（da_adaptive_step_integrate）的常数项轨迹
!      与实数版变步长积分器（adaptive_step_integrate）计算结果一致。
!   2. 验证 DA 引擎在复杂的自适应步长逻辑与循环调用下无内存泄露。
!
! 测试流程：
!   Phase 0: 系统初始化（开启太阳、月球、地球高阶引力场）
!   Phase 1: 设置地月转移轨道 (LTO) 初始状态，无量纲化
!   Phase 2: RKF45 变步长一致性测试 (积分 1 天)
!   Phase 3: RKF78 变步长一致性测试 (积分 1 天)
!   Phase 4: DA 内存泄露测试（短时变步长积分循环）
!   Phase 5: 清理与结果汇总
! ==============================================================================

program test_da_integrator
    use pod_global, only: DP
    use pod_engine_module, only: pod_engine_init
    use pod_config, only: config
    use pod_spice, only: str2et

    ! ---- 实数版力模型与积分器 ----
    use pod_force_model_module, only: &
        set_propagation_epoch_real => set_propagation_epoch, &
        current_epoch0_real => current_epoch0
    use pod_integrator_module, only: &
        adaptive_step_integrate, METHOD_RKF45, METHOD_RKF78

    ! ---- DA 版力模型与积分器 ----
    use pod_da_force_model_module, only: &
        set_propagation_epoch_da => set_propagation_epoch, &
        current_epoch0_da => current_epoch0, &
        init_gravity_network
    use pod_da_integrator_module, only: &
        da_adaptive_step_integrate

    use pod_dace_classes

    implicit none

    ! ---- 测试参数 ----
    character(len=*), parameter :: CONFIG_FILE = 'config/dummy_test_config.txt'
    character(len=*), parameter :: TEST_EPOCH  = '2024-03-09T12:00:00'
    integer,          parameter :: DA_ORDER    = 4
    integer,          parameter :: DA_NVARS    = 6

    ! ---- 物理状态（J2000 惯性系，有量纲） ----
    real(DP), dimension(6) :: state_physical

    ! ---- 无量纲状态 ----
    real(DP), dimension(6) :: state_unit
    real(DP) :: tdb_epoch, t_start_unit, t_end_unit, t_end_leak_unit

    ! ---- 积分结果缓存 ----
    real(DP), allocatable :: times_real(:), times_da(:)
    real(DP), allocatable :: states_real(:,:), nominal_states_da(:,:)
    integer :: n_steps_real, n_steps_da

    real(DP), dimension(6) :: state_rkf_real, state_rkf_da_cons
    type(AlgebraicVector) :: state_da, final_state_da

    ! ---- 杂项变量 ----
    integer :: i
    integer :: count_before, count_after, count_final
    integer :: n_pass, n_fail
    real(DP), parameter :: TOL = 1.0e-10_DP

    n_pass = 0
    n_fail = 0

    write(*,*) '============================================================'
    write(*,*) '  POD Integrator — 变步长积分 一致性 & 内存泄露测试'
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
    ! Phase 1: 设置地月转移轨道 (LTO)
    ! ============================================================
    write(*,*) ''
    write(*,*) '>>> [Phase 1] 设置测试状态 (地月转移轨道 LTO)...'

    ! 近地点高度 ~300 km，速度 ~10.83 km/s
    ! state_physical = [6678.137_DP, 0.0_DP, 0.0_DP, &   ! 位置 (km)
    !                   0.0_DP, 10.83_DP, 0.0_DP]        ! 速度 (km/s)
    state_physical = [-402779.291910_DP,-181111.455576_DP,-109249.167033_DP, &  ! 位置标称值 (km)
                                        0.291707_DP,-0.504100_DP,-0.263272_DP]        ! 速度 (km/s)

    ! 无量纲化
    state_unit(1:3) = state_physical(1:3) / config%LU
    state_unit(4:6) = state_physical(4:6) / config%VU

    ! 设定积分时长为 1 天 (86400 秒)
    t_start_unit = 0.0_DP
    t_end_unit   = 86400.0_DP*4 / config%TU

    write(*,*) '  初始状态 (km, km/s): ', state_physical
    write(*,*) '  积分总时长: 4 天 (86400*4 秒)'

    ! 初始化 DA 状态变量（附带一阶导数展开）
    call state_da%init(6)
    do i = 1, 6
        state_da%elements(i) = state_unit(i) + da_var(i)
    end do

    ! ============================================================
    ! Phase 2: RKF45 变步长一致性测试
    ! ============================================================
    write(*,*) ''
    write(*,*) '>>> [Phase 2] RKF45 变步长积分测试 (4 天)...'

    ! ---- 实数版 RKF45 ----
    call adaptive_step_integrate(state_unit, t_start_unit, t_end_unit, METHOD_RKF45, &
                                 times_real, states_real, n_steps_real)
    state_rkf_real = states_real(n_steps_real, :)

    ! ---- DA 版 RKF45 ----
    call da_adaptive_step_integrate(state_da, t_start_unit, t_end_unit, METHOD_RKF45, &
                                    times_da, nominal_states_da, final_state_da, n_steps_da)
    state_rkf_da_cons = final_state_da%cons()

    ! ---- 比对结果 ----
    write(*,*) '  实数版步数: ', n_steps_real, '步'
    write(*,*) '  DA 版步数 : ', n_steps_da, '步'
    
    if (n_steps_real == n_steps_da) then
        write(*,*) '  ✓ 变步长截断序列完全一致！'
    else
        write(*,*) '  *** 警告：变步长截断序列存在差异！'
    end if

    do i = 1, 6
        if (abs(state_rkf_real(i) - state_rkf_da_cons(i)) < TOL) then
            n_pass = n_pass + 1
        else
            n_fail = n_fail + 1
            write(*,*) '  *** RKF45 不一致! 分量 ', i, &
                       ' 差异 = ', abs(state_rkf_real(i) - state_rkf_da_cons(i))
        end if
    end do
    
    ! 释放 DA 返回的最终状态
    call final_state_da%destroy()

    ! ============================================================
    ! Phase 3: RKF78 变步长一致性测试
    ! ============================================================
    write(*,*) ''
    write(*,*) '>>> [Phase 3] RKF78 变步长积分测试 (4 天)...'

    ! ---- 实数版 RKF78 ----
    call adaptive_step_integrate(state_unit, t_start_unit, t_end_unit, METHOD_RKF78, &
                                 times_real, states_real, n_steps_real)
    state_rkf_real = states_real(n_steps_real, :)

    ! ---- DA 版 RKF78 ----
    call da_adaptive_step_integrate(state_da, t_start_unit, t_end_unit, METHOD_RKF78, &
                                    times_da, nominal_states_da, final_state_da, n_steps_da)
    state_rkf_da_cons = final_state_da%cons()

    ! ! ---- 比对结果 ----
    write(*,*) '  实数版步数: ', n_steps_real, '步'
    write(*,*) '  DA 版步数 : ', n_steps_da, '步'

    do i = 1, 6
        if (abs(state_rkf_real(i) - state_rkf_da_cons(i)) < TOL) then
            n_pass = n_pass + 1
        else
            n_fail = n_fail + 1
            write(*,*) '  *** RKF78 不一致! 分量 ', i, &
                       ' 差异 = ', abs(state_rkf_real(i) - state_rkf_da_cons(i))
        end if
    end do
    
    ! 释放 DA 返回的最终状态
    call final_state_da%destroy()

    ! ============================================================
    ! Phase 4: DA 内存泄露测试
    ! ============================================================
    write(*,*) ''
    write(*,*) '>>> [Phase 4] DA 内存泄露测试 (循环调用变步长积分)...'

    count_before = active_da_count()
    t_end_leak_unit = 3600.0_DP / config%TU  ! 只积 1 个小时
    
    ! 循环 10 次变步长积分
    do i = 1, 10
        call da_adaptive_step_integrate(state_da, t_start_unit, t_end_leak_unit, METHOD_RKF78, &
                                        times_da, nominal_states_da, final_state_da, n_steps_da)
        ! 必须手动销毁输出参数！
        call final_state_da%destroy()
    end do
    
    count_after = active_da_count()
    write(*,*) '  RKF78 变步长循环 10 次后句柄变化: ', count_after - count_before

    if (count_after == count_before) then
        write(*,*) '  ✓ 无内存泄露（复杂步长控制下句柄数未增长）'
        n_pass = n_pass + 1
    else
        write(*,*) '  *** 内存泄露! 句柄增加了 ', count_after - count_before
        n_fail = n_fail + 1
    end if

    ! ============================================================
    ! Phase 5: 最终清理与汇总
    ! ============================================================
    ! 清理自己的临时变量
    call state_da%destroy()
    
    ! ! 清理全局引力场模块中挂载的 DA 变量 (极其重要!)
    ! call earth_grav%dr_da%destroy()
    ! call moon_grav%dr_da%destroy()

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
        write(*,*) '  ✓ DA 变步长截断与实数版丝毫不差！'
        write(*,*) '  ✓ DA 引擎在高强度自适应循环中 0 泄露！'
    else
        write(*,*) ''
        write(*,*) '  *** 存在测试失败，请检查上述输出 ***'
        stop 1
    end if

    write(*,*) '============================================================'

end program test_da_integrator