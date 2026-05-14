! ==============================================================================
! test_da_force_model_consistency.f90
!
! 测试目标：
!   1. 验证 DA 版力模型（da_compute_acceleration）的常数项与实数版
!      （compute_acceleration）计算结果一致。
!   2. 验证 DA 版力模型在多次调用后无内存泄露。
!
! 测试流程：
!   Phase 0: 系统初始化（SPICE + 配置 + 引力网 + DACE）
!   Phase 1: 设置测试状态（J2000 惯性系位置/速度）
!   Phase 2: DA vs Real 一致性测试（da_compute_acceleration vs compute_acceleration）
!   Phase 3: DA 内存泄露测试（循环调用 + active_da_count 检测）
!   Phase 4: 清理与结果汇总
! ==============================================================================

program test_da_force_model_consistency
    use pod_global, only: DP
    use pod_engine_module, only: pod_engine_init
    use pod_config, only: config
    use pod_spice, only: str2et
    use pod_force_model_module, only: compute_acceleration
    use pod_da_force_model_module, only: da_compute_acceleration, &
        set_propagation_epoch, da_var  ! <--- 加上这两个
    use pod_dace_classes
    implicit none

    ! ---- 测试参数 ----
    character(len=*), parameter :: CONFIG_FILE = 'dummy_test_config.txt'
    character(len=*), parameter :: TEST_EPOCH  = '2024-03-09T12:00:00'
    integer,          parameter :: DA_ORDER    = 2
    integer,          parameter :: DA_NVARS    = 6

    ! ---- 测试状态（J2000 惯性系） ----
    real(DP), dimension(3) :: pos_real, vel_real
    real(DP) :: tdb_epoch

    ! ---- 实数版结果 ----
    real(DP), dimension(3) :: acc_real

    ! ---- DA 版结果 ----
    type(AlgebraicVector) :: pos_da, vel_da, acc_da
    real(DP), dimension(3) :: acc_da_cons

    ! ---- DA 临时变量 ----
    integer :: i

    ! ---- 内存泄露检测 ----
    integer :: count_before, count_after, count_final

    ! ---- 测试统计 ----
    integer :: n_pass, n_fail
    real(DP), parameter :: TOL = 1.0e-10_DP

    n_pass = 0
    n_fail = 0

    write(*,*) '============================================================'
    write(*,*) '  POD Force Model — DA vs Real 一致性 & 内存泄露测试'
    write(*,*) '============================================================'

    ! ============================================================
    ! Phase 0: 系统初始化
    ! ============================================================
    write(*,*) ''
    write(*,*) '>>> [Phase 0] 系统初始化...'

    ! 0.1 一键初始化（SPICE + 配置 + 实数版引力网）
    call pod_engine_init(CONFIG_FILE)

    ! 0.4 初始化 DACE 引擎
    call dace_initialize(DA_ORDER, DA_NVARS)

    ! 0.5 获取测试历元
    call str2et(TEST_EPOCH, tdb_epoch)
    call set_propagation_epoch(tdb_epoch)
    write(*,*) '  测试历元 (TDB): ', tdb_epoch, ' 秒'

    ! ============================================================
    ! Phase 1: 设置测试状态（J2000 惯性系）
    ! ============================================================
    write(*,*) ''
    write(*,*) '>>> [Phase 1] 设置测试状态...'

    ! 选择一个典型的地月空间轨道状态
    pos_real = [100000.0_DP, 50000.0_DP, 20000.0_DP]  ! km
    vel_real = [1.5_DP, 2.5_DP, 0.5_DP]               ! km/s

    ! ---- DA 版：用常数 DA 赋值（导数 = 0） ----
    if (pos_da%size /= 3) call pos_da%init(3)
    if (vel_da%size /= 3) call vel_da%init(3)
    if (acc_da%size /= 3) call acc_da%init(3)
    do i = 1, 3
        pos_da%elements(i) = pos_real(i) + da_var(i)
        vel_da%elements(i) = vel_real(i) + da_var(i+3)
    end do

    write(*,*) '  位置 (km): ', pos_real
    write(*,*) '  速度 (km/s): ', vel_real

    ! ============================================================
    ! Phase 2: DA vs Real 一致性测试
    ! ============================================================
    write(*,*) ''
    write(*,*) '>>> [Phase 2] DA vs Real 一致性测试...'

    ! ---- 实数版 ----
    call compute_acceleration(pos_real, vel_real, tdb_epoch, acc_real)

    ! ---- DA 版 ----
    call da_compute_acceleration(pos_da, vel_da, tdb_epoch, acc_da)
    acc_da_cons = acc_da%cons()

    ! ---- 比对 ----
    write(*,*) '  实数版 acc (km/s^2): ', acc_real
    write(*,*) '  DA 版   acc (km/s^2): ', acc_da_cons
    do i = 1, 3
        if (abs(acc_real(i) - acc_da_cons(i)) < TOL) then
            n_pass = n_pass + 1
        else
            n_fail = n_fail + 1
            write(*,*) '  *** 不一致! 分量 ', i, &
                       ' 差异 = ', abs(acc_real(i) - acc_da_cons(i))
        end if
    end do

    ! ============================================================
    ! Phase 3: DA 内存泄露测试
    ! ============================================================
    write(*,*) ''
    write(*,*) '>>> [Phase 3] DA 内存泄露测试...'

    ! ---- 3a: da_compute_acceleration 循环测试 ----
    count_before = active_da_count()
    do i = 1, 50
        call da_compute_acceleration(pos_da, vel_da, tdb_epoch, acc_da)
    end do
    count_after = active_da_count()
    write(*,*) '  da_compute_acceleration 循环 50 次后句柄变化: ', &
               count_after - count_before

    if (count_after == count_before) then
        write(*,*) '  ✓ 无内存泄露（句柄数未增长）'
        n_pass = n_pass + 1
    else
        write(*,*) '  *** 内存泄露! 句柄增加了 ', count_after - count_before
        n_fail = n_fail + 1
    end if

    ! ---- 3b: 清理后检查 ----
    call pos_da%destroy()
    call vel_da%destroy()
    call acc_da%destroy()

    ! call earth_grav%dr_da%destroy()  ! 释放 3 个
    ! call moon_grav%dr_da%destroy()   ! 释放 3 个

    count_final = active_da_count()
    write(*,*) '  最终活动句柄数: ', count_final
    if (count_final == 0) then
        write(*,*) '  ✓ 所有 DA 句柄已释放'
        n_pass = n_pass + 1
    else
        write(*,*) '  *** 仍有 ', count_final, ' 个句柄未释放'
        n_fail = n_fail + 1
    end if

    ! ============================================================
    ! Phase 4: 测试结果汇总
    ! ============================================================
    write(*,*) ''
    write(*,*) '============================================================'
    write(*,*) '  测试结果汇总'
    write(*,*) '============================================================'
    write(*,*) '  通过: ', n_pass, ' / ', n_pass + n_fail
    write(*,*) '  失败: ', n_fail

    if (n_fail == 0) then
        write(*,*) ''
        write(*,*) '  ✓ 全部测试通过!'
        write(*,*) '  ✓ DA 版力模型与实数版计算结果一致'
        write(*,*) '  ✓ DA 版力模型无内存泄露'
    else
        write(*,*) ''
        write(*,*) '  *** 存在测试失败，请检查上述输出 ***'
        stop 1
    end if

    write(*,*) '============================================================'

end program test_da_force_model_consistency
