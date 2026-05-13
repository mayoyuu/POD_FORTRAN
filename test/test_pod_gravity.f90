! ==============================================================================
! test_pod_gravity.f90
!
! 测试目标：
!   1. 验证 DA 版引力场（f_zonal_da / f_tesseral_da）的常数项与实数版
!      （f_zonal / f_tesseral）计算结果一致。
!   2. 验证 DA 版引力场在多次调用后无内存泄露。
!
! 测试流程：
!   Phase 0: 系统初始化（SPICE + 配置 + 引力网）
!   Phase 1: 在 J2000 惯性系中设置测试位置，通过 pxform 转到地固系
!   Phase 2: 带谐项一致性测试（f_zonal vs f_zonal_da）
!   Phase 3: 田谐项一致性测试（f_tesseral vs f_tesseral_da）
!   Phase 4: DA 内存泄露测试（循环调用 + active_da_count 检测）
!   Phase 5: 清理
! ==============================================================================

program test_pod_gravity
    use pod_global, only: DP
    use pod_engine_module, only: pod_engine_init
    use pod_config, only: config
    use pod_spice, only: str2et, pxform
    use pod_da_force_model_module, only: earth_grav, init_gravity_network
    use pod_dace_classes
    implicit none

    ! ---- 测试参数 ----
    character(len=*), parameter :: CONFIG_FILE = 'dummy_test_config.txt'
    character(len=*), parameter :: TEST_EPOCH  = '2024-03-09T12:00:00'
    integer,          parameter :: DA_ORDER    = 2
    integer,          parameter :: DA_NVARS    = 3

    ! ---- 位置向量（J2000 惯性系，单位 km） ----
    real(DP), dimension(3) :: pos_j2000

    ! ---- 坐标转换 ----
    real(DP) :: tdb_epoch
    real(DP), dimension(3,3) :: rot_to_body, rot_to_inertial

    ! ---- 实数版结果 ----
    real(DP), dimension(3) :: acc_z_real, acc_t_real, acc_grav_real

    ! ---- DA 版结果 ----
    type(AlgebraicVector) :: acc_z_da, acc_t_da
    real(DP), dimension(3) :: acc_z_da_cons, acc_t_da_cons

    ! ---- DA 临时变量 ----
    type(DA) :: tmp_da
    integer :: i

    ! ---- 内存泄露检测 ----
    integer :: count_before, count_after, count_final

    ! ---- 测试统计 ----
    integer :: n_pass, n_fail
    real(DP), parameter :: TOL = 1.0e-12_DP

    n_pass = 0
    n_fail = 0

    write(*,*) '============================================================'
    write(*,*) '  POD Gravity Model — DA vs Real 一致性 & 内存泄露测试'
    write(*,*) '============================================================'

    ! ============================================================
    ! Phase 0: 系统初始化
    ! ============================================================
    write(*,*) ''
    write(*,*) '>>> [Phase 0] 系统初始化...'

    ! 0.1 一键初始化（SPICE + 配置 + 实数版引力网）
    call pod_engine_init(CONFIG_FILE)

    ! 0.2 启用地球非球形引力（dummy_test_config.txt 中默认为 false）
    config%use_earth_nspheric = .true.
    config%earth_degree       = 10

    ! 0.3 初始化 DA 版引力网（earth_grav 在 pod_da_force_model_module 中）
    call init_gravity_network()

    ! 0.4 初始化 DACE 引擎
    call dace_initialize(DA_ORDER, DA_NVARS)

    ! 0.5 获取测试历元
    call str2et(TEST_EPOCH, tdb_epoch)
    write(*,*) '  测试历元 (TDB): ', tdb_epoch, ' 秒'

    ! ============================================================
    ! Phase 1: 设置测试位置（J2000 → 地固系）
    ! ============================================================
    write(*,*) ''
    write(*,*) '>>> [Phase 1] 设置测试位置...'

    ! 选择一个典型的地球轨道位置（J2000 惯性系）
    pos_j2000 = [6678.137_DP, 0.0_DP, 0.0_DP]  ! 近地轨道 ~300km 高度

    ! 通过 SPICE pxform 转到地固系
    call pxform('J2000', 'IAU_EARTH', tdb_epoch, rot_to_body)
    rot_to_inertial = transpose(rot_to_body)

    ! ---- 实数版：设置 earth_grav%dr（地固系） ----
    earth_grav%dr = matmul(rot_to_body, pos_j2000)

    ! ---- DA 版：设置 earth_grav%dr_da（地固系，常数 DA） ----
    ! 先初始化 dr_da 为 3 维向量
    if (earth_grav%dr_da%size /= 3) call earth_grav%dr_da%init(3)
    do i = 1, 3
        ! 用常数 DA 赋值：常数 = earth_grav%dr(i)，导数 = 0
        earth_grav%dr_da%elements(i) = earth_grav%dr(i)
    end do

    write(*,*) '  地固系位置 (km): ', earth_grav%dr

    ! ============================================================
    ! Phase 2: 带谐项一致性测试（f_zonal vs f_zonal_da）
    ! ============================================================
    write(*,*) ''
    write(*,*) '>>> [Phase 2] 带谐项一致性测试...'

    ! ---- 实数版 ----
    call earth_grav%f_zonal(acc_z_real)

    ! ---- DA 版 ----
    call acc_z_da%init(3)
    call earth_grav%f_zonal_da(acc_z_da)
    do i = 1, 3
        acc_z_da_cons(i) = acc_z_da%elements(i)%cons()
    end do

    ! ---- 比对 ----
    write(*,*) '  实数版 acc_z (km/s^2): ', acc_z_real
    write(*,*) '  DA 版   acc_z (km/s^2): ', acc_z_da_cons
    do i = 1, 3
        if (abs(acc_z_real(i) - acc_z_da_cons(i)) < TOL) then
            n_pass = n_pass + 1
        else
            n_fail = n_fail + 1
            write(*,*) '  *** 不一致! 分量 ', i, &
                       ' 差异 = ', abs(acc_z_real(i) - acc_z_da_cons(i))
        end if
    end do

    ! ============================================================
    ! Phase 3: 田谐项一致性测试（f_tesseral vs f_tesseral_da）
    ! ============================================================
    write(*,*) ''
    write(*,*) '>>> [Phase 3] 田谐项一致性测试...'

    ! ---- 实数版 ----
    call earth_grav%f_tesseral(acc_t_real)

    ! ---- DA 版 ----
    call acc_t_da%init(3)
    call earth_grav%f_tesseral_da(acc_t_da)
    do i = 1, 3
        acc_t_da_cons(i) = acc_t_da%elements(i)%cons()
    end do

    ! ---- 比对 ----
    write(*,*) '  实数版 acc_t (km/s^2): ', acc_t_real
    write(*,*) '  DA 版   acc_t (km/s^2): ', acc_t_da_cons
    do i = 1, 3
        if (abs(acc_t_real(i) - acc_t_da_cons(i)) < TOL) then
            n_pass = n_pass + 1
        else
            n_fail = n_fail + 1
            write(*,*) '  *** 不一致! 分量 ', i, &
                       ' 差异 = ', abs(acc_t_real(i) - acc_t_da_cons(i))
        end if
    end do

    ! ============================================================
    ! Phase 4: DA 内存泄露测试
    ! ============================================================
    write(*,*) ''
    write(*,*) '>>> [Phase 4] DA 内存泄露测试...'

    ! ---- 4a: f_zonal_da 循环测试 ----
    count_before = active_da_count()
    do i = 1, 50
        call earth_grav%f_zonal_da(acc_z_da)
        call earth_grav%f_tesseral_da(acc_t_da)
    end do
    count_after = active_da_count()
    write(*,*) '  f_zonal_da + f_tesseral_da 循环 50 次后句柄变化: ', &
               count_after - count_before

    if (count_after == count_before) then
        write(*,*) '  ✓ 无内存泄露（句柄数未增长）'
        n_pass = n_pass + 1
    else
        write(*,*) '  *** 内存泄露! 句柄增加了 ', count_after - count_before
        n_fail = n_fail + 1
    end if

    ! ---- 4b: 清理后检查 ----
    call acc_z_da%destroy()
    call acc_t_da%destroy()
    call earth_grav%dr_da%destroy()

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
    ! Phase 5: 测试结果汇总
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
        write(*,*) '  ✓ DA 版引力场与实数版计算结果一致'
        write(*,*) '  ✓ DA 版引力场无内存泄露'
    else
        write(*,*) ''
        write(*,*) '  *** 存在测试失败，请检查上述输出 ***'
        stop 1
    end if

    write(*,*) '============================================================'

end program test_pod_gravity
