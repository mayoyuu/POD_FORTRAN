! ==============================================================================
! test_measurement_da_consistency.f90
!
! 测试目标：
!   1. 验证 DA 版测量模型（compute_measurement_da）的常数项与实数版
!      （compute_measurement）计算结果一致（光学 RA/DEC + 雷达 Range/Az/El）。
!   2. 验证 DA 版测量模型在多次调用后无内存泄露。
!
! 测试流程：
!   Phase 0: 系统初始化（SPICE + 配置 + DACE + 测站）
!   Phase 1: 设置测试状态（J2000 惯性系位置/速度）
!   Phase 2: 光学测量一致性测试（OPTICAL: RA/DEC）
!   Phase 3: 雷达测量一致性测试（RADAR: Range/Az/El）
!   Phase 4: DA 内存泄露测试（循环调用 + active_da_count 检测）
!   Phase 5: 清理与结果汇总
! ==============================================================================

program test_measurement_da_consistency
    use pod_global, only: DP
    use pod_engine_module, only: pod_engine_init
    use pod_spice, only: str2et
    use pod_measurement_base_module, only: observation_station
    use pod_measurement_model_module, only: set_station_from_geodetic, compute_measurement
    use pod_measurement_da_module, only: compute_measurement_da
    use pod_dace_classes
    implicit none

    ! ---- 测试参数 ----
    character(len=*), parameter :: CONFIG_FILE = 'dummy_test_config.txt'
    character(len=*), parameter :: TEST_EPOCH  = '2024-03-09T12:00:00'
    integer,          parameter :: DA_ORDER    = 2
    integer,          parameter :: DA_NVARS    = 6

    ! ---- 测试状态（J2000 惯性系） ----
    real(DP), dimension(3) :: pos_real, vel_real
    real(DP), dimension(6) :: state_real
    real(DP) :: tdb_epoch

    ! ---- 测站 ----
    type(observation_station) :: station

    ! ---- 实数版结果 ----
    real(DP), dimension(2) :: optical_real
    real(DP), dimension(3) :: radar_real

    ! ---- DA 版结果 ----
    type(AlgebraicVector) :: state_da, optical_da, radar_da
    real(DP), dimension(2) :: optical_da_cons
    real(DP), dimension(3) :: radar_da_cons

    ! ---- DA 临时变量 ----
    integer :: i

    ! ---- 内存泄露检测 ----
    integer :: count_before, count_after, count_final

    ! ---- 测试统计 ----
    integer :: n_pass, n_fail
    real(DP), parameter :: TOL = 1.0e-12_DP

    n_pass = 0
    n_fail = 0

    write(*,*) '============================================================'
    write(*,*) '  POD Measurement Model — DA vs Real 一致性 & 内存泄露测试'
    write(*,*) '============================================================'

    ! ============================================================
    ! Phase 0: 系统初始化
    ! ============================================================
    write(*,*) ''
    write(*,*) '>>> [Phase 0] 系统初始化...'

    ! 0.1 一键初始化（SPICE + 配置）
    call pod_engine_init(CONFIG_FILE)

    ! 0.2 初始化 DACE 引擎
    call dace_initialize(DA_ORDER, DA_NVARS)

    ! 0.3 获取测试历元
    call str2et(TEST_EPOCH, tdb_epoch)
    write(*,*) '  测试历元 (TDB): ', tdb_epoch, ' 秒'

    ! 0.4 创建测试测站（北京站，光学+雷达）
    call set_station_from_geodetic(station, 'BEIJING', &
                                   40.0_DP, 116.0_DP, 0.1_DP, 'BOTH')
    write(*,*) '  测站: ', trim(station%name)
    write(*,*) '  测站 ECEF (km): ', station%ecef_position

    ! ============================================================
    ! Phase 1: 设置测试状态（J2000 惯性系）
    ! ============================================================
    write(*,*) ''
    write(*,*) '>>> [Phase 1] 设置测试状态...'

    ! 选择一个典型的地球轨道位置（J2000 惯性系）
    pos_real = [6678.137_DP, 0.0_DP, 0.0_DP]   ! km，近地轨道 ~300km
    vel_real = [0.0_DP, 7.5_DP, 0.0_DP]         ! km/s，圆轨道速度

    state_real(1:3) = pos_real
    state_real(4:6) = vel_real

    ! ---- DA 版：用常数 DA 赋值（导数 = 0） ----
    if (state_da%size /= 6) call state_da%init(6)
    do i = 1, 3
        call state_da%set(i, pos_real(i) + da_var(i))
        call state_da%set(i+3, vel_real(i) + da_var(i+3))
    end do

    write(*,*) '  位置 (km): ', pos_real
    write(*,*) '  速度 (km/s): ', vel_real

    ! ============================================================
    ! Phase 2: 光学测量一致性测试（OPTICAL: RA/DEC）
    ! ============================================================
    write(*,*) ''
    write(*,*) '>>> [Phase 2] 光学测量一致性测试...'

    ! ---- 实数版 ----
    call compute_measurement(state_real, tdb_epoch, station, 'OPTICAL', optical_real)

    ! ---- DA 版 ----
    call compute_measurement_da(state_da, tdb_epoch, station, 'OPTICAL', optical_da)
    optical_da_cons(1) = optical_da%elements(1)%cons()
    optical_da_cons(2) = optical_da%elements(2)%cons()

    ! ---- 比对 ----
    write(*,*) '  实数版 RA (rad): ', optical_real(1)
    write(*,*) '  DA 版   RA (rad): ', optical_da_cons(1)
    write(*,*) '  实数版 DEC (rad): ', optical_real(2)
    write(*,*) '  DA 版   DEC (rad): ', optical_da_cons(2)
    do i = 1, 2
        if (abs(optical_real(i) - optical_da_cons(i)) < TOL) then
            n_pass = n_pass + 1
        else
            n_fail = n_fail + 1
            write(*,*) '  *** 光学不一致! 分量 ', i, &
                       ' 差异 = ', abs(optical_real(i) - optical_da_cons(i))
        end if
    end do

    ! ============================================================
    ! Phase 3: 雷达测量一致性测试（RADAR: Range/Az/El）
    ! ============================================================
    write(*,*) ''
    write(*,*) '>>> [Phase 3] 雷达测量一致性测试...'

    ! ---- 实数版 ----
    call compute_measurement(state_real, tdb_epoch, station, 'RADAR', radar_real)

    ! ---- DA 版 ----
    call compute_measurement_da(state_da, tdb_epoch, station, 'RADAR', radar_da)
    radar_da_cons(1) = radar_da%elements(1)%cons()
    radar_da_cons(2) = radar_da%elements(2)%cons()
    radar_da_cons(3) = radar_da%elements(3)%cons()

    ! ---- 比对 ----
    write(*,*) '  实数版 Range (km): ', radar_real(1)
    write(*,*) '  DA 版   Range (km): ', radar_da_cons(1)
    write(*,*) '  实数版 Azimuth (rad): ', radar_real(2)
    write(*,*) '  DA 版   Azimuth (rad): ', radar_da_cons(2)
    write(*,*) '  实数版 Elevation (rad): ', radar_real(3)
    write(*,*) '  DA 版   Elevation (rad): ', radar_da_cons(3)
    do i = 1, 3
        if (abs(radar_real(i) - radar_da_cons(i)) < TOL) then
            n_pass = n_pass + 1
        else
            n_fail = n_fail + 1
            write(*,*) '  *** 雷达不一致! 分量 ', i, &
                       ' 差异 = ', abs(radar_real(i) - radar_da_cons(i))
        end if
    end do

    ! ============================================================
    ! Phase 4: DA 内存泄露测试
    ! ============================================================
    write(*,*) ''
    write(*,*) '>>> [Phase 4] DA 内存泄露测试...'

    ! ---- 4a: compute_measurement_da 循环测试（光学） ----
    count_before = active_da_count()
    do i = 1, 50
        call compute_measurement_da(state_da, tdb_epoch, station, 'OPTICAL', optical_da)
        call compute_measurement_da(state_da, tdb_epoch, station, 'RADAR', radar_da)
    end do
    count_after = active_da_count()
    write(*,*) '  compute_measurement_da 循环 50 次后句柄变化: ', &
               count_after - count_before

    if (count_after == count_before) then
        write(*,*) '  ✓ 无内存泄露（句柄数未增长）'
        n_pass = n_pass + 1
    else
        write(*,*) '  *** 内存泄露! 句柄增加了 ', count_after - count_before
        n_fail = n_fail + 1
    end if

    ! ---- 4b: 清理后检查 ----
    call state_da%destroy()
    call optical_da%destroy()
    call radar_da%destroy()

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
        write(*,*) '  ✓ DA 版测量模型与实数版计算结果一致'
        write(*,*) '  ✓ DA 版测量模型无内存泄露'
    else
        write(*,*) ''
        write(*,*) '  *** 存在测试失败，请检查上述输出 ***'
        stop 1
    end if

    write(*,*) '============================================================'

end program test_measurement_da_consistency
