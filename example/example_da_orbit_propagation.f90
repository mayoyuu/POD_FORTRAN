! ==============================================================================
! example_da_orbit_propagation.f90
!
! 功能：DA (Differential Algebra) 轨道传播示例
!
! 本示例演示如何使用 POD 框架的微分代数 (DA) 轨道传播器。
! DA 技术能够将轨道状态表示为泰勒多项式，在一次传播中同时获得：
!   1. 名义轨道（常数项）
!   2. 状态转移矩阵 (STM) — 状态对初始偏差的偏导数
!   3. 高阶非线性映射信息
!
! 运行前请确保：
!   1. 已正确配置 SPICE 星历内核路径 (config/pod_config.txt)
!   2. DACE 库已正确链接
!
! 编译与运行 (使用 fpm):
!   fpm run --example example_da_orbit_propagation
!
! 或使用 cmake:
!   ./build/example/example_da_orbit_propagation
! ==============================================================================

program example_da_orbit_propagation
    use pod_global, only: DP
    use pod_engine_module, only: pod_engine_init
    use pod_config, only: config
    use pod_spice, only: str2et

    ! ---- DA 版力模型与传播器 ----
    use pod_da_force_model_module, only: &
        set_propagation_epoch_da => set_propagation_epoch, &
        current_epoch0_da => current_epoch0, &
        init_gravity_network
    use pod_da_orbit_propagation, only: &
        da_orbit_state, da_propagation_result, &
        propagate_da_orbit, display_da_propagation_results, &
        save_da_propagation_results, cleanup_da_propagation_result, &
        extract_stm_from_result

    use pod_dace_classes
    use pod_basicmath_module, only: inverse_and_determinant

    implicit none

    ! ---- 示例参数 ----
    character(len=*), parameter :: CONFIG_FILE = 'dummy_test_config.txt'
    character(len=*), parameter :: TEST_EPOCH  = '2025-12-15T00:00:01.000000'
    integer,          parameter :: DA_ORDER    = 3   ! DA 泰勒展开阶数
    integer,          parameter :: DA_NVARS    = 6   ! 状态变量数 (x,y,z,vx,vy,vz)

    ! ---- 物理状态（J2000 惯性系，有量纲） ----
    real(DP), dimension(6) :: state_physical

    ! ---- DA 版传播 ----
    type(da_orbit_state) :: da_initial_state
    type(da_propagation_result) :: da_result

    ! ---- STM 提取 ----
    real(DP), dimension(6, 6) :: stm, stm_inv
    real(DP) :: stm_det
    integer :: info

    ! ---- 杂项 ----
    real(DP) :: tdb_epoch, propagation_time
    integer :: i, j

    write(*, *) '============================================================'
    write(*, *) '   DA 轨道传播示例'
    write(*, *) '============================================================'

    ! ============================================================
    ! 1. 系统初始化
    ! ============================================================
    write(*, *) ''
    write(*, *) '>>> [Step 1] 系统初始化...'

    ! 1.1 一键初始化（加载 SPICE 内核、时间系统等）
    call pod_engine_init(CONFIG_FILE)

    ! 1.2 开启复杂引力环境
    !     地球 10 阶非球形引力 + 月球 10 阶非球形引力 + 太阳第三体
    config%use_earth_nspheric = .true.
    config%earth_degree       = 10
    config%use_moon_nspheric  = .true.
    config%moon_degree        = 10
    config%use_third_body     = .true.
    config%use_planet(10)     = .true.  ! 激活月球
    config%use_planet(11)     = .true.  ! 激活太阳

    ! 1.3 初始化 DA 版引力网络
    call init_gravity_network()

    ! 1.4 初始化 DACE 引擎
    call dace_initialize(DA_ORDER, DA_NVARS)

    ! 1.5 获取并设置测试历元
    call str2et(TEST_EPOCH, tdb_epoch)
    call set_propagation_epoch_da(tdb_epoch)
    write(*, *) '   测试历元 (TDB): ', tdb_epoch, ' 秒'

    ! ============================================================
    ! 2. 设置 DRO 初始状态
    ! ============================================================
    write(*, *) ''
    write(*, *) '>>> [Step 2] 设置初始轨道状态 (DRO 轨道)...'

    ! DRO (Distant Retrograde Orbit) 初始状态
    ! 单位: 位置 (km), 速度 (km/s)
    state_physical = [-402779.291910_DP, -181111.455576_DP, -109249.167033_DP, &  ! 位置
                             0.291707_DP,    -0.504100_DP,    -0.263272_DP]        ! 速度

    ! 构建 DA 初始状态
    da_initial_state%nominal_state = state_physical
    da_initial_state%epoch         = tdb_epoch
    da_initial_state%da_order      = DA_ORDER

    propagation_time = 86400.0_DP * 4  ! 积分 4 天

    write(*, *) '   初始状态 (km, km/s): ', state_physical
    write(*, *) '   积分时长: 4 天 (', propagation_time, ' 秒)'
    write(*, *) '   DA 阶数: ', DA_ORDER

    ! ============================================================
    ! 3. 执行 DA 轨道传播
    ! ============================================================
    write(*, *) ''
    write(*, *) '>>> [Step 3] 执行 DA 轨道传播 (RKF78 积分器)...'

    ! propagate_da_orbit 使用 DA 版力模型和积分器进行轨道传播
    ! 第二个参数: 积分时长 (秒)
    ! 第三个参数: 积分器类型 (1=RKF45, 2=RKF78)
    call propagate_da_orbit(da_initial_state, propagation_time, 2, da_result)

    write(*, *) '   DA 轨道传播完成！'
    write(*, *) '   最终状态 (常数项, km, km/s): ', da_result%final_state%cons()

    ! ============================================================
    ! 4. 提取状态转移矩阵 (STM)
    ! ============================================================
    write(*, *) ''
    write(*, *) '>>> [Step 4] 提取 6×6 状态转移矩阵 (STM)...'

    ! 从 DA 传播结果中提取 STM
    ! STM 描述了最终状态对初始状态偏差的线性映射关系
    call extract_stm_from_result(da_result, stm)

    write(*, *) '   提取的 STM (6×6):'
    do i = 1, 6
        write(*, "(4X, 6(ES15.6, 2X))") (stm(i, j), j = 1, 6)
    end do

    ! 验证 STM 的辛结构（行列式应 ≈ 1.0）
    call inverse_and_determinant(stm, stm_inv, stm_det, info)
    write(*, *) ''
    write(*, *) '   STM 行列式: ', stm_det
    if (info == 0) then
        write(*, *) '   ✓ STM 矩阵可逆，辛结构保持良好'
    end if

    ! ============================================================
    ! 5. 显示与保存结果
    ! ============================================================
    write(*, *) ''
    write(*, *) '>>> [Step 5] 显示与保存传播结果...'

    ! 显示详细的 DA 传播结果
    call display_da_propagation_results(da_result)

    ! 保存结果到文件
    call save_da_propagation_results(da_result, 'output/da_propagation_result.txt')

    ! ============================================================
    ! 6. 清理
    ! ============================================================
    call cleanup_da_propagation_result(da_result)
    write(*, *) ''
    write(*, *) '============================================================'
    write(*, *) '   DA 轨道传播示例完成！'
    write(*, *) '   结果已保存至: output/da_propagation_result.txt'
    write(*, *) '============================================================'

end program example_da_orbit_propagation
