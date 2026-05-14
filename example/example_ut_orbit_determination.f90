! ==============================================================================
! example_ut_orbit_determination.f90
!
! 功能：UT (Unscented Transform) 轨道确定示例
!
! 本示例演示如何使用 POD 框架进行基于无迹变换 (UT) 的轨道确定。
! UT 方法通过确定性采样（Sigma 点）来近似非线性系统的状态分布，
! 相比传统 EKF 具有更高的非线性近似精度，且无需计算雅可比矩阵。
!
! 运行前请确保：
!   1. 已正确配置 SPICE 星历内核路径 (config/pod_config.txt)
!   2. 观测数据文件与测站配置文件存在
!
! 编译与运行 (使用 fpm):
!   fpm run --example example_ut_orbit_determination
!
! 或使用 cmake:
!   ./build/example/example_ut_orbit_determination
! ==============================================================================

program example_ut_orbit_determination
    use pod_global, only: DP, MAX_STRING_LEN
    use pod_engine_module, only: pod_engine_init
    use pod_ut_runner_module, only: run_ut_orbit_determination

    implicit none

    ! ---- 文件路径配置 ----
    character(len=MAX_STRING_LEN) :: config_file
    character(len=MAX_STRING_LEN) :: obs_file
    character(len=MAX_STRING_LEN) :: site_json_file
    character(len=MAX_STRING_LEN) :: initial_json_file
    character(len=MAX_STRING_LEN) :: output_file_name

    ! ===================================================================
    ! 1. 全局物理环境初始化（最先执行！）
    ! ===================================================================
    ! 该步骤会加载 SPICE 内核、初始化时间系统、配置力模型参数等。
    config_file = 'dummy_test_config.txt'
    write(*, *) '>>> 正在初始化 POD 物理引擎与星历环境...'
    call pod_engine_init(trim(config_file))
    write(*, *) '>>> 物理引擎初始化完成！'

    ! ===================================================================
    ! 2. 配置输入输出文件路径
    ! ===================================================================
    ! 请根据实际数据位置修改以下路径
    obs_file          = 'input/DROB_20251210_20260111_cor.obs'  ! 观测数据文件
    site_json_file    = 'config/site.json'                       ! 测站配置文件
    initial_json_file = 'input/DROb_20251210_9.opm'              ! 初始轨道 (OPM 格式)
    output_file_name  = 'output/DROb_202601_ut_result'           ! 定轨结果输出路径

    ! ===================================================================
    ! 3. 打印配置信息
    ! ===================================================================
    write(*, *) '=================================================='
    write(*, *) '   UT 轨道确定示例'
    write(*, *) '=================================================='
    write(*, *) '观测数据      : ', trim(obs_file)
    write(*, *) '测站配置      : ', trim(site_json_file)
    write(*, *) '初始轨道      : ', trim(initial_json_file)
    write(*, *) '结果输出      : ', trim(output_file_name)
    write(*, *) '--------------------------------------------------'

    ! ===================================================================
    ! 4. 执行 UT 轨道确定
    ! ===================================================================
    ! run_ut_orbit_determination 是基于无迹变换的定轨 API，内部封装了：
    !   - 观测数据读取与预处理
    !   - Sigma 点生成与 UT 传播
    !   - 观测更新与状态估计
    !   - 结果输出与后处理
    !
    ! 与 EMDAC-N 相比，UT 方法计算量更小，适合中等非线性程度的定轨场景。
    call run_ut_orbit_determination( &
        obs_file          = obs_file, &
        site_json_file    = site_json_file, &
        initial_json_file = initial_json_file, &
        output_file_name  = output_file_name)

    write(*, *) '=================================================='
    write(*, *) '   UT 轨道确定完成！'
    write(*, *) '   结果已保存至: ', trim(output_file_name)
    write(*, *) '=================================================='

end program example_ut_orbit_determination
