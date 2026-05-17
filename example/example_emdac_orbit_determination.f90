! ==============================================================================
! example_emdac_orbit_determination.f90
!
! 功能：EMDAC-N (Extended-Multi-Derivative-Algebra-Conditioned) 轨道确定示例
!
! 本示例演示如何使用 POD 框架进行基于微分代数 (DA) 的高阶非线性定轨。
! EMDAC-N 滤波器结合了 UKF 的粒子思想与 DA 的高阶泰勒展开能力，
! 能够高效处理非线性误差传播与复杂动力学环境下的高精度定轨任务。
!
! 运行前请确保：
!   1. 已正确配置 SPICE 星历内核路径 (config/pod_config.txt)
!   2. 观测数据文件与测站配置文件存在
!   3. DACE 库已正确链接
!
! 编译与运行 (使用 fpm):
!   fpm run --example example_emdac_orbit_determination
! ==============================================================================

program example_emdac_orbit_determination
    use pod_global, only: DP, MAX_STRING_LEN
    use pod_engine_module, only: pod_engine_init
    use pod_emdac_runner_module, only: run_emdac_orbit_determination

    implicit none

    ! ---- 文件路径配置 ----
    character(len=MAX_STRING_LEN) :: config_file
    character(len=MAX_STRING_LEN) :: obs_file, site_json_file, initial_json_file, output_file_name
    character(len=MAX_STRING_LEN) :: ref_orbit_file, output_residual_file, output_error_file

    ! ---- 算法参数 ----
    integer  :: opt_particles  = 100000   ! 粒子总数（用于 GMM 近似）
    integer  :: opt_da_order   = 4        ! DA 泰勒展开阶数
    integer  :: opt_em_max_iter = 50      ! EM 算法最大迭代次数
    real(DP) :: opt_em_tol     = 1.0e-4_DP  ! EM 收敛容差
    integer  :: n_components   = 3        ! 高斯混合模型 (GMM) 分量数
    logical  :: gmm_in_switch  = .false.  ! 是否启用 GMM 初始化

    ! ===================================================================
    ! 1. 全局物理环境初始化（最先执行！）
    ! ===================================================================
    config_file = 'config/dummy_test_config.txt'
    write(*, *) '>>> 正在初始化 POD 物理引擎与星历环境...'
    call pod_engine_init(trim(config_file))
    write(*, *) '>>> 物理引擎初始化完成！'

    ! ===================================================================
    ! 2. 配置输入输出文件路径
    ! ===================================================================
    obs_file            = 'OBS/DRO/DRO_single_R91_1h.obs'
    site_json_file      = 'config/site.json'
    initial_json_file   = 'OPM/DRO/DRO_init.opm.json'
    ref_orbit_file      = 'ORBITS_REF/DRO/DRO_single_R91_1h.ref'
    output_file_name    = 'OPM/DRO/DRO_single_R91_1h_emdac'
    output_residual_file = 'OPM/DRO/DRO_single_R91_1h_emdac'
    output_error_file   = 'OPM/DRO/DRO_single_R91_1h_emdac'

    ! ===================================================================
    ! 3. 打印配置信息
    ! ===================================================================
    write(*, *) '=================================================='
    write(*, *) '   EMDAC-N 轨道确定示例'
    write(*, *) '=================================================='
    write(*, *) '粒子总数      : ', opt_particles
    write(*, *) 'GMM 分量数    : ', n_components
    write(*, *) 'DA 阶数       : ', opt_da_order
    write(*, *) 'EM 最大迭代   : ', opt_em_max_iter
    write(*, *) 'EM 收敛容差   : ', opt_em_tol
    write(*, *) 'GMM 初始化    : ', gmm_in_switch
    write(*, *) '--------------------------------------------------'
    write(*, *) '观测数据      : ', trim(obs_file)
    write(*, *) '参考轨道      : ', trim(ref_orbit_file)
    write(*, *) '测站配置      : ', trim(site_json_file)
    write(*, *) '初始轨道      : ', trim(initial_json_file)
    write(*, *) '结果输出      : ', trim(output_file_name)
    write(*, *) '残差输出      : ', trim(output_residual_file)
    write(*, *) '误差输出      : ', trim(output_error_file)
    write(*, *) '--------------------------------------------------'

    ! ===================================================================
    ! 4. 执行 EMDAC-N 轨道确定
    ! ===================================================================
    call run_emdac_orbit_determination( &
        obs_file             = obs_file, &
        site_json_file       = site_json_file, &
        ref_orbit_file       = ref_orbit_file, &
        initial_json_file    = initial_json_file, &
        output_opm_file      = output_file_name, &
        output_residual_file = output_residual_file, &
        output_error_file    = output_error_file, &
        gmm_in_switch        = gmm_in_switch, &
        opt_particles        = opt_particles, &
        max_da_order         = opt_da_order, &
        opt_em_max_iter      = opt_em_max_iter, &
        opt_em_tol           = opt_em_tol, &
        n_components         = n_components)

    write(*, *) '=================================================='
    write(*, *) '   EMDAC-N 轨道确定完成！'
    write(*, *) '   结果已保存至: ', trim(output_file_name)
    write(*, *) '=================================================='

end program example_emdac_orbit_determination
