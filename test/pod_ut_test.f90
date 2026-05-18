!> @file pod_ut_test.f90
!> @brief ut 轨道改进集成测试入口
program pod_ut_test
    use pod_global, only: DP, MAX_STRING_LEN
    use pod_engine_module, only: pod_engine_init
    ! 仅引入 Runner 暴露的 API
    use pod_ut_runner_module, only: run_ut_orbit_determination

    implicit none

    ! 命令行与配置参数
    character(len=MAX_STRING_LEN) :: config_file
    character(len=MAX_STRING_LEN) :: obs_file, initial_json_file, output_file_name, site_json_file
    character(len=MAX_STRING_LEN) :: ref_orbit_file, output_residual_file, output_error_file
    character(len=32) :: arg_str
    integer :: i, num_args


    ! ===================================================================
    ! 0. 全局物理环境初始化 (最先执行！)
    ! ===================================================================
    config_file = 'config/dummy_test_config.txt'
    write(*,*) '>>> 正在初始化 CAT POD 物理引擎与星历环境...'
    call pod_engine_init(trim(config_file))
    write(*,*) '>>> 物理引擎初始化完成！'
    
    ! 1. 默认文件路径
    obs_file          = 'OBS/DRO/DRO_single_R91_1h.obs'
    site_json_file    = 'config/site.json'
    initial_json_file = 'OPM/DRO/DRO_init.opm.json'
    ref_orbit_file       = 'ORBITS_REF/DRO/DRO_single_R91_1h.ref'
    output_file_name     = 'OPM/DRO/DRO_single_R91_1h_ut'
    output_residual_file = 'OPM/DRO/DRO_single_R91_1h_ut'
    output_error_file    = 'OPM/DRO/DRO_single_R91_1h_ut'

    ! 3. 打印配置并启动集成 API
    write(*,*) '=================================================='
    write(*,*) '        CAT POD ut-N 轨道改进测试终端          '
    write(*,*) '=================================================='
    write(*,*) '观测数据输入 : ', trim(obs_file)
    write(*,*) '参考轨道输入 : ', trim(ref_orbit_file)
    write(*,*) '定轨结果输出 : ', trim(output_file_name)
    write(*,*) '残差输出     : ', trim(output_residual_file)
    write(*,*) '误差输出     : ', trim(output_error_file)
    write(*,*) '--------------------------------------------------'

    call run_ut_orbit_determination( &
        obs_file             = obs_file, &
        site_json_file       = site_json_file, &
        ref_orbit_file       = ref_orbit_file, &
        initial_json_file    = initial_json_file, &
        output_opm_file      = output_file_name, &
        output_residual_file = output_residual_file, &
        output_error_file    = output_error_file) 
                                       
    write(*,*) '✅ 测试任务圆满完成！'
    write(*,*) '=================================================='

end program pod_ut_test