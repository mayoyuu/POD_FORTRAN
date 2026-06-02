!> @file run_pod_ut.f90
!> @brief CAT POD 系统 UT 轨道改进命令行应用入口 (支持批量调度)
! nohup fpm run run_pod_ut -- \
!   -obs OBS/L1Halo-1/L1Halo-1_supp_single_R91_1h.obs \
!   -init OPM/L1Halo-1/L1Halo-1_init.opm.json \
!   -ref ORBITS_REF/L1Halo-1/L1Halo-1_supp_single_R91_1h.ref \
!   -out OPM/L1Halo-1/0_noise_Q_L1Halo-1_supp_single_R91_1h_ut \
!   -res OPM/L1Halo-1/0_noise_Q_L1Halo-1_supp_single_R91_1h_ut \
!   -err OPM/L1Halo-1/0_noise_Q_L1Halo-1_supp_single_R91_1h_ut \
! !   -p 100000 -o 4 -gmm > ut_single.log 2>&1 &
program run_pod_ut
    use pod_global, only: DP, MAX_STRING_LEN
    use pod_engine_module, only: pod_engine_init
    use pod_ut_runner_module, only: run_ut_orbit_determination

    implicit none

    ! 命令行与配置参数字符串
    character(len=MAX_STRING_LEN) :: config_file
    character(len=MAX_STRING_LEN) :: obs_file, initial_json_file, output_file_name, site_json_file
    character(len=MAX_STRING_LEN) :: ref_orbit_file, output_residual_file, output_error_file
    character(len=MAX_STRING_LEN) :: arg_str
    integer :: i, num_args
    integer :: ext_pos

    ! ===================================================================
    ! 1. 初始化路径变量 (解除硬编码)
    ! ===================================================================
    config_file          = 'config/ut_config.txt'  ! 默认引擎配置文件
    site_json_file       = 'config/site.json'    ! 默认地面站配置文件
    obs_file             = ''
    initial_json_file    = ''
    ref_orbit_file       = ''
    output_file_name     = ''
    output_residual_file = ''
    output_error_file    = ''

    ! ===================================================================
    ! 2. 灵活解析命令行参数
    ! ===================================================================
    num_args = command_argument_count()
    i = 1
    do while (i <= num_args)
        call get_command_argument(i, arg_str)
        select case (trim(arg_str))

            ! --- 核心文件输入/输出接口 ---
            case ('-obs', '--observation')
                call get_command_argument(i+1, arg_str)
                obs_file = trim(arg_str)
                i = i + 1
            case ('-init', '--initial')
                call get_command_argument(i+1, arg_str)
                initial_json_file = trim(arg_str)
                i = i + 1
            case ('-site', '--site')
                call get_command_argument(i+1, arg_str)
                site_json_file = trim(arg_str)
                i = i + 1
            case ('-out', '--output')
                call get_command_argument(i+1, arg_str)
                output_file_name = trim(arg_str)
                i = i + 1
            case ('-ref', '--ref_orbit')
                call get_command_argument(i+1, arg_str)
                ref_orbit_file = trim(arg_str)
                i = i + 1
            case ('-res', '--residual')
                call get_command_argument(i+1, arg_str)
                output_residual_file = trim(arg_str)
                i = i + 1
            case ('-err', '--error')
                call get_command_argument(i+1, arg_str)
                output_error_file = trim(arg_str)
                i = i + 1
            case ('-cfg', '--config')
                call get_command_argument(i+1, arg_str)
                config_file = trim(arg_str)
                i = i + 1

            case default
                write(*,*) '⚠️ 警告: 忽略未知的命令行参数: ', trim(arg_str)
        end select
        i = i + 1
    end do

    ! ===================================================================
    ! 3. 防御性校验与路径自动补全 (确保批量运行时行为安全)
    ! ===================================================================
    if (len_trim(obs_file) == 0) then
        write(*,*) '❌ 错误: 缺漏必需参数 -obs <观测文件路径>'
        stop 1
    end if

    if (len_trim(initial_json_file) == 0) then
        write(*,*) '❌ 错误: 缺漏必需参数 -init <先验状态JSON路径>'
        stop 1
    end if

    ! 如果批量脚本没有指定具体的输出文件名，基于观测文件名自动剔除后缀并生成
    if (len_trim(output_file_name) == 0) then
        ext_pos = index(obs_file, '.', back=.true.)
        if (ext_pos > 1) then
            output_file_name = obs_file(1:ext_pos-1) // '_ut'
        else
            output_file_name = trim(obs_file) // '_ut'
        end if
    end if

    ! 若残差文件路径未显式指定，默认与主输出文件名保持一致
    if (len_trim(output_residual_file) == 0) then
        output_residual_file = trim(output_file_name)
    end if

    ! 若误差文件路径未显式指定且参考轨道已指定，默认与主输出文件名保持一致
    if (len_trim(output_error_file) == 0 .and. len_trim(ref_orbit_file) > 0) then
        output_error_file = trim(output_file_name)
    end if

    ! ===================================================================
    ! 4. 全局物理环境初始化 (最先执行)
    ! ===================================================================
    write(*,*) '>>> 正在初始化 CAT POD 物理引擎与星历环境...'
    write(*,*) '    加载配置文件: ', trim(config_file)
    call pod_engine_init(trim(config_file))
    write(*,*) '>>> 物理引擎初始化完成！'

    ! ===================================================================
    ! 5. 打印任务配置快照
    ! ===================================================================
    write(*,*) '=================================================='
    write(*,*) '         CAT POD UT 轨道改进 CLI 工具              '
    write(*,*) '=================================================='
    write(*,*) '观测数据输入 : ', trim(obs_file)
    write(*,*) '先验状态输入 : ', trim(initial_json_file)
    write(*,*) '测站配置输入 : ', trim(site_json_file)
    if (len_trim(ref_orbit_file) > 0) then
        write(*,*) '参考轨道输入 : ', trim(ref_orbit_file)
    else
        write(*,*) '参考轨道输入 : [未指定]'
    end if
    write(*,*) '定轨结果输出 : ', trim(output_file_name)
    write(*,*) '残差输出     : ', trim(output_residual_file)
    if (len_trim(output_error_file) > 0) then
        write(*,*) '误差输出     : ', trim(output_error_file)
    else
        write(*,*) '误差输出     : [未指定]'
    end if
    write(*,*) '--------------------------------------------------'

    ! ===================================================================
    ! 6. 调用 UT 核心算法接口
    ! ===================================================================
    write(*,*) '>>> 正在启动 UT 无迹卡尔曼滤波计算...'

    if (len_trim(ref_orbit_file) > 0) then
        call run_ut_orbit_determination( &
            obs_file             = obs_file, &
            site_json_file       = site_json_file, &
            ref_orbit_file       = ref_orbit_file, &
            initial_json_file    = initial_json_file, &
            output_opm_file      = output_file_name, &
            output_residual_file = output_residual_file, &
            output_error_file    = output_error_file)
    else
        call run_ut_orbit_determination( &
            obs_file             = obs_file, &
            site_json_file       = site_json_file, &
            initial_json_file    = initial_json_file, &
            output_opm_file      = output_file_name, &
            output_residual_file = output_residual_file)
    end if

    write(*,*) '✅ 该组定轨任务处理完成！'
    write(*,*) '=================================================='

end program run_pod_ut
