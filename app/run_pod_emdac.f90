!> @file run_pod_emdac.f90
!> @brief CAT POD 系统 EMDAC-N 轨道改进命令行应用入口 (支持批量调度)
! nohup fpm run run_pod_emdac -- \
!   -obs OBS/L1Halo/L1Halo_single_R91_1h_mag20_par.obs \
!   -init OPM/L1Halo/L1Halo_init.opm.json \
!   -ref ORBITS_REF/L1Halo/L1Halo_single_R91_1h_mag20_par.ref \
!   -out OPM/L1Halo/0NoiseR_L1Halo_single_R91_1h_mag20_par_emdac_n3 \
!   -res OPM/L1Halo/0NoiseR_L1Halo_single_R91_1h_mag20_par_emdac_n3 \
!   -err OPM/L1Halo/0NoiseR_L1Halo_single_R91_1h_mag20_par_emdac_n3 \
! !   -p 100000 -o 4 -gmm > emdac_single.log 2>&1 &
program run_pod_emdac
    use pod_global, only: DP, MAX_STRING_LEN
    use pod_engine_module, only: pod_engine_init
    use pod_emdac_runner_module, only: run_emdac_orbit_determination

    implicit none

    ! 命令行与配置参数字符串
    character(len=MAX_STRING_LEN) :: config_file
    character(len=MAX_STRING_LEN) :: obs_file, initial_json_file, output_file_name, site_json_file
    character(len=MAX_STRING_LEN) :: ref_orbit_file, output_residual_file, output_error_file
    character(len=MAX_STRING_LEN) :: arg_str
    integer :: i, num_args
    integer :: ext_pos
    
    ! 算法控制参数 (赋予默认业务值)
    integer :: opt_particles = 100000
    integer :: opt_da_order = 5
    integer :: opt_em_max_iter = 50
    real(DP) :: opt_em_tol = 1.0e-4_DP
    integer :: n_components = 3
    logical :: gmm_in_switch = .false.

    ! ===================================================================
    ! 1. 初始化路径变量 (解除硬编码)
    ! ===================================================================
    config_file          = 'config/config.txt'  ! 默认引擎配置文件
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

            ! --- 滤波器与算法控制接口 ---
            case ('-p', '--particles')
                call get_command_argument(i+1, arg_str)
                read(arg_str, *) opt_particles
                i = i + 1
            case ('-o', '--order')
                call get_command_argument(i+1, arg_str)
                read(arg_str, *) opt_da_order
                i = i + 1
            case ('-tol', '--em_tol')
                call get_command_argument(i+1, arg_str)
                read(arg_str, *) opt_em_tol
                i = i + 1
            case ('-iter', '--em_iter')
                call get_command_argument(i+1, arg_str)
                read(arg_str, *) opt_em_max_iter
                i = i + 1
            case ('-n', '--ncomp')
                call get_command_argument(i+1, arg_str)
                read(arg_str, *) n_components
                i = i + 1
            case ('-gmm', '--use_gmm')
                gmm_in_switch = .true.

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
            output_file_name = obs_file(1:ext_pos-1) // '_emdac'
        else
            output_file_name = trim(obs_file) // '_emdac'
        end if
    end if

    ! 若残差与误差文件路径未显式指定，默认与主输出文件名保持一致
    if (len_trim(output_residual_file) == 0) then
        output_residual_file = trim(output_file_name)
    end if
    if (len_trim(output_error_file) == 0) then
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
    write(*,*) '         CAT POD EMDAC-N 轨道改进 CLI 工具        '
    write(*,*) '=================================================='
    write(*,*) '粒子总数     : ', opt_particles
    write(*,*) 'GMM 分量数   : ', n_components
    write(*,*) 'DA 阶数      : ', opt_da_order
    write(*,*) 'EM 最大迭代  : ', opt_em_max_iter
    write(*,*) 'EM 收敛容差  : ', opt_em_tol
    write(*,*) 'GMM 初始化   : ', gmm_in_switch
    write(*,*) '--------------------------------------------------'
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
    write(*,*) '误差输出     : ', trim(output_error_file)
    write(*,*) '--------------------------------------------------'
    
    ! ===================================================================
    ! 6. 调用 EMDAC-N 核心算法接口
    ! ===================================================================
    write(*,*) '>>> 正在启动 EMDAC-N 高阶非线性滤波计算...'
    
    if (len_trim(ref_orbit_file) > 0) then
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
    else
        call run_emdac_orbit_determination( &
            obs_file             = obs_file, &
            site_json_file       = site_json_file, &
            initial_json_file    = initial_json_file, &
            output_opm_file      = output_file_name, &
            output_residual_file = output_residual_file, &
            gmm_in_switch        = gmm_in_switch, &
            opt_particles        = opt_particles, &
            max_da_order         = opt_da_order, &
            opt_em_max_iter      = opt_em_max_iter, &
            opt_em_tol           = opt_em_tol, &
            n_components         = n_components)
    end if
    
    write(*,*) '✅ 该组定轨任务处理完成！'
    write(*,*) '=================================================='

end program run_pod_emdac