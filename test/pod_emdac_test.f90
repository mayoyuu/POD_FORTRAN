!> @file pod_emdac_test.f90
!> @brief EMDAC-N 轨道改进集成测试入口
program pod_emdac_test
    use pod_global, only: DP, MAX_STRING_LEN
    use pod_engine_module, only: pod_engine_init
    ! 仅引入 Runner 暴露的 API
    use pod_emdac_runner_module, only: run_emdac_orbit_determination

    implicit none

    ! 命令行与配置参数
    character(len=MAX_STRING_LEN) :: config_file
    character(len=MAX_STRING_LEN) :: obs_file, initial_json_file, output_json_file, site_json_file
    character(len=32) :: arg_str
    integer :: i, num_args
    
    ! 算法控制参数 (赋予默认值)
    integer :: opt_particles = 100000
    integer :: opt_da_order = 2
    integer :: opt_em_max_iter = 50
    real(DP) :: opt_em_tol = 1.0e-4_DP
    integer :: n_components = 5
    logical :: gmm_in_switch = .false.  ! 补充了 GMM 初始化开关的默认值

    ! ===================================================================
    ! 0. 全局物理环境初始化 (最先执行！)
    ! ===================================================================
    config_file = 'config/pod_config.txt'
    write(*,*) '>>> 正在初始化 CAT POD 物理引擎与星历环境...'
    call pod_engine_init(trim(config_file))
    write(*,*) '>>> 物理引擎初始化完成！'
    
    ! 1. 默认文件路径
    obs_file          = 'data/DROB_20251210_20260111_cor.obs'
    site_json_file    = 'config/site.json'
    initial_json_file = 'data/DROb_20251210_9.opm'
    output_json_file  = 'output/emdac_result.opm'
    
    ! 2. 灵活解析命令行可选参数
    num_args = command_argument_count()
    i = 1
    do while (i <= num_args)
        call get_command_argument(i, arg_str)
        select case (trim(arg_str))
            case ('-p', '--particles')
                call get_command_argument(i+1, arg_str); read(arg_str, *) opt_particles
                i = i + 1
            case ('-o', '--order')
                call get_command_argument(i+1, arg_str); read(arg_str, *) opt_da_order
                i = i + 1
            case ('-tol', '--em_tol')
                call get_command_argument(i+1, arg_str); read(arg_str, *) opt_em_tol
                i = i + 1
            case ('-iter', '--em_iter')
                call get_command_argument(i+1, arg_str); read(arg_str, *) opt_em_max_iter
                i = i + 1
            case ('-n', '--ncomp')
                call get_command_argument(i+1, arg_str); read(arg_str, *) n_components
                i = i + 1
            case ('-gmm', '--use_gmm')
                ! 增加一个无参开关，碰到这个参数就开启 GMM 初始化
                gmm_in_switch = .true.
        end select
        i = i + 1
    end do
    
    ! 3. 打印配置并启动集成 API
    write(*,*) '=================================================='
    write(*,*) '        CAT POD EMDAC-N 轨道改进测试终端          '
    write(*,*) '=================================================='
    write(*,*) '粒子总数     : ', opt_particles
    write(*,*) 'GMM 分量数   : ', n_components
    write(*,*) 'DA 阶数      : ', opt_da_order
    write(*,*) 'EM 最大迭代  : ', opt_em_max_iter
    write(*,*) 'EM 收敛容差  : ', opt_em_tol
    write(*,*) 'GMM 初始化   : ', gmm_in_switch
    write(*,*) '--------------------------------------------------'
    write(*,*) '观测数据输入 : ', trim(obs_file)
    write(*,*) '定轨结果输出 : ', trim(output_json_file)
    write(*,*) '--------------------------------------------------'
    
    ! 直接调用核心 API 进行“黑盒”执行，采用极其安全的【全关键字传参】模式
    call run_emdac_orbit_determination( &
        obs_file          = obs_file, &
        site_json_file    = site_json_file, &
        gmm_in_switch     = gmm_in_switch, &
        initial_json_file = initial_json_file, &
        output_json_file  = output_json_file, &
        opt_particles     = opt_particles, &
        opt_da_order      = opt_da_order, &
        opt_em_max_iter   = opt_em_max_iter, &
        opt_em_tol        = opt_em_tol, &
        n_components      = n_components)
                                       
    write(*,*) '✅ 测试任务圆满完成！'
    write(*,*) '=================================================='

end program pod_emdac_test