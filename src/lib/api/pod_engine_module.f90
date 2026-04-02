!> POD System 核心引擎入口模块
module pod_engine_module
    use pod_spice, only: spice_init
    use pod_config, only: load_config, resolve_config_dependencies, print_config
    use pod_force_model_module, only: init_gravity_network
    
    implicit none
    
contains

    !> 系统一键点火接口
    subroutine pod_engine_init(config_file)
        character(len=*), intent(in) :: config_file
        
        write(*,*) ">>> 正在启动 POD 计算引擎..."
        
        ! 1. 唤醒 SPICE 星历内核 (时间与坐标系基准)
        call spice_init()
        
        ! 2. 注入外部配置文件
        call load_config(config_file)
        
        ! 3. 结算衍生参数与依赖检查
        call resolve_config_dependencies()
        
        ! 4. 构建多体引力网络图 (物理场基准)
        call init_gravity_network()

        call print_config()
        
        write(*,*) ">>> 引擎点火完成，已挂载配置: ", trim(config_file)
        write(*,*) "--------------------------------------------------"
        
    end subroutine pod_engine_init

end module pod_engine_module