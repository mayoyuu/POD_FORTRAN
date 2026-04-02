module pod_uq_propagation
    use pod_global, only: DP, MAX_STRING_LEN, output_directory
    use pod_utils,  only: print_separator, print_vector, print_matrix
    
    ! 引入底层计算库
    use pod_uq_base_module,  only: uq_propagator_base
    use pod_uq_da_module,    only: uq_da_propagator
    use pod_uq_mc_module,    only: uq_mc_propagator
    use pod_uq_state_module, only: uq_state_type
    use pod_random_module,   only: init_random_seed, generate_multivariate_normal
    
    implicit none

    ! 传播算法常量定义 (对外暴露作为开关)
    integer, parameter, public :: METHOD_MC = 1
    integer, parameter, public :: METHOD_DA = 2
    
contains

    ! ====================================================================
    ! 核心 API：执行不确定性传播
    ! 调用者只需传入标称状态、协方差矩阵以及控制开关，即可获得传播前后的完整分布状态
    ! ====================================================================
    subroutine run_uq_propagation(nominal_state, initial_cov, epoch0, t_start, t_end,&
                                  method_switch, n_particles, &
                                  save_results_to_file, initial_state_out, final_state_out, &
                                  integrator_switch, file_prefix)
        
        real(DP), intent(in) :: nominal_state(6)   ! 标称轨道 (作为均值)
        real(DP), intent(in) :: initial_cov(6,6)   ! 初始协方差矩阵
        real(DP), intent(in) :: epoch0             ! 【新增】物理历元基准 (TDB 秒)
        real(DP), intent(in) :: t_start, t_end     ! 积分起止相对时间
        integer,  intent(in) :: method_switch      ! 方法开关 (METHOD_MC 或 METHOD_DA)
        ! integer,  intent(in) :: integrator_switch  ! 积分器开关 (INTEG_RK4, INTEG_RKF45 等)
        integer,  intent(in) :: n_particles        ! 生成的粒子数量
        logical,  intent(in) :: save_results_to_file ! 是否将结果落盘
        integer,  intent(in), optional :: integrator_switch ! 可选的积分器开关参数
        character(len=*), intent(in), optional :: file_prefix ! 【新增】声明
        
        ! 输出参数 (将内存生命周期交还给上层，方便上层继续用于粒子滤波更新)
        type(uq_state_type), intent(out) :: initial_state_out
        type(uq_state_type), intent(out) :: final_state_out
        
        ! 局部变量
        class(uq_propagator_base), allocatable :: propagator
        character(len=MAX_STRING_LEN) :: actual_prefix
        real(DP) :: cpu_start, cpu_end
        
        call print_separator('不确定性/误差传播 (UQ API)')
        write(*,*) '[INFO] 正在初始化分布参数, 粒子数: ', n_particles
        
        ! 1. 根据传入的标称轨道和协方差，生成初始粒子分布
        call initial_state_out%allocate_memory(6, n_particles)
        call init_random_seed(.true.) ! 如果你希望每次运行轨迹不同，可以改为 .false.
        call generate_multivariate_normal(nominal_state, initial_cov, initial_state_out%samples)
        call initial_state_out%compute_moments()
        
        ! 2. 实例化对应的传播器
        select case(method_switch)
            case(METHOD_MC)
                allocate(uq_mc_propagator :: propagator)
            case(METHOD_DA)
                allocate(uq_da_propagator :: propagator)
            case default
                write(*,*) '[ERROR] UQ API: 未知的传播方法开关!'
                return
        end select
        
        ! =======================================================
        ! 3. 配置传播器参数 (核心改动区域)
        ! =======================================================
        propagator%epoch0 = epoch0 ! 【新增】将历元基准注入给动力学传播引擎
        ! 如果用户显式传了积分器开关，才去覆盖默认值；否则使用基类默认的 METHOD_RKF78
        if (present(integrator_switch)) then
            call propagator%set_integrator(integrator_switch)
        end if
        call propagator%set_verbosity(.true.)
        
        ! 4. 核心计算
        write(*,*) '--------------------------------------------------'
        write(*,*) ' 开始执行传播计算...'
        write(*,*) ' 使用算法: ', trim(propagator%get_method_name())
        
        call cpu_time(cpu_start)
        call propagator%propagate(t_start, t_end, initial_state_out, final_state_out)
        call cpu_time(cpu_end)
        
        write(*,*) ' ✅ 传播计算完成！耗时: ', cpu_end - cpu_start, ' 秒'
        write(*,*) '--------------------------------------------------'
        
        ! 5. 计算并打印统计信息 (此时 final_state_out 内已经由传播器自动调用过了 compute_moments)
        call display_uq_results(initial_state_out, final_state_out)
        
        ! 6. 文件保存控制
        if (save_results_to_file) then
            ! 如果用户传了前缀就用用户的，没传就给个安全默认值
            if (present(file_prefix)) then
                actual_prefix = file_prefix
            else
                actual_prefix = './output/uq_default' 
            end if
            
            ! 传入前缀进行保存
            call save_uq_results(final_state_out, actual_prefix)
        end if
        
    end subroutine run_uq_propagation
    
    ! ====================================================================
    ! 以下为内部辅助例程 
    ! ====================================================================
    subroutine display_uq_results(initial_state, final_state)
        type(uq_state_type), intent(in) :: initial_state, final_state
        real(DP) :: std_dev(6)
        integer :: i
        
        call print_separator('传播结果统计摘要')
        call print_vector(initial_state%mean, '初始均值 (Mean_0):', '(6(ES14.5, 1X))')
        
        do i = 1, 6
            std_dev(i) = sqrt(initial_state%cov(i,i))
        end do
        call print_vector(std_dev, '初始标准差 (1-Sigma_0):', '(6(ES14.5, 1X))')
        
        write(*,*) '--------------------------------------------------'
        call print_vector(final_state%mean, '最终均值 (Mean_f):', '(6(ES14.5, 1X))')
        
        do i = 1, 6
            std_dev(i) = sqrt(final_state%cov(i,i))
        end do
        call print_vector(std_dev, '最终标准差 (1-Sigma_f):', '(6(ES14.5, 1X))')
    end subroutine display_uq_results

    subroutine save_uq_results(final_state, file_prefix)
        type(uq_state_type), intent(in) :: final_state
        character(len=*), intent(in) :: file_prefix ! 【新增】用户决定的文件前缀路径
        
        integer :: file_unit, i, dim, n_particles
        character(len=MAX_STRING_LEN) :: filepath_particles, filepath_stats
        
        dim = size(final_state%samples, 1)
        n_particles = size(final_state%samples, 2)
        
        ! 智能拼接路径 (使用用户传入的前缀，避开全局变量中的隐藏换行符)
        filepath_particles = trim(file_prefix) // '_particles.csv'
        filepath_stats = trim(file_prefix) // '_stats.csv'
        
        ! 1. 存粒子
        open(newunit=file_unit, file=trim(filepath_particles), status='replace', action='write')
        write(file_unit, '(A)') 'x,y,z,vx,vy,vz'
        do i = 1, n_particles
            write(file_unit, '(*(ES22.14, :, ","))') final_state%samples(:, i)
        end do
        close(file_unit)
        
        ! 2. 存统计矩
        open(newunit=file_unit, file=trim(filepath_stats), status='replace', action='write')
        write(file_unit, '(A)') '# Mean'
        write(file_unit, '(*(ES22.14, :, ","))') final_state%mean(:)
        write(file_unit, '(A)') '# Covariance Matrix'
        do i = 1, dim
            write(file_unit, '(*(ES22.14, :, ","))') final_state%cov(i, :)
        end do
        close(file_unit)
        
        write(*,*) '✅ 结果已保存至:'
        write(*,*) '  ', trim(filepath_particles)
        write(*,*) '  ', trim(filepath_stats)
    end subroutine save_uq_results

end module pod_uq_propagation