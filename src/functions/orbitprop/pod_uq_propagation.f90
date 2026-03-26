module pod_uq_propagation
    use pod_global, only: DP, MAX_STRING_LEN, output_directory
    use pod_utils,  only: print_separator, get_user_choice, get_user_real, confirm_action
    
    ! 引入底层计算库 (核心黑盒)
    use pod_uq_base_module,  only: uq_propagator_base
    use pod_uq_da_module,    only: uq_da_propagator
    use pod_uq_mc_module,    only: uq_mc_propagator
    use pod_uq_state_module, only: uq_state_type

    use pod_random_module, only: generate_multivariate_normal
    
    implicit none
    
contains
    ! 位于 pod_uq_propagation.f90 或专门的 sampling 模块中
    subroutine get_initial_uq_state(initial_state)
        type(uq_state_type), intent(inout) :: initial_state
        integer :: input_choice
        character(len=MAX_STRING_LEN) :: filepath
        real(DP), allocatable :: user_mean(:), user_cov(:,:)
        integer :: n_particles, dim
        
        dim = 6 ! 默认轨道状态维度
        
        write(*,*) '请选择初始粒子分布的录入方式:'
        write(*,*) '1. 从文件中读取大规模粒子点云'
        write(*,*) '2. 手动输入/读取均值和协方差，并随机生成粒子'
        input_choice = get_user_choice('请选择 (1-2): ', 1, 2)
        
        select case (input_choice)
            case (1)
                ! 用户输入文件路径
                call get_user_string('请输入粒子文件路径: ', filepath)
                
                ! 解析文件获取粒子数，然后分配内存
                ! call count_lines_in_file(filepath, n_particles)
                ! call initial_state%allocate_memory(dim, n_particles)
                
                ! 将文件数据填入 initial_state%samples
                ! call read_particles_from_csv(filepath, initial_state%samples)
                
            case (2)
                n_particles = int(get_user_real('请输入要生成的粒子数量: ', 10.0_DP, 1000000.0_DP))
                
                ! 1. 在这里获取 user_mean 和 user_cov (可以通过终端输入或读取一个小配置文件)
                ! ...
                
                ! 2. 分配状态容器内存
                call initial_state%allocate_memory(dim, n_particles)
                
                ! 3. 调用你的随机数生成引擎 (例如 Cholesky分解 + Box-Muller)
                call  generate_multivariate_normal(user_mean, user_cov, initial_state%samples)
                
        end select
        
        ! 生成或读取完样本后，立刻计算一次初始矩，以备后续比对
        call initial_state%compute_moments()
        
    end subroutine get_initial_uq_state

    subroutine run_uq_propagation()
        type(uq_state_type) :: initial_state, final_state
        class(uq_propagator_base), allocatable :: propagator
        real(DP) :: t_start, t_end, dt
        integer :: method_choice
        logical :: save_to_file
        
        call print_separator('不确定性/误差传播 (UP/UQ)')
        
        ! 1. [I/O] 获取初始分布 (比如输入均值和协方差，并在内部生成粒子)
        call get_initial_uq_state(initial_state)
        
        ! 2. [I/O] 获取传播参数和选择方法 (DA, MC, UT等)
        call get_uq_parameters(t_start, t_end, dt, method_choice)
        
        ! 3. [业务逻辑] 多态实例化具体的传播器
        select case(method_choice)
            case(1)
                allocate(uq_mc_propagator :: propagator)
            case(2)
                allocate(uq_da_propagator :: propagator)
                ! 如果是DA，还可以进一步获取阶数等参数
        end select
        
        ! 4. [核心计算] 调用纯粹的底层数学库，不涉及任何屏幕输出或文件读写
        write(*,*) '开始执行不确定性传播计算...'
        call propagator%propagate(t_start, t_end, initial_state, final_state)
        write(*,*) '传播计算完成！'
        
        ! 5. [I/O] 计算统计矩并显示结果
        call final_state%compute_moments()
        call display_uq_results(initial_state, final_state)
        
        ! 6. [I/O] 询问并保存结果到文件
        save_to_file = confirm_action('是否保存粒子分布和协方差矩阵到文件')
        if (save_to_file) then
            call save_uq_results(final_state)
        end if
        
        ! 7. 清理内存
        ! ...
    end subroutine run_uq_propagation
    
    ! ================= 以下为拆分出来的具体 I/O 函数 =================

    subroutine get_initial_uq_state(initial_state)
        type(uq_state_type), intent(out) :: initial_state
        ! 在这里提示用户输入均值向量和协方差矩阵
        ! 然后调用相关的数学工具库生成正态分布的初始粒子云
    end subroutine get_initial_uq_state
    
    subroutine display_uq_results(initial_state, final_state)
        type(uq_state_type), intent(in) :: initial_state, final_state
        ! 在屏幕上打印传播前后的均值和协方差矩阵的对比
    end subroutine display_uq_results
    
    subroutine save_uq_results(final_state)
        type(uq_state_type), intent(in) :: final_state
        ! 将 final_state%samples 写成 particles.csv
        ! 将 final_state%mean 和 final_state%cov 写成 statistics.csv
    end subroutine save_uq_results

end module pod_uq_propagation