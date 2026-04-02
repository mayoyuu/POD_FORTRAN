module pod_orbit_propagation
    use pod_global, only: DP, MAX_STRING_LEN, output_directory
    use pod_config, only: config
    use pod_utils, only: print_separator, print_vector, &
                        get_user_choice, get_user_real, get_user_string, &
                        confirm_action, pause_execution
    use pod_spice, only: str2et
    use pod_integrator_module, only: adaptive_step_integrate, METHOD_RKF45, METHOD_RKF78
    use pod_force_model_module, only: set_propagation_epoch,current_epoch0
    
    implicit none
    
    ! =========================================================
    ! 轨道状态类型 (面向用户的真实物理量)
    ! =========================================================
    type orbit_state
        real(DP), dimension(6) :: state  ! 真实物理状态 [x, y, z, vx, vy, vz] (km, km/s)
        real(DP) :: epoch                ! 真实的物理历元基准 (TDB 秒)
    end type orbit_state
    
    ! =========================================================
    ! 传播结果类型 (面向用户的真实物理量)
    ! =========================================================
    type propagation_result
        integer :: n_steps
        real(DP), allocatable, dimension(:) :: times    ! 传播经历的时间 (秒)
        real(DP), allocatable, dimension(:,:) :: states ! 轨道历经的真实物理状态 (km, km/s)
    end type propagation_result
    
contains

    !>交互主程序：只需调用此函数，跟着终端提示走即可
    subroutine run_orbit_propagation()
        type(orbit_state) :: initial_state
        type(propagation_result) :: result
        real(DP) :: final_time
        integer :: integrator_choice
        logical :: save_to_file
        
        call print_separator('高精度轨道传播 (OP)')
        
        ! 1. 获取初始状态
        call get_initial_orbit_state(initial_state)
        
        ! 2. 获取传播参数 (已删去手动设定步长)
        call get_propagation_parameters(final_time, integrator_choice)
        
        ! 3. 核心计算：执行轨道传播 (内部自动处理量纲与历元)
        call propagate_orbit(initial_state, final_time, integrator_choice, result)
        
        ! 4. 显示结果
        call display_propagation_results(result)
        
        ! 5. 询问并保存结果
        save_to_file = confirm_action('是否将笛卡尔传播结果保存到 CSV 文件')
        if (save_to_file) then
            call save_propagation_results(result)
        end if
        
        ! 6. 清理内存
        call cleanup_propagation_result(result)
        
        call pause_execution()
    end subroutine run_orbit_propagation
    
    
    subroutine get_initial_orbit_state(initial_state)
        type(orbit_state), intent(out) :: initial_state
        character(len=MAX_STRING_LEN) :: utc_string
        real(DP), dimension(6) :: cartesian_state
        
        write(*, *) '请输入初始笛卡尔状态 (真实物理量 km, km/s):'
        cartesian_state(1) = get_user_real('  x (km): ', -1.0e8_DP, 1.0e8_DP)
        cartesian_state(2) = get_user_real('  y (km): ', -1.0e8_DP, 1.0e8_DP)
        cartesian_state(3) = get_user_real('  z (km): ', -1.0e8_DP, 1.0e8_DP)
        cartesian_state(4) = get_user_real(' vx (km/s): ', -50.0_DP, 50.0_DP)
        cartesian_state(5) = get_user_real(' vy (km/s): ', -50.0_DP, 50.0_DP)
        cartesian_state(6) = get_user_real(' vz (km/s): ', -50.0_DP, 50.0_DP)
        initial_state%state = cartesian_state
        
        write(*, *) '请输入历元时间 (UTC):'
        call get_user_string('  UTC (例如 2024-03-09T12:00:00): ', utc_string)
        call str2et(trim(utc_string), initial_state%epoch)
        
        write(*, *) '初始轨道状态已记录。'
    end subroutine get_initial_orbit_state
    
    
    subroutine get_propagation_parameters(final_time, integrator_choice)
        real(DP), intent(out) :: final_time
        integer, intent(out) :: integrator_choice
        
        write(*, *) '传播参数设置:'
        ! 默认最大支持传播 30 天
        final_time = get_user_real('  总传播时长 (秒): ', 1.0_DP, 2592000.0_DP)
        
        write(*, *) '选择核心积分器:'
        write(*, *) '  1. RKF 4(5) [推荐中低精度，较密集的输出点]'
        write(*, *) '  2. RKF 7(8) [推荐高精度，极大的步长跨度]'
        
        integrator_choice = get_user_choice('请选择积分器 (1-2): ', 1, 2)
    end subroutine get_propagation_parameters
    
    
    !> 核心封装引擎：对外部隐藏所有数学/物理处理的脏活累活
    subroutine propagate_orbit(initial_state, final_time, integrator_choice, result)
        type(orbit_state), intent(in) :: initial_state
        real(DP), intent(in) :: final_time
        integer, intent(in) :: integrator_choice
        type(propagation_result), intent(out) :: result
        
        real(DP), dimension(6) :: nondim_state
        real(DP) :: t_start_nondim, t_end_nondim
        integer :: actual_method
        
        write(*, *) '-----------------------------------------'
        write(*, *) '正在准备底层物理环境与数值积分器...'
        
        ! 1. 历元基准 (Context Injection)
        call set_propagation_epoch(initial_state%epoch)
        
        ! 2. 初始状态去量纲化 (Nondimensionalize)
        nondim_state(1:3) = initial_state%state(1:3) / config%LU
        nondim_state(4:6) = initial_state%state(4:6) / config%VU
        t_start_nondim = 0.0_DP
        t_end_nondim   = final_time / config%TU
        
        ! 3. 映射积分器选项
        if (integrator_choice == 1) then
            actual_method = METHOD_RKF45
        else
            actual_method = METHOD_RKF78
        end if
        
        write(*, *) '积分器已启动，正在自适应计算中...'
        
        ! 4. 调用纯数学积分器 (自动分配 result%times 和 result%states 的内存)
        call adaptive_step_integrate( &
            state = nondim_state, &
            t_start = t_start_nondim, &
            t_end = t_end_nondim, &
            integrator_method = actual_method, &
            times = result%times, &
            states = result%states, &
            n_steps = result%n_steps &
        )
        
        ! 5. 结果原地还原物理量纲 (Dimensionalize In-Place)
        ! 积分器输出的是无量纲数组，需要乘回物理单位
        result%times = result%times * config%TU
        result%states(:, 1:3) = result%states(:, 1:3) * config%LU
        result%states(:, 4:6) = result%states(:, 4:6) * config%VU
        
        write(*, *) '轨道传播计算圆满完成!'
        write(*, *) '-----------------------------------------'
    end subroutine propagate_orbit
    
    
    subroutine display_propagation_results(result)
        type(propagation_result), intent(in) :: result
        
        write(*, *) ''
        write(*, *) '传播结果摘要:'
        write(*, *) '  自适应总步数: ', result%n_steps
        write(*, *) '  传播物理时长: ', result%times(result%n_steps), ' 秒'
        
        write(*, *) '  初始状态 (km, km/s):'
        call print_vector(result%states(1, :), '    ')
        
        write(*, *) '  最终状态 (km, km/s):'
        call print_vector(result%states(result%n_steps, :), '    ')
    end subroutine display_propagation_results
    
    
    subroutine save_propagation_results(result)
        type(propagation_result), intent(in) :: result
        character(len=MAX_STRING_LEN) :: filename
        integer :: unit, i
        
        ! 保存笛卡尔状态
        ! filename = trim(output_directory) // 'orbit_propagation_cartesian.csv'
        filename = './output/orbit_propagation_cartesian.csv'
        open(newunit=unit, file=filename, status='replace', action='write')
        
        write(unit, '(A)') '# Time(s), X(km), Y(km), Z(km), VX(km/s), VY(km/s), VZ(km/s)'
        do i = 1, result%n_steps
            write(unit, '(F16.6, 6(",",F20.12))') result%times(i), result%states(i, :)
        end do
        
        close(unit)
        
        write(*, *) '✅ 传播结果已成功保存到:'
        write(*, *) '   -> ', trim(filename)
    end subroutine save_propagation_results
    
    
    subroutine cleanup_propagation_result(result)
        type(propagation_result), intent(inout) :: result
        
        if (allocated(result%times)) deallocate(result%times)
        if (allocated(result%states)) deallocate(result%states)
    end subroutine cleanup_propagation_result

end module pod_orbit_propagation