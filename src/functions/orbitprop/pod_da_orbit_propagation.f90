module pod_da_orbit_propagation
    use pod_global, only: DP, MAX_STRING_LEN, output_directory
    use pod_config, only: config
    use pod_utils, only: print_separator, print_vector, &
                        get_user_choice, get_user_real, get_user_string, &
                        confirm_action, pause_execution
    use pod_spice, only: str2et
    
    ! 引入 DA 积分器和 DACE 核心类
    use pod_da_integrator_module, only: da_adaptive_step_integrate, METHOD_RKF45, METHOD_RKF78
    use pod_dace_classes
    use pod_da_force_model_module, only: set_propagation_epoch
    
    implicit none
    
    ! =========================================================
    ! 轨道状态类型 (DA版 - 面向用户的物理参数)
    ! =========================================================
    type da_orbit_state
        real(DP), dimension(6) :: nominal_state  ! 真实的物理标称状态 [x, y, z, vx, vy, vz] (km, km/s)
        real(DP) :: epoch                        ! 真实的物理历元基准 (TDB 秒)
        integer :: da_order                      ! 用户选择的 DA 截断阶数
    end type da_orbit_state
    
    ! =========================================================
    ! 传播结果类型 (DA版 - 承载高维代数多项式)
    ! =========================================================
    type da_propagation_result
        integer :: n_steps
        real(DP), allocatable, dimension(:) :: times              ! 传播经历的时间 (秒)
        type(AlgebraicVector), allocatable, dimension(:) :: states ! 轨道历经的 DA 状态 (自带物理量纲)
    end type da_propagation_result
    
contains

    !> 交互主程序：DA 高精度轨道与不确定性传播
    subroutine run_da_orbit_propagation()
        type(da_orbit_state) :: initial_state
        type(da_propagation_result) :: result
        real(DP) :: final_time
        integer :: integrator_choice
        logical :: save_to_file
        
        call print_separator('高精度微分代数 (DA) 轨道传播')
        
        ! 1. 获取初始状态与 DA 阶数
        call get_initial_da_orbit_state(initial_state)
        
        ! 2. 获取传播参数
        call get_da_propagation_parameters(final_time, integrator_choice)
        
        ! 3. 核心计算：执行 DA 轨道传播 (自动处理量纲、历元、DA注入)
        call propagate_da_orbit(initial_state, final_time, integrator_choice, result)
        
        ! 4. 显示结果 (标称状态 + STM 提取)
        call display_da_propagation_results(result)

        ! 额外展示：打印最终 DA 状态多项式（包含不确定性信息）
        call print_da_vector(result%states(result%n_steps), "积分后的最终 DA 状态多项式")
        
        ! 5. 询问并保存结果
        save_to_file = confirm_action('是否将标称轨迹结果保存到 CSV 文件')
        if (save_to_file) then
            call save_da_propagation_results(result)
        end if
        
        ! 6. 清理内存 (释放 DA 句柄和数组)
        call cleanup_da_propagation_result(result)
        
        call pause_execution()
    end subroutine run_da_orbit_propagation
    
    
    subroutine get_initial_da_orbit_state(initial_state)
        type(da_orbit_state), intent(out) :: initial_state
        character(len=MAX_STRING_LEN) :: utc_string
        real(DP), dimension(6) :: cartesian_state
        
        write(*, *) '请输入初始标称笛卡尔状态 (真实物理量 km, km/s):'
        cartesian_state(1) = get_user_real('  x (km): ', -1.0e8_DP, 1.0e8_DP)
        cartesian_state(2) = get_user_real('  y (km): ', -1.0e8_DP, 1.0e8_DP)
        cartesian_state(3) = get_user_real('  z (km): ', -1.0e8_DP, 1.0e8_DP)
        cartesian_state(4) = get_user_real(' vx (km/s): ', -50.0_DP, 50.0_DP)
        cartesian_state(5) = get_user_real(' vy (km/s): ', -50.0_DP, 50.0_DP)
        cartesian_state(6) = get_user_real(' vz (km/s): ', -50.0_DP, 50.0_DP)
        
        initial_state%nominal_state = cartesian_state
        
        write(*, *) '请输入历元时间 (UTC):'
        call get_user_string('  UTC (例如 2024-03-09T12:00:00): ', utc_string)
        call str2et(trim(utc_string), initial_state%epoch)
        
        write(*, *) '请配置微分代数 (DA) 引擎参数:'
        ! 允许用户选择 DA 截断阶数 (通常 1 阶用于 STM，2~5阶用于高阶不确定性传播)
        initial_state%da_order = get_user_choice('  选择 DA 截断阶数 (1-15): ', 1, 15)
        
        write(*, *) '初始轨道状态与 DA 配置已记录。'
    end subroutine get_initial_da_orbit_state
    
    
    subroutine get_da_propagation_parameters(final_time, integrator_choice)
        real(DP), intent(out) :: final_time
        integer, intent(out) :: integrator_choice
        
        write(*, *) '传播参数设置:'
        final_time = get_user_real('  总传播时长 (秒): ', 1.0_DP, 2592000.0_DP)
        
        write(*, *) '选择核心积分器:'
        write(*, *) '  1. RKF 4(5) [推荐中低精度，较密集的输出点]'
        write(*, *) '  2. RKF 7(8) [推荐高精度，极大的步长跨度]'
        
        integrator_choice = get_user_choice('请选择积分器 (1-2): ', 1, 2)
    end subroutine get_da_propagation_parameters
    
    
    !> 核心 DA 封装引擎
    subroutine propagate_da_orbit(initial_state, final_time, integrator_choice, result)
        type(da_orbit_state), intent(in) :: initial_state
        real(DP), intent(in) :: final_time
        integer, intent(in) :: integrator_choice
        type(da_propagation_result), intent(out) :: result
        
        type(AlgebraicVector) :: nondim_state
        real(DP) :: t_start_nondim, t_end_nondim
        integer :: actual_method, i, j
        
        write(*, *) '-----------------------------------------'
        write(*, *) '正在准备底层物理环境与 DA 数值积分器...'
        
        ! 1. 初始化 DACE 引擎 (阶数为用户输入，自变量固定为 6 个轨道根数)
        call dace_initialize(initial_state%da_order, 6)
        
        ! 2. 历元基准注入
        call set_propagation_epoch(initial_state%epoch)
        
        ! 3. 初始状态: DA 变量注入与无量纲化
        call nondim_state%init(6)
        do i = 1, 3
            ! 物理值 + da_var 注入独立微小偏差，随后除以尺度归一化
            nondim_state%elements(i)   = (initial_state%nominal_state(i) + da_var(i)) / config%LU

            nondim_state%elements(i+3) = (initial_state%nominal_state(i+3) + da_var(i+3)) / config%VU
        end do
        
        t_start_nondim = 0.0_DP
        t_end_nondim   = final_time / config%TU
        
        ! 4. 映射积分器选项
        if (integrator_choice == 1) then
            actual_method = METHOD_RKF45
        else
            actual_method = METHOD_RKF78
        end if
        
        write(*, *) 'DA 积分器已启动，正在自适应计算中 (过程可能较长)...'
        
        ! 5. 调用 DA 积分器
        call da_adaptive_step_integrate( &
            state = nondim_state, &
            t_start = t_start_nondim, &
            t_end = t_end_nondim, &
            integrator_method = actual_method, &
            times = result%times, &
            states = result%states, &
            n_steps = result%n_steps &
        )
        
        ! 6. 结果原地还原物理量纲
        result%times = result%times * config%TU
        do i = 1, result%n_steps
            do j = 1, 3
                ! 调用底层的 DA * 实数 重载运算
                result%states(i)%elements(j)   = result%states(i)%elements(j) * config%LU
                result%states(i)%elements(j+3) = result%states(i)%elements(j+3) * config%VU
            end do
        end do
        
        write(*, *) 'DA 轨道传播计算圆满完成!'
        write(*, *) '-----------------------------------------'
    end subroutine propagate_da_orbit
    
    
    subroutine display_da_propagation_results(result)
        type(da_propagation_result), intent(in) :: result
        integer :: i
        
        write(*, *) ''
        write(*, *) 'DA 传播结果摘要:'
        write(*, *) '  自适应总步数: ', result%n_steps
        write(*, *) '  传播物理时长: ', result%times(result%n_steps), ' 秒'
        
        write(*, *) '  最终标称状态 (常数项, km, km/s):'
        ! cons() 会提取出 0 阶标称轨迹
        call print_vector(result%states(result%n_steps)%cons(), '    ')
        
        write(*, *) '  -----------------------------------------'
        write(*, *) '  ✨ 最终状态的 STM (状态转移矩阵) 左上角 3x3 ✨'
        write(*, *) '  -----------------------------------------'
        do i = 1, 3
            write(*, "(4X, 3(ES15.6, 2X))") &
                result%states(result%n_steps)%elements(i)%get_deriv_value(1), &
                result%states(result%n_steps)%elements(i)%get_deriv_value(2), &
                result%states(result%n_steps)%elements(i)%get_deriv_value(3)
        end do
    end subroutine display_da_propagation_results
    
    
    subroutine save_da_propagation_results(result)
        type(da_propagation_result), intent(in) :: result
        character(len=MAX_STRING_LEN) :: filename
        integer :: unit, i
        real(DP), dimension(6) :: nominal_vec
        
        filename = './output/da_orbit_propagation_nominal.csv'
        open(newunit=unit, file=filename, status='replace', action='write')
        
        write(unit, '(A)') '# Time(s), X(km), Y(km), Z(km), VX(km/s), VY(km/s), VZ(km/s)'
        do i = 1, result%n_steps
            ! 提取当前时刻的标称实数轨迹
            nominal_vec = result%states(i)%cons()
            write(unit, '(F16.6, 6(",",F20.12))') result%times(i), nominal_vec
        end do
        
        close(unit)
        
        write(*, *) '✅ DA 标称传播结果已成功保存到:'
        write(*, *) '   -> ', trim(filename)
    end subroutine save_da_propagation_results

    !> 傻瓜式 DA 向量打印工具
    !> 允许用户传入一个标题，自动美化输出 AlgebraicVector 的多项式内容
    subroutine print_da_vector(da_vec, title)
        use pod_dace_classes, only: AlgebraicVector
        implicit none
        
        class(AlgebraicVector), intent(in) :: da_vec
        character(len=*), intent(in), optional :: title
        
        write(*, *) ''
        ! 如果用户传了标题，就打个漂亮的表头
        if (present(title)) then
            write(*, *) '=================================================='
            write(*, *) ' 📊 ', trim(title)
            write(*, *) '=================================================='
        end if
        
        ! 安全检查：防止空向量崩溃
        if (da_vec%size <= 0 .or. .not. allocated(da_vec%elements)) then
            write(*, *) '  [警告] 该 DA 向量为空或未初始化！'
            write(*, *) '=================================================='
            return
        end if
        
        ! 核心魔法：直接触发 AlgebraicVector 内置的面向对象 print 方法
        call da_vec%print()
        
        if (present(title)) then
            write(*, *) '=================================================='
        end if
        
    end subroutine print_da_vector
    
    
    subroutine cleanup_da_propagation_result(result)
        type(da_propagation_result), intent(inout) :: result
        integer :: i
        
        if (allocated(result%times)) deallocate(result%times)
        
        if (allocated(result%states)) then
            ! 显式销毁 DA 句柄防内存泄漏
            do i = 1, size(result%states)
                call result%states(i)%destroy()
            end do
            deallocate(result%states)
        end if
    end subroutine cleanup_da_propagation_result

end module pod_da_orbit_propagation