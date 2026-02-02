module cat_orbit_propagation
    use cat_global, only: DP, MAX_STRING_LEN, output_directory
    use cat_config, only: config
    use cat_utils, only: print_separator, print_vector, save_matrix_to_file, &
                        get_user_choice, get_user_real, get_user_string, &
                        confirm_action, pause_execution
    use cat_frame_module, only: cartesian_to_keplerian, keplerian_to_cartesian
    use cat_integrator_module, only: rk4_integrate, rkf45_integrate
    use cat_force_model_module, only: compute_acceleration
    
    implicit none
    
    ! 轨道状态类型
    type orbit_state
        real(DP), dimension(6) :: state  ! [x, y, z, vx, vy, vz]
        real(DP) :: epoch  ! 历元时间（秒）
        character(len=MAX_STRING_LEN) :: frame  ! 参考坐标系
    end type orbit_state
    
    ! 传播结果类型
    type propagation_result
        integer :: n_steps
        real(DP), allocatable, dimension(:) :: times
        real(DP), allocatable, dimension(:,:) :: states
        real(DP), allocatable, dimension(:,:) :: keplerian_elements
    end type propagation_result
    
contains

    subroutine run_orbit_propagation()
        type(orbit_state) :: initial_state
        type(propagation_result) :: result
        real(DP) :: final_time, step_size
        integer :: integrator_choice
        logical :: save_to_file
        
        call print_separator('轨道传播 (OP)')
        
        ! 获取初始轨道状态
        call get_initial_orbit_state(initial_state)
        
        ! 获取传播参数
        call get_propagation_parameters(final_time, step_size, integrator_choice)
        
        ! 执行轨道传播
        call propagate_orbit(initial_state, final_time, step_size, integrator_choice, result)
        
        ! 显示结果
        call display_propagation_results(result)
        
        ! 保存结果
        save_to_file = confirm_action('是否保存传播结果到文件')
        if (save_to_file) then
            call save_propagation_results(result)
        end if
        
        ! 清理内存
        call cleanup_propagation_result(result)
        
        call pause_execution()
    end subroutine run_orbit_propagation
    
    subroutine get_initial_orbit_state(initial_state)
        type(orbit_state), intent(out) :: initial_state
        integer :: input_choice
        real(DP), dimension(6) :: cartesian_state
        real(DP), dimension(6) :: keplerian_elements
        
        write(*, *) '选择初始轨道状态输入方式:'
        write(*, *) '1. 笛卡尔坐标 (x, y, z, vx, vy, vz)'
        write(*, *) '2. 开普勒轨道根数 (a, e, i, Ω, ω, ν)'
        
        input_choice = get_user_choice('请选择 (1-2): ', 1, 2)
        
        select case (input_choice)
            case (1)
                ! 输入笛卡尔坐标
                write(*, *) '请输入初始笛卡尔状态 (km, km/s):'
                cartesian_state(1) = get_user_real('x (km): ', -100000.0_DP, 100000.0_DP)
                cartesian_state(2) = get_user_real('y (km): ', -100000.0_DP, 100000.0_DP)
                cartesian_state(3) = get_user_real('z (km): ', -100000.0_DP, 100000.0_DP)
                cartesian_state(4) = get_user_real('vx (km/s): ', -20.0_DP, 20.0_DP)
                cartesian_state(5) = get_user_real('vy (km/s): ', -20.0_DP, 20.0_DP)
                cartesian_state(6) = get_user_real('vz (km/s): ', -20.0_DP, 20.0_DP)
                
                initial_state%state = cartesian_state
                
            case (2)
                ! 输入开普勒轨道根数
                write(*, *) '请输入初始开普勒轨道根数:'
                keplerian_elements(1) = get_user_real('半长轴 a (km): ', 1000.0_DP, 100000.0_DP)
                keplerian_elements(2) = get_user_real('偏心率 e: ', 0.0_DP, 0.99_DP)
                keplerian_elements(3) = get_user_real('轨道倾角 i (度): ', 0.0_DP, 180.0_DP)
                keplerian_elements(4) = get_user_real('升交点赤经 Ω (度): ', 0.0_DP, 360.0_DP)
                keplerian_elements(5) = get_user_real('近地点幅角 ω (度): ', 0.0_DP, 360.0_DP)
                keplerian_elements(6) = get_user_real('真近点角 ν (度): ', 0.0_DP, 360.0_DP)
                
                ! 转换为笛卡尔坐标
                call keplerian_to_cartesian(keplerian_elements, cartesian_state)
                initial_state%state = cartesian_state
        end select
        
        initial_state%epoch = 0.0_DP
        initial_state%frame = trim(config%reference_frame)
        
        write(*, *) '初始轨道状态:'
        call print_vector(initial_state%state, '笛卡尔状态 (km, km/s)')
    end subroutine get_initial_orbit_state
    
    subroutine get_propagation_parameters(final_time, step_size, integrator_choice)
        real(DP), intent(out) :: final_time, step_size
        integer, intent(out) :: integrator_choice
        
        write(*, *) '传播参数设置:'
        final_time = get_user_real('传播时间 (秒): ', 1.0_DP, 86400.0_DP)
        step_size = get_user_real('传播步长 (秒): ', 1.0_DP, final_time)
        
        write(*, *) '选择积分器:'
        write(*, *) '1. RK4 (四阶龙格库塔)'
        write(*, *) '2. RKF45 (自适应步长)'
        
        integrator_choice = get_user_choice('请选择积分器 (1-2): ', 1, 2)
    end subroutine get_propagation_parameters
    
    subroutine propagate_orbit(initial_state, final_time, step_size, integrator_choice, result)
        type(orbit_state), intent(in) :: initial_state
        real(DP), intent(in) :: final_time, step_size
        integer, intent(in) :: integrator_choice
        type(propagation_result), intent(out) :: result
        
        integer :: n_steps, i
        real(DP), dimension(6) :: current_state
        real(DP), dimension(6) :: keplerian_elements
        
        ! 计算步数
        n_steps = int(final_time / step_size) + 1
        
        ! 分配内存
        allocate(result%times(n_steps))
        allocate(result%states(n_steps, 6))
        allocate(result%keplerian_elements(n_steps, 6))
        result%n_steps = n_steps
        
        ! 初始化
        result%times(1) = 0.0_DP
        result%states(1, :) = initial_state%state
        
        ! 计算初始开普勒轨道根数
        call cartesian_to_keplerian(initial_state%state, keplerian_elements)
        result%keplerian_elements(1, :) = keplerian_elements
        
        write(*, *) '开始轨道传播...'
        write(*, *) '总步数: ', n_steps
        
        ! 执行传播
        current_state = initial_state%state
        do i = 2, n_steps
            result%times(i) = (i - 1) * step_size
            
            ! 选择积分器
            select case (integrator_choice)
                case (1)
                    call rk4_integrate(current_state, step_size, result%times(i), current_state)
                case (2)
                    call rkf45_integrate(current_state, step_size, result%times(i), current_state)
            end select
            
            result%states(i, :) = current_state
            
            ! 计算开普勒轨道根数
            call cartesian_to_keplerian(current_state, keplerian_elements)
            result%keplerian_elements(i, :) = keplerian_elements
            
            ! 显示进度
            if (mod(i, max(1, n_steps/10)) == 0) then
                write(*, '(A,I3,A)') '传播进度: ', int(100.0_DP * i / n_steps), '%'
            end if
        end do
        
        write(*, *) '轨道传播完成!'
    end subroutine propagate_orbit
    
    subroutine display_propagation_results(result)
        type(propagation_result), intent(in) :: result
        
        write(*, *) '传播结果摘要:'
        write(*, *) '总步数: ', result%n_steps
        write(*, *) '初始时间: ', result%times(1), ' 秒'
        write(*, *) '结束时间: ', result%times(result%n_steps), ' 秒'
        write(*, *) '传播时长: ', result%times(result%n_steps) - result%times(1), ' 秒'
        
        write(*, *) '初始状态:'
        call print_vector(result%states(1, :), '笛卡尔状态 (km, km/s)')
        
        write(*, *) '最终状态:'
        call print_vector(result%states(result%n_steps, :), '笛卡尔状态 (km, km/s)')
        
        write(*, *) '初始开普勒轨道根数:'
        call print_vector(result%keplerian_elements(1, :), '开普勒根数 (km, -, 度)')
        
        write(*, *) '最终开普勒轨道根数:'
        call print_vector(result%keplerian_elements(result%n_steps, :), '开普勒根数 (km, -, 度)')
    end subroutine display_propagation_results
    
    subroutine save_propagation_results(result)
        type(propagation_result), intent(in) :: result
        character(len=MAX_STRING_LEN) :: filename
        integer :: unit, i
        
        ! 保存笛卡尔状态
        filename = trim(output_directory) // 'orbit_propagation_cartesian.csv'
        open(newunit=unit, file=filename, status='replace', action='write')
        
        write(unit, '(A)') '# 时间(s), x(km), y(km), z(km), vx(km/s), vy(km/s), vz(km/s)'
        do i = 1, result%n_steps
            write(unit, '(F12.3,6(",",F15.6))') result%times(i), result%states(i, :)
        end do
        
        close(unit)
        
        ! 保存开普勒轨道根数
        filename = trim(output_directory) // 'orbit_propagation_keplerian.csv'
        open(newunit=unit, file=filename, status='replace', action='write')
        
        write(unit, '(A)') '# 时间(s), a(km), e, i(deg), Omega(deg), omega(deg), nu(deg)'
        do i = 1, result%n_steps
            write(unit, '(F12.3,6(",",F15.6))') result%times(i), result%keplerian_elements(i, :)
        end do
        
        close(unit)
        
        write(*, *) '传播结果已保存到:'
        write(*, *) '  笛卡尔状态: ', trim(output_directory) // 'orbit_propagation_cartesian.csv'
        write(*, *) '  开普勒根数: ', trim(output_directory) // 'orbit_propagation_keplerian.csv'
    end subroutine save_propagation_results
    
    subroutine cleanup_propagation_result(result)
        type(propagation_result), intent(inout) :: result
        
        if (allocated(result%times)) deallocate(result%times)
        if (allocated(result%states)) deallocate(result%states)
        if (allocated(result%keplerian_elements)) deallocate(result%keplerian_elements)
    end subroutine cleanup_propagation_result

end module cat_orbit_propagation
