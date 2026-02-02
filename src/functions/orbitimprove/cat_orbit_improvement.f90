module cat_orbit_improvement
    use cat_global, only: DP, MAX_STRING_LEN, output_directory
    use cat_config, only: config
    use cat_utils, only: print_separator, print_vector, print_matrix, &
                        get_user_choice, get_user_real, get_user_string, &
                        confirm_action, pause_execution, save_matrix_to_file
    use cat_frame_module, only: cartesian_to_keplerian, keplerian_to_cartesian
    use cat_statistics_module, only: least_squares_fit
    use cat_measurement_model_module, only: compute_measurement, compute_measurement_jacobian
    
    implicit none
    
    ! 观测数据类型
    type observation_data
        integer :: n_obs
        real(DP), allocatable, dimension(:) :: times
        real(DP), allocatable, dimension(:,:) :: measurements  ! [range, az, el] 或 [x, y, z]
        real(DP), allocatable, dimension(:,:) :: measurement_weights
        character(len=MAX_STRING_LEN) :: measurement_type  ! 'RADAR', 'OPTICAL', 'GPS'
    end type observation_data
    
    ! 轨道改进结果类型
    type orbit_improvement_result
        real(DP), dimension(6) :: initial_state
        real(DP), dimension(6) :: improved_state
        real(DP), dimension(6,6) :: covariance_matrix
        real(DP) :: convergence_criteria
        integer :: iterations
        logical :: converged
        real(DP) :: final_residual
    end type orbit_improvement_result
    
contains

    subroutine run_orbit_improvement()
        type(observation_data) :: obs_data
        type(orbit_improvement_result) :: result
        real(DP), dimension(6) :: initial_guess
        integer :: method_choice
        logical :: save_to_file
        
        call print_separator('轨道改进 (OM)')
        
        ! 获取观测数据
        call get_observation_data(obs_data)
        
        ! 获取初始轨道猜测
        call get_initial_orbit_guess(initial_guess)
        
        ! 选择改进方法
        call select_improvement_method(method_choice)
        
        ! 执行轨道改进
        call improve_orbit(obs_data, initial_guess, method_choice, result)
        
        ! 显示结果
        call display_improvement_results(result)
        
        ! 保存结果
        save_to_file = confirm_action('是否保存改进结果到文件')
        if (save_to_file) then
            call save_improvement_results(result)
        end if
        
        ! 清理内存
        call cleanup_observation_data(obs_data)
        
        call pause_execution()
    end subroutine run_orbit_improvement
    
    subroutine get_observation_data(obs_data)
        type(observation_data), intent(out) :: obs_data
        integer :: input_choice, n_obs, i
        character(len=MAX_STRING_LEN) :: filename
        
        write(*, *) '选择观测数据输入方式:'
        write(*, *) '1. 手动输入'
        write(*, *) '2. 从文件读取'
        
        input_choice = get_user_choice('请选择 (1-2): ', 1, 2)
        
        select case (input_choice)
            case (1)
                ! 手动输入观测数据
                call get_manual_observations(obs_data)
                
            case (2)
                ! 从文件读取观测数据
                call get_user_string('请输入观测数据文件名: ', filename)
                call load_observations_from_file(filename, obs_data)
        end select
    end subroutine get_observation_data
    
    subroutine get_manual_observations(obs_data)
        type(observation_data), intent(out) :: obs_data
        integer :: n_obs, obs_type, i
        real(DP) :: time, range, az, el, x, y, z
        
        write(*, *) '选择观测类型:'
        write(*, *) '1. 雷达观测 (距离、方位角、俯仰角)'
        write(*, *) '2. 光学观测 (赤经、赤纬)'
        write(*, *) '3. GPS观测 (x, y, z)'
        
        obs_type = get_user_choice('请选择观测类型 (1-3): ', 1, 3)
        
        n_obs = get_user_choice('请输入观测数量: ', 1, 1000)
        
        ! 分配内存
        allocate(obs_data%times(n_obs))
        allocate(obs_data%measurements(n_obs, 3))
        allocate(obs_data%measurement_weights(n_obs, 3))
        obs_data%n_obs = n_obs
        
        select case (obs_type)
            case (1)
                obs_data%measurement_type = 'RADAR'
                write(*, *) '请输入雷达观测数据 (时间(s), 距离(km), 方位角(度), 俯仰角(度)):'
            case (2)
                obs_data%measurement_type = 'OPTICAL'
                write(*, *) '请输入光学观测数据 (时间(s), 赤经(度), 赤纬(度)):'
            case (3)
                obs_data%measurement_type = 'GPS'
                write(*, *) '请输入GPS观测数据 (时间(s), x(km), y(km), z(km)):'
        end select
        
        do i = 1, n_obs
            write(*, '(A,I3,A)') '观测 ', i, ':'
            obs_data%times(i) = get_user_real('时间 (秒): ', 0.0_DP, 86400.0_DP)
            
            select case (obs_type)
                case (1)
                    range = get_user_real('距离 (km): ', 0.0_DP, 100000.0_DP)
                    az = get_user_real('方位角 (度): ', 0.0_DP, 360.0_DP)
                    el = get_user_real('俯仰角 (度): ', -90.0_DP, 90.0_DP)
                    obs_data%measurements(i, :) = [range, az, el]
                case (2)
                    x = get_user_real('赤经 (度): ', 0.0_DP, 360.0_DP)
                    y = get_user_real('赤纬 (度): ', -90.0_DP, 90.0_DP)
                    obs_data%measurements(i, :) = [x, y, 0.0_DP]
                case (3)
                    x = get_user_real('x (km): ', -100000.0_DP, 100000.0_DP)
                    y = get_user_real('y (km): ', -100000.0_DP, 100000.0_DP)
                    z = get_user_real('z (km): ', -100000.0_DP, 100000.0_DP)
                    obs_data%measurements(i, :) = [x, y, z]
            end select
            
            ! 设置测量权重（简化处理）
            obs_data%measurement_weights(i, :) = [1.0_DP, 1.0_DP, 1.0_DP]
        end do
    end subroutine get_manual_observations
    
    subroutine load_observations_from_file(filename, obs_data)
        character(len=*), intent(in) :: filename
        type(observation_data), intent(out) :: obs_data
        integer :: unit, ios, i, n_obs
        character(len=MAX_STRING_LEN) :: line
        
        open(newunit=unit, file=filename, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            write(*, *) '无法打开文件: ', trim(filename)
            return
        end if
        
        ! 计算观测数量
        n_obs = 0
        do
            read(unit, '(A)', iostat=ios) line
            if (ios /= 0) exit
            if (len_trim(line) > 0 .and. line(1:1) /= '#') n_obs = n_obs + 1
        end do
        
        rewind(unit)
        
        ! 分配内存
        allocate(obs_data%times(n_obs))
        allocate(obs_data%measurements(n_obs, 3))
        allocate(obs_data%measurement_weights(n_obs, 3))
        obs_data%n_obs = n_obs
        
        ! 读取数据
        i = 1
        do
            read(unit, '(A)', iostat=ios) line
            if (ios /= 0) exit
            if (len_trim(line) > 0 .and. line(1:1) /= '#') then
                read(line, *, iostat=ios) obs_data%times(i), obs_data%measurements(i, :)
                if (ios == 0) i = i + 1
            end if
        end do
        
        close(unit)
        
        ! 设置默认权重
        obs_data%measurement_weights = 1.0_DP
        obs_data%measurement_type = 'RADAR'  ! 默认类型
        
        write(*, *) '从文件加载了 ', n_obs, ' 个观测数据'
    end subroutine load_observations_from_file
    
    subroutine get_initial_orbit_guess(initial_guess)
        real(DP), dimension(6), intent(out) :: initial_guess
        integer :: input_choice
        real(DP), dimension(6) :: keplerian_elements
        
        write(*, *) '选择初始轨道猜测输入方式:'
        write(*, *) '1. 笛卡尔坐标 (x, y, z, vx, vy, vz)'
        write(*, *) '2. 开普勒轨道根数 (a, e, i, Ω, ω, ν)'
        
        input_choice = get_user_choice('请选择 (1-2): ', 1, 2)
        
        select case (input_choice)
            case (1)
                write(*, *) '请输入初始笛卡尔状态猜测 (km, km/s):'
                initial_guess(1) = get_user_real('x (km): ', -100000.0_DP, 100000.0_DP)
                initial_guess(2) = get_user_real('y (km): ', -100000.0_DP, 100000.0_DP)
                initial_guess(3) = get_user_real('z (km): ', -100000.0_DP, 100000.0_DP)
                initial_guess(4) = get_user_real('vx (km/s): ', -20.0_DP, 20.0_DP)
                initial_guess(5) = get_user_real('vy (km/s): ', -20.0_DP, 20.0_DP)
                initial_guess(6) = get_user_real('vz (km/s): ', -20.0_DP, 20.0_DP)
                
            case (2)
                write(*, *) '请输入初始开普勒轨道根数猜测:'
                keplerian_elements(1) = get_user_real('半长轴 a (km): ', 1000.0_DP, 100000.0_DP)
                keplerian_elements(2) = get_user_real('偏心率 e: ', 0.0_DP, 0.99_DP)
                keplerian_elements(3) = get_user_real('轨道倾角 i (度): ', 0.0_DP, 180.0_DP)
                keplerian_elements(4) = get_user_real('升交点赤经 Ω (度): ', 0.0_DP, 360.0_DP)
                keplerian_elements(5) = get_user_real('近地点幅角 ω (度): ', 0.0_DP, 360.0_DP)
                keplerian_elements(6) = get_user_real('真近点角 ν (度): ', 0.0_DP, 360.0_DP)
                
                ! 转换为笛卡尔坐标
                call keplerian_to_cartesian(keplerian_elements, initial_guess)
        end select
        
        write(*, *) '初始轨道猜测:'
        call print_vector(initial_guess, '笛卡尔状态 (km, km/s)')
    end subroutine get_initial_orbit_guess
    
    subroutine select_improvement_method(method_choice)
        integer, intent(out) :: method_choice
        
        write(*, *) '选择轨道改进方法:'
        write(*, *) '1. 最小二乘法 (Least Squares)'
        write(*, *) '2. 扩展卡尔曼滤波 (Extended Kalman Filter)'
        write(*, *) '3. 无迹卡尔曼滤波 (Unscented Kalman Filter)'
        
        method_choice = get_user_choice('请选择方法 (1-3): ', 1, 3)
    end subroutine select_improvement_method
    
    subroutine improve_orbit(obs_data, initial_guess, method_choice, result)
        type(observation_data), intent(in) :: obs_data
        real(DP), dimension(6), intent(in) :: initial_guess
        integer, intent(in) :: method_choice
        type(orbit_improvement_result), intent(out) :: result
        
        real(DP), dimension(6) :: current_state, state_update
        real(DP), dimension(6,6) :: covariance
        real(DP), dimension(3) :: predicted_measurement, measurement_residual
        real(DP), dimension(3,6) :: jacobian
        real(DP) :: residual_norm, prev_residual_norm
        integer :: i, iteration
        
        ! 初始化
        result%initial_state = initial_guess
        current_state = initial_guess
        covariance = 0.0_DP
        do i = 1, 6
            covariance(i, i) = 1.0_DP  ! 单位矩阵作为初始协方差
        end do
        
        write(*, *) '开始轨道改进...'
        write(*, *) '观测数量: ', obs_data%n_obs
        write(*, *) '最大迭代次数: ', config%max_iterations
        
        ! 迭代改进
        prev_residual_norm = huge(1.0_DP)
        do iteration = 1, config%max_iterations
            residual_norm = 0.0_DP
            
            ! 对每个观测进行处理
            do i = 1, obs_data%n_obs
                ! 计算预测观测值
                call compute_measurement(current_state, obs_data%times(i), &
                                       obs_data%measurement_type, predicted_measurement)
                
                ! 计算观测残差
                measurement_residual = obs_data%measurements(i, :) - predicted_measurement
                
                ! 计算雅可比矩阵
                call compute_measurement_jacobian(current_state, obs_data%times(i), &
                                                obs_data%measurement_type, jacobian)
                
                ! 根据方法选择更新策略
                select case (method_choice)
                    case (1)
                        ! 最小二乘法
                        call least_squares_update(current_state, covariance, &
                                                measurement_residual, jacobian, &
                                                obs_data%measurement_weights(i, :), &
                                                state_update)
                    case (2)
                        ! 扩展卡尔曼滤波
                        call extended_kalman_update(current_state, covariance, &
                                                   measurement_residual, jacobian, &
                                                   obs_data%measurement_weights(i, :), &
                                                   state_update)
                    case (3)
                        ! 无迹卡尔曼滤波
                        call unscented_kalman_update(current_state, covariance, &
                                                    measurement_residual, jacobian, &
                                                    obs_data%measurement_weights(i, :), &
                                                    state_update)
                end select
                
                ! 更新状态
                current_state = current_state + state_update
                residual_norm = residual_norm + sum(measurement_residual**2)
            end do
            
            residual_norm = sqrt(residual_norm / obs_data%n_obs)
            
            ! 检查收敛性
            if (abs(residual_norm - prev_residual_norm) < config%convergence_tolerance) then
                result%converged = .true.
                exit
            end if
            
            prev_residual_norm = residual_norm
            
            ! 显示进度
            if (mod(iteration, 5) == 0) then
                write(*, '(A,I3,A,F12.6)') '迭代 ', iteration, ', 残差: ', residual_norm
            end if
        end do
        
        ! 保存结果
        result%improved_state = current_state
        result%covariance_matrix = covariance
        result%iterations = iteration
        result%final_residual = residual_norm
        result%convergence_criteria = config%convergence_tolerance
        
        if (iteration > config%max_iterations) then
            result%converged = .false.
            write(*, *) '警告: 达到最大迭代次数，可能未收敛'
        end if
        
        write(*, *) '轨道改进完成!'
    end subroutine improve_orbit
    
    subroutine least_squares_update(state, covariance, residual, jacobian, weights, update)
        real(DP), dimension(6), intent(in) :: state
        real(DP), dimension(6,6), intent(inout) :: covariance
        real(DP), dimension(3), intent(in) :: residual
        real(DP), dimension(3,6), intent(in) :: jacobian
        real(DP), dimension(3), intent(in) :: weights
        real(DP), dimension(6), intent(out) :: update
        
        real(DP), dimension(6,6) :: HtWH, HtWH_inv
        real(DP), dimension(6,3) :: HtW
        real(DP), dimension(3,3) :: W
        integer :: i
        
        ! 构建权重矩阵
        W = 0.0_DP
        do i = 1, 3
            W(i, i) = weights(i)
        end do
        
        ! 计算 H^T * W * H
        HtW = matmul(transpose(jacobian), W)
        HtWH = matmul(HtW, jacobian)
        
        ! 计算 (H^T * W * H)^(-1)
        call matrix_inverse(HtWH, HtWH_inv)
        
        ! 计算状态更新
        update = matmul(matmul(HtWH_inv, HtW), residual)
        
        ! 更新协方差矩阵
        covariance = HtWH_inv
    end subroutine least_squares_update
    
    subroutine extended_kalman_update(state, covariance, residual, jacobian, weights, update)
        real(DP), dimension(6), intent(in) :: state
        real(DP), dimension(6,6), intent(inout) :: covariance
        real(DP), dimension(3), intent(in) :: residual
        real(DP), dimension(3,6), intent(in) :: jacobian
        real(DP), dimension(3), intent(in) :: weights
        real(DP), dimension(6), intent(out) :: update
        
        real(DP), dimension(3,3) :: R, S
        real(DP), dimension(6,3) :: K
        real(DP), dimension(6,6) :: I
        
        ! 构建测量噪声协方差矩阵
        R = 0.0_DP
        R(1,1) = 1.0_DP / weights(1)
        R(2,2) = 1.0_DP / weights(2)
        R(3,3) = 1.0_DP / weights(3)
        
        ! 计算创新协方差
        S = matmul(matmul(jacobian, covariance), transpose(jacobian)) + R
        
        ! 计算卡尔曼增益
        K = matmul(matmul(covariance, transpose(jacobian)), S)
        
        ! 计算状态更新
        update = matmul(K, residual)
        
        ! 更新协方差矩阵
        I = 0.0_DP
        I(1, 1) = 1.0_DP
        I(2, 2) = 1.0_DP
        I(3, 3) = 1.0_DP
        I(4, 4) = 1.0_DP
        I(5, 5) = 1.0_DP
        I(6, 6) = 1.0_DP
        covariance = matmul(matmul(I - matmul(K, jacobian), covariance), &
                           transpose(I - matmul(K, jacobian))) + &
                    matmul(matmul(K, R), transpose(K))
    end subroutine extended_kalman_update
    
    subroutine unscented_kalman_update(state, covariance, residual, jacobian, weights, update)
        real(DP), dimension(6), intent(in) :: state
        real(DP), dimension(6,6), intent(inout) :: covariance
        real(DP), dimension(3), intent(in) :: residual
        real(DP), dimension(3,6), intent(in) :: jacobian
        real(DP), dimension(3), intent(in) :: weights
        real(DP), dimension(6), intent(out) :: update
        
        ! 简化实现，实际应用中需要更复杂的无迹变换
        call extended_kalman_update(state, covariance, residual, jacobian, weights, update)
    end subroutine unscented_kalman_update
    
    subroutine matrix_inverse(A, A_inv)
        real(DP), dimension(:,:), intent(in) :: A
        real(DP), dimension(:,:), intent(out) :: A_inv
        integer :: n, info, i
        real(DP), allocatable, dimension(:) :: work
        real(DP), allocatable, dimension(:,:) :: A_copy
        
        n = size(A, 1)
        allocate(A_copy(n, n), work(n))
        
        A_copy = A
        A_inv = 0.0_DP
        do i = 1, n
            A_inv(i, i) = 1.0_DP
        end do
        
        ! 简化实现，不使用LAPACK
        ! 这里应该实现矩阵求逆算法
        A_inv = A_copy  ! 临时简化
        
        deallocate(A_copy, work)
    end subroutine matrix_inverse
    
    subroutine display_improvement_results(result)
        type(orbit_improvement_result), intent(in) :: result
        real(DP), dimension(6) :: keplerian_initial, keplerian_improved
        
        write(*, *) '轨道改进结果:'
        write(*, *) '收敛状态: ', result%converged
        write(*, *) '迭代次数: ', result%iterations
        write(*, *) '最终残差: ', result%final_residual
        write(*, *) '收敛容差: ', result%convergence_criteria
        
        write(*, *) '初始状态:'
        call print_vector(result%initial_state, '笛卡尔状态 (km, km/s)')
        
        write(*, *) '改进后状态:'
        call print_vector(result%improved_state, '笛卡尔状态 (km, km/s)')
        
        ! 转换为开普勒轨道根数进行比较
        call cartesian_to_keplerian(result%initial_state, keplerian_initial)
        call cartesian_to_keplerian(result%improved_state, keplerian_improved)
        
        write(*, *) '初始开普勒轨道根数:'
        call print_vector(keplerian_initial, '开普勒根数 (km, -, 度)')
        
        write(*, *) '改进后开普勒轨道根数:'
        call print_vector(keplerian_improved, '开普勒根数 (km, -, 度)')
        
        write(*, *) '协方差矩阵:'
        call print_matrix(result%covariance_matrix, '状态协方差')
    end subroutine display_improvement_results
    
    subroutine save_improvement_results(result)
        type(orbit_improvement_result), intent(in) :: result
        character(len=MAX_STRING_LEN) :: filename
        
        ! 保存改进后的状态
        filename = trim(output_directory) // 'orbit_improvement_state.csv'
        call save_matrix_to_file(reshape(result%improved_state, [6, 1]), filename, '改进后轨道状态')
        
        ! 保存协方差矩阵
        filename = trim(output_directory) // 'orbit_improvement_covariance.csv'
        call save_matrix_to_file(result%covariance_matrix, filename, '状态协方差矩阵')
        
        write(*, *) '改进结果已保存到:'
        write(*, *) '  状态: ', trim(output_directory) // 'orbit_improvement_state.csv'
        write(*, *) '  协方差: ', trim(output_directory) // 'orbit_improvement_covariance.csv'
    end subroutine save_improvement_results
    
    subroutine cleanup_observation_data(obs_data)
        type(observation_data), intent(inout) :: obs_data
        
        if (allocated(obs_data%times)) deallocate(obs_data%times)
        if (allocated(obs_data%measurements)) deallocate(obs_data%measurements)
        if (allocated(obs_data%measurement_weights)) deallocate(obs_data%measurement_weights)
    end subroutine cleanup_observation_data

end module cat_orbit_improvement
