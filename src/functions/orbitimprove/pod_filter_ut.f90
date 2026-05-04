!-----------------------------------------------------------------------------------------
! 模块: pod_ut_filter_module
! 功能: UT 顶层滤波器集成
! 职责: 调度时间更新(UT)与测量更新(UT+卡尔曼融合+权重更新)
!-----------------------------------------------------------------------------------------
module pod_ut_filter_module
    use global, only: DP
    use pod_orbit_propagation_module, only: propagate_orbit,  orbit_state, propagation_result
    use pod_basicmath_module, only: inverse_and_determinant, dpotrf
    use pod_random_module, only: generate_multivariate_normal
    use pod_integrator_module, only: METHOD_RKF45, METHOD_RKF78
    use pod_measurement_model_module, only: compute_measurement
    implicit none

    private

    public :: ut_filter

    type :: ut_filter
        real(DP) :: current_epoch  ! 当前时间
        real(DP), allocatable :: state_mean(:)  ! 状态向量
        real(DP), allocatable :: state_cov(:,:)  ! 协方差矩阵
        real(DP), allocatable :: sigma_points(:,:)  ! Sigma 点矩阵
        real(DP), allocatable :: Weights(:)  ! Sigma 点权重

        
        real(DP), allocatable :: process_noise(:,:)  ! 过程噪声矩阵
        real(DP), allocatable :: measurement_noise(:,:)  ! 测量噪声矩阵

        real(DP) :: kappa  ! UT 参数

    contains
        procedure :: filter_time_update
        procedure :: filter_measurement_update
        procedure :: get_noise_vecs_from_noise
        procedure :: filter_init
        procedure :: calculate_sigma_points
    end type ut_filter

contains

    subroutine filter_init(this, initial_state, initial_cov, process_noise, measurement_noise, kappa)
        class(ut_filter), intent(inout) :: this
        real(DP), intent(in) :: initial_state(:)
        real(DP), intent(in) :: initial_cov(:,:)
        real(DP), intent(in) :: process_noise(:,:)
        real(DP), intent(in) :: measurement_noise(:,:)
        real(DP), intent(in) :: kappa

        this%state_mean = initial_state
        this%state_cov = initial_cov
        this%process_noise = process_noise
        this%measurement_noise = measurement_noise
        this%kappa = kappa
    end subroutine filter_init
    !> 时间更新：使用轨道传播器进行状态预测
    subroutine filter_time_update(this, et)
        class(ut_filter), intent(inout) :: this
        real(DP), intent(in) :: et  ! 更新时间

        integer :: dimension
        real(DP), allocatable :: sigma_points(:,:), Weights(:), propagated_points(:,:), noise_vec(:,:)

        type(orbit_state) :: initial_state
        type(propagation_result) :: result
        real(DP) :: t_end

        initilal_state%state = this%state_mean
        initial_state%epoch = this%current_epoch

        t_end = et - this%current_epoch

        dimension = size(this%state_mean)
        allocate(sigma_points(dimension, 2*dimension+1), propagated_points(dimension, 2*dimension+1))
        allocate(Weights(2*dimension+1))
        allocate(noise_vec(dimension, 2*dimension+1))
        ! 权重
        Weights(1) = this%kappa / (dimension + this%kappa)
        Weights(2:2*dimension+1) = 1.0_DP / (2.0_DP * (dimension + this%kappa))
        ! 根据状态向量以及协方差矩阵计算sigma点：
        call calculate_sigma_points(this%state_mean, this%state_cov, this%kappa, sigma_points)
        

        ! 调用轨道传播器进行状态预测
        !! 增加过噪
        call get_noise_vecs_from_noise(this, 2*dimension+1, noise_vec)
        
        for i = 1, size(sigma_points, 2)
            call propagate_orbit(sigma_points(:, i), t_end, METHOD_RKF78, result)
            propagated_points(:, i) = result%states(:, end) + noise_vec(:,i) ! 使用传播结果的最终状态
        end do
       
        ! 根据传播结果计算预测的状态均值和协方差
        this%state_mean = 0.0_DP
        for i = 1, size(sigma_points, 2)
            this%state_mean = this%state_mean + Weights(i) * propagated_points(:, i)  ! 使用传播结果的最终状态
        end do    
        
        this%state_cov = 0.0_DP
        for i = 1, size(sigma_points, 2)
            this%state_cov = this%state_cov + Weights(i) * matmul(reshape(propagated_points(:, i) - &
            this%state_mean, (/dimension, 1/)), reshape(propagated_points(:, i) - this%state_mean, (/1, dimension/)))
        end do

        ! 文章未增加过程噪声(Julier), 实际上UKF也不需要
        
        ! 更新时间
        this%current_epoch = et
        ! 更新sigma点
        this%sigma_points = propagated_points
        this%Weights = Weights

        ! 释放临时变量
        deallocate(sigma_points, propagated_points, Weights)

    end subroutine filter_time_update

    !> 测量更新：融合测量信息，更新状态和协方差
    subroutine filter_measurement_update(this, y_meas, et, station)
        class(ut_filter), intent(inout) :: this
        real(DP), dimension(:), intent(in) :: y_meas  ! 真实测量值
        real(DP), intent(in) :: et  ! 测量时间
        type(observation_station), intent(in) :: station  ! 测站信息

        real(DP), allocatable :: points_z(:,:) ! sigma点在测量空间的映射
        real(DP), allocatable :: predicted_measurement(:) ! 预测测量值
        
        real(DP), allocatable :: P_xz(:,:), P_zz(:,:), P_zz_inv(:,:)
        real(DP) :: det_Pzz
        real(DP), allocatable :: Kalman_gain(:,:)
        integer :: i, measurement_dim, state_dim

        measurement_dim = size(y_meas)
        state_dim = size(this%state_mean)

        allocate(points_z(measurement_dim, 2*state_dim+1))
        allocate(predicted_measurement(measurement_dim))

        allocate(P_xz(state_dim, measurement_dim), P_zz(measurement_dim, measurement_dim), &
        P_zz_inv(measurement_dim, measurement_dim))
        ! 1. 计算预测测量值
        predicted_measurement = 0.0_DP
        for i = 1, size(this%sigma_points, 2)
            points_z(:, i) = compute_measurement(this%sigma_points(:, i), et, station, 'OPTICAL') 
            predicted_measurement = predicted_measurement + this%Weights(i) * points_z(:, i) 
        end do

        ! 2. 计算协方差矩阵 P_xz 和 P_zz
        P_xz = 0.0_DP
        P_zz = 0.0_DP
        for i = 1, size(this%sigma_points, 2)
            P_xz = P_xz + this%Weights(i) * matmul(reshape(this%sigma_points(:, i) - this%state_mean, (/state_dim, 1/)), 
            &reshape(points_z(:, i) - predicted_measurement, (/1, measurement_dim/)))
            P_zz = P_zz + this%Weights(i) * matmul(reshape(points_z(:, i) - predicted_measurement, (/measurement_dim, 1/)),
            & reshape(points_z(:, i) - predicted_measurement, (/1, measurement_dim/)))
        end do
        ! 加上测量噪声协方差
        P_zz = P_zz + this%measurement_noise

        ! 3. 计算卡尔曼增益
        call inverse_and_determinant(P_zz, P_zz_inv, det_Pzz)
        Kalman_gain = matmul(P_xz, P_zz_inv)

        ! 4. 更新状态均值和协方差
        this%state_mean = this%state_mean + matmul(Kalman_gain, reshape(y_meas - predicted_measurement, (/measurement_dim, 1/)))
        this%state_cov = this%state_cov - matmul(Kalman_gain, matmul(P_zz, transpose(Kalman_gain)))

        ! 释放临时变量
        deallocate(points_z, predicted_measurement, P_xz, P_zz, P_zz_inv, Kalman_gain)
    end subroutine filter_measurement_update


    subroutine get_noise_vecs_from_noise(this, n_samples, addos_out)
        class(ut_filter), intent(in) :: this
        integer, intent(in) :: n_samples
        real(DP), dimension(:,:), intent(out) :: addos_out
        
        call generate_multivariate_normal(0.0_DP, this%process_noise, addos_out)
    end subroutine get_noise_vecs_from_noise

    subroutine calculate_sigma_points(mean, cov, kappa, sigma_points)
        real(DP), intent(in) :: mean(:)
        real(DP), intent(in) :: cov(:,:)
        real(DP), intent(in) :: kappa
        real(DP), intent(out) :: sigma_points(:,:)

        integer :: dimension
        dimension = size(mean)

        ! 计算协方差矩阵的平方根（Cholesky分解）
        real(DP), allocatable :: sqrt_cov(:,:)
        allocate(sqrt_cov(dimension, dimension))
        call dpotrf('L', dimension, sqrt_cov, dimension, info)

        ! 生成sigma点
        sigma_points(:,1) = mean
        do i = 1, dimension
            sigma_points(:, i+1) = mean + sqrt(dimension + kappa) * sqrt_cov(:, i)
            sigma_points(:, i+1+dimension) = mean - sqrt(dimension + kappa) * sqrt_cov(:, i)
        end do
    end subroutine calculate_sigma_points


end module pod_ut_filter_module

