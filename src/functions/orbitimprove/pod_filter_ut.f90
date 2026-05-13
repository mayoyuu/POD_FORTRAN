!-----------------------------------------------------------------------------------------
! 模块: pod_filter_ut_module
! 功能: UT 顶层滤波器集成
! 职责: 调度时间更新(UT)与测量更新(UT+卡尔曼融合+权重更新)
! cite : The unscented Kalman filter for nonlinear estimation 2000
!-----------------------------------------------------------------------------------------
module pod_filter_ut_module
    use pod_global, only: DP
    use pod_orbit_propagation, only: propagate_orbit,  orbit_state, propagation_result
    use pod_basicmath_module, only: inverse_and_determinant, dpotrf, PI
    use pod_random_module, only: generate_multivariate_normal
    use pod_integrator_module, only: METHOD_RKF45, METHOD_RKF78
    use pod_measurement_model_module, only: compute_measurement
    use pod_measurement_base_module, only: observation_station 
    implicit none

    private

    public :: ut_filter

    type :: ut_filter
        real(DP) :: current_epoch  ! 当前时间
        real(DP), allocatable :: state_mean(:)  ! 状态向量
        real(DP), allocatable :: state_cov(:,:)  ! 协方差矩阵
        real(DP), allocatable :: sigma_points(:,:)  ! Sigma 点矩阵
        ! real(DP), allocatable :: Weights(:)  ! Sigma 点权重

        ! real(DP) :: kappa = 0.5_DP ! UT 参数
        ! 修改：将 Weights 拆分为均值权重和协方差权重
        real(DP), allocatable :: Weights_m(:)  ! Sigma 点均值权重
        real(DP), allocatable :: Weights_c(:)  ! Sigma 点协方差权重

        ! 修改：引入 van der Merwe 参数
        real(DP) :: alpha = 1.0e-3_DP
        real(DP) :: beta  = 2.0_DP
        real(DP) :: kappa = 0.0_DP

        real(DP) :: last_residual(6)    ! 存储最近一次测量更新的 6 列残差
        real(DP) :: last_comp_val(2)   ! 存储最近一次计算出的预测观测值 (Lon, Lat)

    contains
        procedure :: time_update => filter_time_update
        procedure :: measurement_update => filter_measurement_update
        
        procedure :: filter_init
        ! procedure :: set_kappa => filter_set_kappa
        procedure :: set_epoch => filter_set_epoch
        procedure :: get_current_epoch => filter_get_current_epoch
        procedure :: get_current_state => filter_get_current_state
        procedure :: get_current_cov => filter_get_current_cov
        procedure :: get_last_residual => filter_get_last_residual

    end type ut_filter

contains

    subroutine filter_init(this, initial_epoch, initial_state, initial_cov)
        class(ut_filter), intent(inout) :: this
        real(DP), intent(in) :: initial_epoch
        real(DP), intent(in) :: initial_state(:)
        real(DP), intent(in) :: initial_cov(:,:)

        this%current_epoch = initial_epoch
        this%state_mean = initial_state
        this%state_cov = initial_cov
    end subroutine filter_init

    subroutine filter_set_epoch(this, epoch)
        class(ut_filter), intent(inout) :: this
        real(DP), intent(in) :: epoch
        this%current_epoch = epoch
    end subroutine filter_set_epoch

    subroutine filter_get_current_epoch(this, epoch_out)
        class(ut_filter), intent(in) :: this
        real(DP), intent(out) :: epoch_out
        epoch_out = this%current_epoch
    end subroutine filter_get_current_epoch

    subroutine filter_get_current_state(this, mean_out)
        class(ut_filter), intent(in) :: this
        real(DP), dimension(:), intent(out) :: mean_out
        mean_out = this%state_mean
    end subroutine filter_get_current_state

    subroutine filter_get_current_cov(this, cov_out)
        class(ut_filter), intent(in) :: this
        real(DP), dimension(:,:), intent(out) :: cov_out
        cov_out = this%state_cov
    end subroutine filter_get_current_cov

        !> 获取单步残差的接口
    subroutine filter_get_last_residual(this, res_out, comp_out)
        class(ut_filter), intent(in) :: this
        real(DP), intent(out) :: res_out(6), comp_out(2)
        res_out = this%last_residual
        comp_out = this%last_comp_val
    end subroutine filter_get_last_residual

    !> 时间更新：使用轨道传播器进行状态预测
    subroutine filter_time_update(this, et, noise_Q)
        class(ut_filter), intent(inout) :: this
        real(DP), intent(in) :: et  ! 更新时间
        real(DP), dimension(:,:), intent(in), optional:: noise_Q    ! 过程噪声协方差 Q

        integer :: dimension,i
        real(DP) :: lambda
        real(DP), allocatable :: sigma_points(:,:), Weights_m(:), Weights_c(:), propagated_points(:,:)
        
        type(orbit_state) :: initial_state
        type(propagation_result) :: result
        real(DP) :: t_end

        ! allocate(Weights(2*dimension+1))
        ! ! 权重
        ! Weights(1) = this%kappa / (real(dimension, DP) + this%kappa)
        ! Weights(2:2*dimension+1) = 1.0_DP / (2.0_DP * (real(dimension, DP) + this%kappa))
        
        ! ... 分配内存
        dimension = size(this%state_mean)
        allocate(sigma_points(dimension, 2*dimension+1), propagated_points(dimension, 2*dimension+1))
        allocate(Weights_m(2*dimension+1), Weights_c(2*dimension+1))

        !! ... 权重计算
        !   计算lambda
        lambda = (this%alpha**2) * (real(dimension, DP) + this%kappa) - real(dimension, DP)
        ! 权重计算 (Van der Merwe)
        Weights_m(1) = lambda / (real(dimension, DP) + lambda)
        Weights_c(1) = lambda / (real(dimension, DP) + lambda) + (1.0_DP - this%alpha**2 + this%beta)
        
        Weights_m(2:2*dimension+1) = 1.0_DP / (2.0_DP * (real(dimension, DP) + lambda))
        Weights_c(2:2*dimension+1) = Weights_m(2:2*dimension+1)

        !! 生成sigma点
        ! 将 lambda传给 sigma 点生成函数
        call calculate_sigma_points(this%state_mean, this%state_cov, lambda, sigma_points)


        initial_state%epoch = this%current_epoch
        t_end = et - this%current_epoch

        ! 调用轨道传播器进行状态预测
        do i = 1, size(sigma_points, 2)
            initial_state%state = sigma_points(:, i) 
            call propagate_orbit(initial_state, t_end, METHOD_RKF78, result)
            propagated_points(:, i) = result%states(result%n_steps, :) ! 使用传播结果的最终状态
        end do
       
        ! 根据传播结果计算预测的状态均值和协方差
        this%state_mean = 0.0_DP
        do i = 1, size(sigma_points, 2)
            this%state_mean = this%state_mean + Weights_m(i) * propagated_points(:, i)  ! 使用传播结果的最终状态
        end do    
        
        this%state_cov = 0.0_DP
        do i = 1, size(sigma_points, 2)
            this%state_cov = this%state_cov + Weights_c(i) * matmul(reshape(propagated_points(:, i) - &
            this%state_mean, (/dimension, 1/)), reshape(propagated_points(:, i) - this%state_mean, (/1, dimension/)))
        end do
        if(present(noise_Q)) this%state_cov = this%state_cov + noise_Q
        this%state_cov = 0.5_DP * (this%state_cov + transpose(this%state_cov))
        
        ! 更新时间
        this%current_epoch = et
        ! 更新sigma点
        ! call calculate_sigma_points(this%state_mean, this%state_cov, lambda, propagated_points)
        this%sigma_points = propagated_points
        this%Weights_m = Weights_m
        this%Weights_c = Weights_c

        ! 释放临时变量
        deallocate(sigma_points, propagated_points, Weights_m, Weights_c)

    end subroutine filter_time_update

    !> 测量更新：融合测量信息，更新状态和协方差
    subroutine filter_measurement_update(this, y_meas, noise_R, et, station)
        class(ut_filter), intent(inout) :: this
        real(DP), dimension(:), intent(in) :: y_meas  ! 真实测量值
        real(DP), dimension(:,:), intent(in):: noise_R    ! 测量噪声协方差 R
        real(DP), intent(in) :: et  ! 测量时间
        type(observation_station), intent(in) :: station  ! 测站信息

        real(DP), allocatable :: points_z(:,:) ! sigma点在测量空间的映射
        real(DP), allocatable :: predicted_measurement(:) ! 预测测量值
        
        real(DP), allocatable :: P_xz(:,:), P_zz(:,:), P_zz_inv(:,:)
        real(DP) :: det_Pzz
        real(DP), allocatable :: Kalman_gain(:,:)
        integer :: i, info, measurement_dim, state_dim
        real(DP) :: pred_z_pre(2)

        measurement_dim = size(y_meas)
        state_dim = size(this%state_mean)

        allocate(points_z(measurement_dim, 2*state_dim+1))
        allocate(predicted_measurement(measurement_dim))

        allocate(P_xz(state_dim, measurement_dim), P_zz(measurement_dim, measurement_dim), &
        P_zz_inv(measurement_dim, measurement_dim))

        ! ==========================================================
        ! 【核心新增】：计算更新前的残差 (Innovation) 用于输出报告
        ! ==========================================================
        call compute_measurement(this%state_mean, et, station, 'OPTICAL',pred_z_pre) 
        
        this%last_comp_val = pred_z_pre(1:2)
        this%last_residual = 0.0_DP
        ! 1. dA * cos(Dec) 单位：arcsec
        this%last_residual(1) = (y_meas(1) - pred_z_pre(1)) * 3600.0_DP &
        * cos(y_meas(2) * PI / 180.0_DP)
        ! 2. dDec 单位：arcsec
        this%last_residual(2) = (y_meas(2) - pred_z_pre(2)) * 3600.0_DP
        ! 3. Total Angle Residual
        this%last_residual(3) = sqrt(this%last_residual(1)**2 + this%last_residual(2)**2)
        ! 4-6. 占位符 (未来可扩展 RTN 残差)
        this%last_residual(4:6) = 0.0_DP
        ! 1. 计算预测测量值
        predicted_measurement = 0.0_DP

        do i = 1, size(this%sigma_points, 2)
            call compute_measurement(this%sigma_points(:, i), et, station, 'OPTICAL',points_z(:,i)) 
            predicted_measurement = predicted_measurement + this%Weights_m(i) * points_z(:, i) 
        end do

        ! 2. 计算协方差矩阵 P_xz 和 P_zz
        P_xz = 0.0_DP
        P_zz = 0.0_DP
        do i = 1, size(this%sigma_points, 2)
            P_xz = P_xz + this%Weights_c(i) * matmul(reshape(this%sigma_points(:, i) - this%state_mean, (/state_dim, 1/)), &
            reshape(points_z(:, i) - predicted_measurement, (/1, measurement_dim/)))
            P_zz = P_zz + this%Weights_c(i) * matmul(reshape(points_z(:, i) - predicted_measurement, (/measurement_dim, 1/)),& 
            reshape(points_z(:, i) - predicted_measurement, (/1, measurement_dim/)))
        end do
        ! 加上测量噪声协方差
        P_zz = P_zz + noise_R

        ! 3. 计算卡尔曼增益
        call inverse_and_determinant(P_zz, P_zz_inv, det_Pzz, info)
        Kalman_gain = matmul(P_xz, P_zz_inv)

        ! 4. 更新状态均值和协方差
        this%state_mean = this%state_mean + matmul(Kalman_gain, y_meas - predicted_measurement)
        this%state_cov = this%state_cov - matmul(Kalman_gain, matmul(P_zz, transpose(Kalman_gain)))
        this%state_cov = 0.5_DP * (this%state_cov + transpose(this%state_cov))

        ! 释放临时变量
        deallocate(points_z, predicted_measurement, P_xz, P_zz, P_zz_inv, Kalman_gain)
    end subroutine filter_measurement_update


    subroutine calculate_sigma_points(mean, cov, lambda, sigma_points)
        real(DP), intent(in) :: mean(:)
        real(DP), intent(in) :: cov(:,:)
        real(DP), intent(in) :: lambda
        real(DP), intent(out) :: sigma_points(:,:)
        integer :: i,ii, info
        integer :: dimension
        ! 计算协方差矩阵的平方根（Cholesky分解）
        real(DP), allocatable :: sqrt_cov(:,:)

        dimension = size(mean)
        allocate(sqrt_cov(dimension, dimension))
        sqrt_cov = cov
        call dpotrf('L', dimension, sqrt_cov, dimension, info)
        if (info /= 0) then
            ! 协方差不正定，加微小对角扰动后重试
            do ii = 1, dimension
                sqrt_cov(ii,ii) = cov(ii,ii) + 1.0e-10_DP
            end do
            call dpotrf('L', dimension, sqrt_cov, dimension, info)
            if (info /= 0) then
                write(*,*) 'FATAL: Covariance still not PD after regularization'
                return
            end if
        end if

        ! 生成sigma点
        sigma_points(:,1) = mean
        do i = 1, dimension
            sigma_points(:, i+1) = mean + sqrt(real(dimension, DP) + lambda) * sqrt_cov(:, i)
            sigma_points(:, i+1+dimension) = mean - sqrt(real(dimension, DP) + lambda) * sqrt_cov(:, i)
        end do
    end subroutine calculate_sigma_points


end module pod_filter_ut_module

