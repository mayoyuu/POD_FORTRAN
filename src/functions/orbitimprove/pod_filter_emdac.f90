!-----------------------------------------------------------------------------------------
! 模块: pod_filter_emdac_module
! 功能: EMDAC-N 顶层滤波器集成
! 职责: 调度时间更新(粒子传播+EM聚类)与测量更新(UT/DA映射+卡尔曼融合+权重更新)
!-----------------------------------------------------------------------------------------
module pod_filter_emdac_module
    use pod_global, only: DP
    use pod_uq_gmm_state_module, only: uq_gmm_state_type
    use pod_uq_state_module, only: uq_state_type
    use pod_uq_prop_da_module, only: uq_da_propagator
    use pod_gmm_math_module, only: fit_gmm_to_particles
    use pod_uq_transform_ut_module
    use pod_uq_transform_da_module
    use pod_measurement_model_module, only: observation_station
    use pod_basicmath_module, only: inverse_and_determinant

    implicit none
    private
    
    public :: emdac_filter

    !> ======================================================================
    !> 滤波器主类定义
    !> ======================================================================
    type :: emdac_filter
        type(uq_gmm_state_type) :: gmm_state          ! 当前 k 时刻的 GMM 状态
        integer                 :: n_particles = 10000! 粒子总数
        real(DP)                :: em_tol = 1.0e-4_DP ! EM 算法收敛容差
        integer                 :: em_max_iter = 50   ! EM 算法最大迭代次数
        
        type(uq_da_propagator)  :: propagator         ! 绑定的 DA 时间传播器
    contains
        procedure :: init => filter_init
        procedure :: time_update => filter_time_update
        procedure :: measurement_update => filter_measurement_update
    end type emdac_filter

contains

    !> 1. 初始化滤波器
    subroutine filter_init(this, initial_mean, initial_cov, n_comp, n_part)
        class(emdac_filter), intent(inout) :: this
        real(DP), intent(in) :: initial_mean(:), initial_cov(:,:)
        integer, intent(in)  :: n_comp, n_part
        
        this%n_particles = n_part
        ! 利用单高斯分布初始化最初的 GMM
        call this%gmm_state%init_from_single(n_comp, initial_mean, initial_cov)
    end subroutine filter_init


    !> ======================================================================
    !> 2. 时间更新 (Time Update)
    !> 逻辑: k步GMM -> 按权重采样粒子 -> DA动力学传播 -> EM聚类 -> k+1步GMM
    !> ======================================================================
    subroutine filter_time_update(this, t_start, t_end)
        class(emdac_filter), intent(inout) :: this
        real(DP), intent(in) :: t_start, t_end
        
        type(uq_state_type) :: uq_in, uq_out
        integer :: dim
        
        dim = this%gmm_state%state_dim
        call uq_in%allocate_memory(dim, this%n_particles)
        
        ! 步骤 A: 桥接参数空间与粒子空间 —— 根据 GMM 权重生成散布粒子
        call sample_particles_from_gmm(this%gmm_state, this%n_particles, uq_in%samples)
        
        ! 步骤 B: 调用 DA 传播器进行高度非线性的动力学积分
        write(*,*) '[EMDAC] 执行时间更新: ', t_start, ' -> ', t_end
        call this%propagator%propagate(t_start, t_end, uq_in, uq_out)
        
        ! 步骤 C: 重新执行 K-means + EM 算法，提取传播后的 GMM 拓扑结构
        write(*,*) '[EMDAC] 执行高斯混合聚类...'
        call fit_gmm_to_particles(uq_out%samples, this%gmm_state, this%em_max_iter, this%em_tol)
        
    end subroutine filter_time_update


    !> ======================================================================
    !> 3. 测量更新 (Measurement Update)
    !> 逻辑: UT/DA 变换 -> 卡尔曼状态更新 -> 似然计算 -> GMM 权重更新
    !> ======================================================================
    subroutine filter_measurement_update(this, y_meas, noise_R, et, station, is_angle, use_da)
        class(emdac_filter), intent(inout) :: this
        real(DP), dimension(:), intent(in) :: y_meas      ! 真实测量值 y_{k+1}
        real(DP), dimension(:,:), intent(in):: noise_R    ! 测量噪声协方差 R
        real(DP), intent(in)               :: et          ! 历元
        type(observation_station), intent(in):: station   ! 测站信息
        logical, dimension(:), intent(in)  :: is_angle    ! 是否为角度测量 (用于边界平滑)
        logical, intent(in), optional      :: use_da      ! 控制开关
        
        integer :: i, n_comp, ny, info
        real(DP), allocatable :: mean_z(:), P_zz(:,:), P_xz(:,:), K_gain(:,:), P_zz_inv(:,:)
        real(DP), allocatable :: innovation(:), log_likelihood(:)
        real(DP) :: det_Pzz, max_ll, sum_exp, mahalanobis_sq
        logical :: do_da
        real(DP) :: PI = 3.14159265358979323846_DP
        
        do_da = .false.
        if (present(use_da)) do_da = use_da
        
        n_comp = this%gmm_state%n_components
        ny = size(y_meas)
        allocate(log_likelihood(n_comp))
        
        write(*,*) '[EMDAC] 执行测量更新...'
        
        ! 遍历 GMM 中的每一个高斯核
        do i = 1, n_comp
            ! ---------------------------------------------------------
            ! A. 空间映射: 获取预测均值 \hat{z}，预测协方差 P_{zz}，交叉协方差 P_{xz}
            ! ---------------------------------------------------------
            if (do_da) then
                ! TODO: 插入你的 DA 变换调用链 (init_da_expansion -> evaluate_da -> reconstruct_da)
                ! call reconstruct_da_moments(..., mean_z, P_zz, P_xz, is_angle)
            else
                ! TODO: 插入你的 UT 变换调用链 (generate_sigma_points -> 物理方程 -> reconstruct_ut)
                ! call reconstruct_ut_moments(..., mean_z, P_zz, P_xz, is_angle)
            end if
            
            ! ---------------------------------------------------------
            ! B. 标准卡尔曼方程更新
            ! ---------------------------------------------------------
            ! 加入测量噪声
            P_zz = P_zz + noise_R
            
            call inverse_and_determinant(P_zz, P_zz_inv, det_Pzz, info)
            if (info /= 0) write(*,*) '[警告] P_zz 求逆失败！'
            
            ! 计算卡尔曼增益 K = P_{xz} * (P_{zz})^{-1}
            K_gain = matmul(P_xz, P_zz_inv)
            
            ! 计算残差 (Innovation)
            innovation = y_meas - mean_z
            ! [极其关键] 角度残差修正，防止 359度 - 1度 = 358度 的谬误
            call fix_angle_innovation(innovation, is_angle) 
            
            ! 更新当前核的均值和协方差 (Eq. 30, Eq. 31)
            this%gmm_state%components(i)%mean = this%gmm_state%components(i)%mean + matmul(K_gain, innovation)
            this%gmm_state%components(i)%cov  = this%gmm_state%components(i)%cov  - matmul(K_gain, matmul(P_zz, transpose(K_gain)))
            
            ! ---------------------------------------------------------
            ! C. 计算更新权重所需的对数似然概率 (Log-Likelihood)
            ! 公式: \ln(\mathcal{N}(y | \hat{z}, P_{zz})) + \ln(\mu_i^-)
            ! ---------------------------------------------------------
            mahalanobis_sq = dot_product(innovation, matmul(P_zz_inv, innovation))
            log_likelihood(i) = log(this%gmm_state%components(i)%weight) &
                              - 0.5_DP * ny * log(2.0_DP * PI) &
                              - 0.5_DP * log(det_Pzz) &
                              - 0.5_DP * mahalanobis_sq
        end do
        
        ! ---------------------------------------------------------
        ! D. GMM 权重软更新 (使用 Log-Sum-Exp 技巧防止下溢)
        ! ---------------------------------------------------------
        max_ll = maxval(log_likelihood)
        sum_exp = 0.0_DP
        do i = 1, n_comp
            sum_exp = sum_exp + exp(log_likelihood(i) - max_ll)
        end do
        
        ! 归一化写回
        do i = 1, n_comp
            this%gmm_state%components(i)%weight = exp(log_likelihood(i) - max_ll - log(sum_exp))
        end do
        
        deallocate(log_likelihood)
    end subroutine filter_measurement_update

    !> ======================================================================
    !> 辅助子程序: 处理光学观测的残差跨界问题
    !> ======================================================================
    subroutine fix_angle_innovation(innovation, is_angle)
        real(DP), dimension(:), intent(inout) :: innovation
        logical, dimension(:), intent(in)     :: is_angle
        integer :: r
        real(DP) :: PI = 3.14159265358979323846_DP
        
        do r = 1, size(innovation)
            if (is_angle(r)) then
                if (innovation(r) > PI)  innovation(r) = innovation(r) - 2.0_DP * PI
                if (innovation(r) < -PI) innovation(r) = innovation(r) + 2.0_DP * PI
            end if
        end do
    end subroutine fix_angle_innovation

end module pod_filter_emdac_module