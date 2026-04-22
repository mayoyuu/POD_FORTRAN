!-----------------------------------------------------------------------------------------
! 模块: pod_filter_emdac_module
! 功能: EMDAC-N 顶层滤波器集成
! 职责: 调度时间更新(粒子传播+EM聚类)与测量更新(UT/DA映射+卡尔曼融合+权重更新)
!-----------------------------------------------------------------------------------------
module pod_filter_emdac_module
    use pod_global, only: DP
    use pod_dace_classes
    use pod_uq_gmm_state_module, only: uq_gmm_state_type
    use pod_uq_state_module, only: uq_state_type
    use pod_gmm_math_module, only: fit_gmm_to_particles
    use pod_uq_propagation, only: run_particle_propagation
    use pod_measurement_base_module, only: observation_station 
    use pod_measurement_da_module, only: compute_measurement_da
    use pod_basicmath_module, only: inverse_and_determinant
     use pod_random_module, only: generate_multivariate_normal

    implicit none
    private
    
    public :: emdac_filter

    !> ======================================================================
    !> 滤波器主类定义
    !> ======================================================================
    type :: emdac_filter
        private
        real(DP) :: current_epoch ! 当前历元
        real(DP), allocatable :: state_mean(:)          ! 当前 GMM 的全局均值 (动态分配)  
        real(DP), allocatable :: state_cov(:,:)         ! 当前 GMM 的全局协方差 (动态分配)     
        type(uq_gmm_state_type), public :: gmm_state            ! 当前 k 时刻的 GMM 状态

        integer                 :: n_particles = 10000 ! 粒子总数
        real(DP)                :: em_tol = 1.0e-4_DP   ! EM 算法收敛容差
        integer                 :: em_max_iter = 50     ! EM 算法最大迭代次数
        integer                 :: da_order = 2         ! DA 阶数

        type(uq_state_type) :: propagated_particles   ! 保存时间更新后的粒子
        real(DP), allocatable :: current_omega(:,:)   ! 保存 EM 聚类后的责任度矩阵
        real(DP), allocatable :: current_W(:)         ! 保存各核的总责任度
        
    contains

        procedure :: init => filter_init

        procedure :: set_da_order => filter_set_da_order
        procedure :: set_em_parameters => filter_set_em_parameters
        procedure :: set_n_particles => filter_set_n_particles
        
        procedure :: time_update => filter_time_update
        procedure :: measurement_update => filter_measurement_update

        procedure :: sample_particles_from_gmm => sample_particles_from_gmm
        procedure, private :: update_global_mean, update_global_cov

        procedure :: get_current_epoch => filter_get_current_epoch
        procedure :: get_current_state => filter_get_current_state
        procedure :: get_current_cov => filter_get_current_cov
        procedure :: get_current_gmm => filter_get_current_gmm
    end type emdac_filter

contains

    !> 初始化函数，一次性注入所有配置并确立初始历元
    !> 初始化函数，一次性注入所有配置并确立初始历元
    subroutine filter_init(this, initial_epoch, initial_mean, initial_cov, n_comp, &
                           initial_gmm, n_part, opt_da_order, opt_em_tol, opt_em_max_iter)
        class(emdac_filter), intent(inout) :: this
        real(DP), intent(in) :: initial_epoch
        real(DP), intent(in) :: initial_mean(:), initial_cov(:,:)
        integer, intent(in)  :: n_comp  ! 必传参数全部挪到前面
        
        ! 所有的 optional 参数统一放到后面
        type(uq_gmm_state_type), intent(in), optional :: initial_gmm
        integer, intent(in), optional :: n_part, opt_da_order, opt_em_max_iter
        real(DP), intent(in), optional :: opt_em_tol

        ! 基础必传赋值
        this%current_epoch = initial_epoch
        this%state_mean = initial_mean
        this%state_cov = initial_cov

        ! 安全的可选参数赋值：只有外部传了，才覆盖默认值
        if (present(n_part)) this%n_particles = n_part
        if (present(opt_da_order)) this%da_order = opt_da_order
        if (present(opt_em_tol)) this%em_tol = opt_em_tol
        if (present(opt_em_max_iter)) this%em_max_iter = opt_em_max_iter

        ! 初始化 GMM 状态
        if (present(initial_gmm)) then
            this%gmm_state = initial_gmm
        else 
            call this%gmm_state%init_from_single(n_comp, initial_mean, initial_cov)
        end if

        write(*,*) '[EMDAC] 滤波器初始化完成！'
        write(*,*) '  初始历元: ', this%current_epoch
        write(*,*) '  粒子总数: ', this%n_particles
        write(*,*) '  DA 阶数 : ', this%da_order
        write(*,*) '  EM 最大迭代次数: ', this%em_max_iter
        write(*,*) '  EM 收敛容差: ', this%em_tol
    end subroutine filter_init

    subroutine filter_set_da_order(this, order)
        class(emdac_filter), intent(inout) :: this
        integer, intent(in) :: order
        this%da_order = order
    end subroutine filter_set_da_order

    subroutine filter_set_em_parameters(this, max_iter, tol)
        class(emdac_filter), intent(inout) :: this
        integer, intent(in), optional :: max_iter
        real(DP), intent(in), optional :: tol
        if (present(max_iter)) this%em_max_iter = max_iter
        if (present(tol)) this%em_tol = tol
    end subroutine filter_set_em_parameters

    subroutine filter_set_n_particles(this, n_part)
        class(emdac_filter), intent(inout) :: this
        integer, intent(in) :: n_part
        this%n_particles = n_part
    end subroutine filter_set_n_particles

    subroutine filter_get_current_epoch(this, epoch_out)
        class(emdac_filter), intent(in) :: this
        real(DP), intent(out) :: epoch_out
        epoch_out = this%current_epoch
    end subroutine filter_get_current_epoch

    subroutine filter_get_current_state(this, mean_out)
        class(emdac_filter), intent(in) :: this
        real(DP), dimension(:), intent(out) :: mean_out
        mean_out = this%state_mean
    end subroutine filter_get_current_state

    subroutine filter_get_current_cov(this, cov_out)
        class(emdac_filter), intent(in) :: this
        real(DP), dimension(:,:), intent(out) :: cov_out
        cov_out = this%state_cov
    end subroutine filter_get_current_cov

    subroutine filter_get_current_gmm(this, gmm_out)
        class(emdac_filter), intent(in) :: this
        type(uq_gmm_state_type), intent(out) :: gmm_out
        gmm_out = this%gmm_state
    end subroutine filter_get_current_gmm


    !> ======================================================================
    !> 2. 时间更新 (Time Update)
    !> 逻辑: k步GMM -> 按权重采样粒子 -> DA动力学传播 -> EM聚类 -> k+1步GMM
    !> ======================================================================
    subroutine filter_time_update(this, et)
        class(emdac_filter), intent(inout) :: this
        real(DP), intent(in) :: et  ! 目标历元 (绝对时间)
        type(uq_state_type) :: uq_in
        real(DP) :: t_end
        integer :: dim
        
        dim = this%gmm_state%state_dim
        call uq_in%allocate_memory(dim, this%n_particles)
        
        ! 步骤 A: 桥接参数空间与粒子空间 —— 根据 GMM 权重生成散布粒子
        call sample_particles_from_gmm(this%gmm_state, this%n_particles, uq_in%samples)

        t_end = et - this%current_epoch
        ! 步骤 B: 调用 DA 传播器进行高度非线性的动力学积分
        call run_particle_propagation(uq_in, this%state_mean, this%current_epoch, &
                                      0.0, t_end, METHOD_DA,  propagated_particles, &
                                      da_order=this%da_order, &
                                      reference_orbit_out=this%state_mean) 
        
        ! 步骤 C: 重新执行 K-means + EM 算法，提取传播后的 GMM 拓扑结构
        write(*,*) '[EMDAC] 执行高斯混合聚类...'
        call fit_gmm_to_particles(propagated_particles%samples, this%gmm_state, this%em_max_iter, this%em_tol, &
                                 this%current_omega, this%current_W)

        ! 时间更新
        this%current_epoch = et                         
        
        ! 释放临时粒子对象内存，保留 this%propagated_particles 供后续测量更新使用
        call uq_in%deallocate_memory()
        
    end subroutine filter_time_update


    !> ======================================================================
    !> 3. 测量更新 (Measurement Update)
    !> 卡尔曼状态更新 -> 似然计算 -> GMM 权重更新
    !> ======================================================================
    subroutine filter_measurement_update(this, y_meas, noise_R, et, station)
        class(emdac_filter), intent(inout) :: this
        real(DP), dimension(:), intent(in) :: y_meas      ! 真实测量值 y_{k+1}
        real(DP), dimension(:,:), intent(in):: noise_R    ! 测量噪声协方差 R
        real(DP), intent(in)               :: et          ! 历元
        type(observation_station), intent(in):: station   ! 测站信息
        
        class(AlgebraicVector) :: pos_j2000, measurement_da
        integer :: j, n_comp, ny, info
        real(DP), allocatable :: particles_z(:,:) ! 存储每个粒子的预测测量值 (ny, n_particles)
        real(DP), allocatable :: means_z(:,:)   ! 存储每个核的预测测量均值 (ny, n_comp)
        real(DP), allocatable :: P_zz(:,:,:), P_xz(:,:,:), K_gain(:,:,:), P_zz_inv(:,:,:) ! 存储每个核的预测测量协方差、交叉协方差、卡尔曼增益和协方差逆
        ! real(DP), allocatable :: innovation(:), log_likelihood(:)
        real(DP), allocatable :: log_likelihood(:) ! 存储每个核的对数似然概率 (n_comp)
        real(DP), allocatable :: innovation(:, :) ! 存储每个核的测量残差 (ny, n_comp)
        real(DP), allocatable :: det_Pzz(:),  mahalanobis_sq(:) ! 存储每个核的 P_zz 行列式和马氏距离平方 (n_comp)
        real(DP) :: sum_exp
        
        
        n_comp = this%gmm_state%n_components
        ny = size(y_meas)
        dim = this%gmm_state%state_dim ! [修复] 获取状态维度
        
        ! 显式分配所有 allocatable 数组的内存！
        allocate(log_likelihood(n_comp), det_Pzz(n_comp), mahalanobis_sq(n_comp))
        allocate(particles_z(ny, this%n_particles))
        allocate(means_z(ny, n_comp), innovation(ny, n_comp))
        allocate(P_zz(ny, ny, n_comp), P_zz_inv(ny, ny, n_comp))
        allocate(P_xz(dim, ny, n_comp), K_gain(dim, ny, n_comp))
        
        ! 确保此时的滤波器时间与观测时间对齐
        if (abs(this%current_epoch - et) > 1.0e-6_DP) then
            write(*,*) '[警告] 测量更新的历元与滤波器当前历元不匹配！请先调用 time_update。'
        end if

        call dace_initialize(this%da_order, dim)
        call pos_j2000%init(dim)
        do i = 1, 3
            pos_j2000(i) = this%state_mean(i) + da_var(i)
            pos_j2000(i+3) = this%state_mean(i+3) + da_var(i+3)
        end do
    

        call compute_measurement_da(pos_j2000, et, station, 'OPTICAL', measurement_da)

        ! 进一步进行测量更新，更新权重、均值和协方差
        ! 计算每个粒子的预测测量值
        do j = 1, n_particles
            particles_z(:, j) = measurement_da%eval(propagated_particles%samples(:, j) - pos_j2000%elements)
        end do
        
        do i = 1, n_comp
            means_z(:, i) = 0.0_DP
            P_xz(:,:, i) = 0.0_DP
            P_xz(:,:, i) = 0.0_DP
            do j = 1, n_particles
                means_z(:, i) = means_z(:, i) + this%current_omega(i, j) * particles_z(:, j)                                                              
            end do
            means_z(:, i) = means_z(:, i) / this%current_W(i)
            do j = 1, n_particles
                P_zz(:,:, i) = P_zz(:,:, i) + this%current_omega(i, j) * matmul(reshape(particles_z(:, j) - means_z(:, i), (/ny, 1/)), &
                                                                                 reshape(particles_z(:, j) - means_z(:, i), (/1, ny/))) + noise_R
                P_xz(:,:, i) = P_xz(:,:, i) + this%current_omega(i, j) * matmul(reshape(propagated_particles%samples(:, j) - this%state_mean, (/dim, 1/)), &
                                                                                    reshape(particles_z(:, j) - means_z(:, i), (/1, ny/)))   
            end do
            P_zz(:,:, i) = P_zz(:,:, i) / this%current_W(i)
            P_xz(:,:, i) = P_xz(:,:, i) / this%current_W(i)

            ! 计算对数似然概率
            call inverse_and_determinant(P_zz(:,:, i), P_zz_inv(:,:, i), det_Pzz(i), info)
            if (info /= 0) write(*,*) '[警告] P_zz 求逆失败！'

            K_gain(:,:, i) = matmul(P_xz(:,:, i), P_zz_inv(:,:, i))

            this%gmm_state%components(i)%mean = this%gmm_state%components(i)%mean + matmul(K_gain(:,:, i), y_meas - means_z(:, i))
            this%gmm_state%components(i)%cov = this%gmm_state%components(i)%cov - matmul(K_gain(:,:, i), matmul(P_zz(:,:, i), transpose(K_gain(:,:, i))))
            
            innovation(:, i) = y_meas - means_z(:, i)
            mahalanobis_sq(i) = dot_product(innovation(:, i), matmul(P_zz_inv(:,:, i), innovation(:, i)))
            log_likelihood(i) = log(this%gmm_state%components(i)%weight) &
                              - 0.5_DP * ny * log(2.0_DP * PI) &
                              - 0.5_DP * log(det_Pzz(i)) &
                              - 0.5_DP * mahalanobis_sq(i)
        end do

        ! 权重更新：
        ! A. 找到最大似然值防止溢出
        sum_exp = maxval(log_likelihood(1:n_comp))
        
        ! B. 计算平移后的指数和
        det_Pzz(1) = 0.0_DP ! 借用变量作为 LSE 的分母
        do i = 1, n_comp
            det_Pzz(1) = det_Pzz(1) + exp(log_likelihood(i) - sum_exp)
        end do
        
        ! C. 最终权重更新 (注意：log_likelihood 中已含先验权重，不可再乘！)
        do i = 1, n_comp
            this%gmm_state%components(i)%weight = exp(log_likelihood(i) - sum_exp) / det_Pzz(1)
        end do

        ! 更新全局均值
        call this%update_global_mean()
        call this%update_global_cov()

        ! 释放临时数组内存
        deallocate(log_likelihood, particles_z, means_z, P_zz, P_xz, K_gain, P_zz_inv, innovation, det_Pzz, mahalanobis_sq)

    end subroutine filter_measurement_update

    !> ======================================================================
    !> 根据 GMM 的权重和分布，精确采样生成离散粒子云
    !> ======================================================================
    subroutine sample_particles_from_gmm(gmm, n_particles, samples)
        use pod_global, only: DP
        use pod_uq_gmm_state_module, only: uq_gmm_state_type
        ! 引入已有的高效随机数模块
        use pod_random_module, only: generate_multivariate_normal
        
        implicit none
        type(uq_gmm_state_type), intent(in)   :: gmm
        integer, intent(in)                   :: n_particles
        real(DP), dimension(:,:), intent(out) :: samples ! (dim, n_particles)
        
        integer :: n_comp, i
        integer :: remainder, max_weight_idx
        integer, allocatable :: num_per_comp(:)
        real(DP) :: max_weight
        integer :: start_idx, end_idx, p_count
        
        n_comp = gmm%n_components
        
        allocate(num_per_comp(n_comp))
        
        ! ---------------------------------------------------------
        ! 1. 粒子数分配 (解决离散取整误差)
        ! ---------------------------------------------------------
        p_count = 0
        max_weight = -1.0_DP
        max_weight_idx = 1
        
        do i = 1, n_comp
            ! 按权重比例四舍五入分配粒子数
            num_per_comp(i) = nint(gmm%components(i)%weight * real(n_particles, DP))
            p_count = p_count + num_per_comp(i)
            
            ! 记录权重最大的分量，用于吸收分配误差
            if (gmm%components(i)%weight > max_weight) then
                max_weight = gmm%components(i)%weight
                max_weight_idx = i
            end if
        end do
        
        ! 修正总数：将多出或缺少的粒子直接补给权重最大的那个高斯核
        remainder = n_particles - p_count
        num_per_comp(max_weight_idx) = num_per_comp(max_weight_idx) + remainder
        
        ! ---------------------------------------------------------
        ! 2. 利用已有底层接口进行切片批量生成 (调用 pod_random_module)
        ! ---------------------------------------------------------
        start_idx = 1
        
        do i = 1, n_comp
            if (num_per_comp(i) <= 0) cycle
            
            ! 计算当前高斯核应该填充到 samples 矩阵的哪一段列区间
            end_idx = start_idx + num_per_comp(i) - 1
            
            ! 直接复用底层的 OpenMP 多元正态生成器，传入数组切片
            call generate_multivariate_normal(gmm%components(i)%mean, &
                                              gmm%components(i)%cov,  &
                                              samples(:, start_idx:end_idx))
            
            ! 更新下一个核的起始索引
            start_idx = end_idx + 1
        end do
        
        deallocate(num_per_comp)
    end subroutine sample_particles_from_gmm

    subroutine update_global_mean(this)
        class(emdac_filter), intent(inout) :: this
        integer :: i
        this%state_mean = 0.0_DP
        do i = 1, this%gmm_state%n_components
            this%state_mean = this%state_mean + &
                this%gmm_state%components(i)%weight * this%gmm_state%components(i)%mean
        end do
    end subroutine update_global_mean

    subroutine update_global_cov(this)
        class(emdac_filter), intent(inout) :: this
        integer :: i
        real(DP) :: dim
        dim = size(this%state_mean)
        this%state_cov = 0.0_DP
        do i = 1, this%gmm_state%n_components
            this%state_cov = this%state_cov + &
                this%gmm_state%components(i)%weight * &
                (this%gmm_state%components(i)%cov + &
                 matmul(reshape(this%gmm_state%components(i)%mean - this%state_mean, (/dim, 1/)), &
                        reshape(this%gmm_state%components(i)%mean - this%state_mean, (/1, dim/))))
        end do
    end subroutine update_global_cov

end module pod_filter_emdac_module