module pod_gmm_math_module
    use pod_global, only: DP
    use pod_basicmath
    use pod_uq_gmm_state_module, only: uq_gmm_state_type, gaussian_component
    
    implicit none   
    
    private :: kmeans_cluster, EM_step
    public :: fit_gmm_to_particles
contains
    !> 🚀 顶层接口：将粒子集拟合为 GMM，并直接更新 OOP 数据结构
    !> @param[in]    particles  粒子状态矩阵 (dim, n_particles)
    !> @param[inout] gmm        传入时只需分配好 n_components，输出时内部参数全部更新
    !> @param[in]    max_iter   EM 算法最大迭代次数 (建议 50-100)
    !> @param[in]    tol        权重变化的收敛容差 (建议 1e-4)
    subroutine fit_gmm_to_particles(particles, gmm, max_iter, tol)
        real(DP), dimension(:,:), intent(in) :: particles
        type(uq_gmm_state_type), intent(inout) :: gmm
        integer, intent(in) :: max_iter
        real(DP), intent(in) :: tol
        
        integer :: n_clusters, dim, n_particles, i, j, k,c,r, iter
        real(DP), allocatable :: means(:,:), covariances(:,:,:), weights(:), old_weights(:), diff(:)
        integer, allocatable :: assignments(:), counts(:)
        real(DP) :: weight_diff
        
        n_clusters = gmm%n_components
        dim = size(particles, 1)
        n_particles = size(particles, 2)
        
        ! 1. 分配连续内存用于高性能数学计算 (避免直接在 OOP 结构里做繁重的跨步内存访问)
        allocate(means(dim, n_clusters))
        allocate(covariances(dim, dim, n_clusters))
        allocate(weights(n_clusters))
        allocate(old_weights(n_clusters))
        
        
        ! ... 在之前的 allocate 块中加入：
        allocate(assignments(n_particles), counts(n_clusters))

        ! 2. K-means 初始化 
        ! 必须确保你的 kmeans_cluster 将每个粒子的归属输出到 assignments 中
        call kmeans_cluster(particles, n_clusters, means, assignments) 
        
        ! ---------------------------------------------------------
        ! ✨ 严格遵循论文逻辑：根据 K-means 硬聚类结果初始化 EM 参数
        ! ---------------------------------------------------------
        counts = 0
        covariances = 0.0_DP
        
        ! 统计每个簇的粒子数量，并累加协方差的外积
        do j = 1, n_particles
            k = assignments(j)
            counts(k) = counts(k) + 1
            
            diff = particles(:, j) - means(:, k)
            ! 展开外积加速
            do c = 1, dim
                do r = 1, dim
                    covariances(r, c, k) = covariances(r, c, k) + diff(r) * diff(c)
                end do
            end do
        end do
        
        ! 计算初始权重与样本协方差
        do i = 1, n_clusters
            if (counts(i) > 1) then
                ! 权重 = 该簇粒子数 / 总粒子数
                weights(i) = real(counts(i), DP) / real(n_particles, DP)
                ! 样本协方差 = 累加外积 / (N - 1) 
                covariances(:,:,i) = covariances(:,:,i) / real(counts(i) - 1, DP)
            else
                ! 极端退化保护（如果某个簇极其不幸没分到粒子，给它一个极小的初始值避免除零）
                weights(i) = 1.0e-8_DP
                do j = 1, dim
                    covariances(j,j,i) = 1.0e-6_DP
                end do
            end if
        end do
        
        ! 3. EM 算法迭代循环
        do iter = 1, max_iter
            old_weights = weights
            
            ! 调用核心单步 EM (也就是你刚才写的那个完美优化的版本)
            call EM_step(particles, n_clusters, means, covariances, weights)
            
            ! 收敛性检查：比较前后两次权重的变化 (论文 3.3 节末尾的推荐准则)
            weight_diff = maxval(abs(weights - old_weights))
            if (weight_diff < tol) exit
        end do
        
        ! 4. 将计算好的连续内存数据，塞回 OOP 的 gmm 对象中
        do i = 1, n_clusters
            ! 这里调用我们在 pod_uq_gmm_state_module 里写的 init 绑定过程
            call gmm%components(i)%init(weights(i), means(:, i), covariances(:, :, i))
        end do
        
        deallocate(means, covariances, weights, old_weights)
    end subroutine fit_gmm_to_particles
    !> kmeans
    !> 核心的 K-means 聚类算法 (Lloyd's Algorithm)
    !> @param[in]  particles   所有粒子的状态矩阵 (dim, n_particles)
    !> @param[in]  n_clusters  聚类中心数量 (K)
    !> @param[out] means       输出的聚类中心矩阵 (dim, n_clusters)
    !> @param[out] assignments 记录每个粒子属于哪个中心 (1, n_particles)
    subroutine kmeans_cluster(particles, n_clusters, means, assignments)
        real(DP), dimension(:,:), intent(in)  :: particles
        integer, intent(in)                   :: n_clusters
        real(DP), dimension(:,:), intent(out) :: means
        integer, dimension(:), intent(out)    :: assignments
        
        integer :: dim, n_particles, iter, i, j, best_k
        real(DP) :: min_dist, dist
        integer, dimension(n_clusters) :: counts
        logical :: converged
        
        dim = size(particles, 1)
        n_particles = size(particles, 2)
        
        ! 1. 简单初始化：随机挑选前 n_clusters 个粒子作为初始中心
        ! (你也可以调用 pod_random 模块来随机抽取索引，这是最基础的写法)
        means = particles(:, 1:n_clusters)
        
        ! 最大迭代次数保护
        do iter = 1, 100
            converged = .true.
            
            ! ---------------------------------------------------------
            ! 步骤 A: 期望步 (E-step) - 为每个粒子寻找最近的中心点
            ! ---------------------------------------------------------
            !$omp parallel do default(none) &
            !$omp shared(n_particles, n_clusters, dim, particles, means, assignments, converged) &
            !$omp private(i, j, best_k, min_dist, dist)
            do i = 1, n_particles
                min_dist = huge(1.0_DP)
                best_k = 1
                
                ! 遍历所有中心，计算欧氏距离的平方 (省去开平方运算以加速)
                do j = 1, n_clusters
                    dist = sum((particles(:, i) - means(:, j))**2)
                    if (dist < min_dist) then
                        min_dist = dist
                        best_k = j
                    end if
                end do
                
                ! 检查是否收敛 (如果有任何一个粒子改变了阵营，说明还没收敛)
                if (assignments(i) /= best_k) converged = .false.
                assignments(i) = best_k
            end do
            !$omp end parallel do
            
            if (converged) exit
            
            ! ---------------------------------------------------------
            ! 步骤 B: 最大化步 (M-step) - 重新计算每个聚类的均值
            ! ---------------------------------------------------------
            means = 0.0_DP
            counts = 0
            
            ! 累加同组的粒子坐标
            do i = 1, n_particles
                best_k = assignments(i)
                means(:, best_k) = means(:, best_k) + particles(:, i)
                counts(best_k) = counts(best_k) + 1
            end do
            
            ! 计算平均值
            do j = 1, n_clusters
                if (counts(j) > 0) then
                    means(:, j) = means(:, j) / real(counts(j), DP)
                else
                    ! 处理空簇异常：如果某个簇空了，随机拿一个粒子顶替 (K-means 常见边缘情况)
                    means(:, j) = particles(:, j)
                end if
            end do
            
        end do
        
    end subroutine kmeans_cluster

    subroutine EM_step(particles, n_clusters, means, covariances, weights)
        real(DP), dimension(:,:), intent(in)      :: particles
        integer, intent(in)                       :: n_clusters
        real(DP), dimension(:,:), intent(inout)   :: means
        real(DP), dimension(:,:,:), intent(inout) :: covariances
        real(DP), dimension(:), intent(inout)     :: weights

        integer :: n_particles, dim, i, j, h, r, c
        integer :: info
        real(DP) :: PI = 3.14159265358979323846_DP
        real(DP) :: det_cov, mahalanobis_sq, sum_exp
        
        ! 利用自动数组声明维度，彻底干掉没用的冗余中间变量
        real(DP), dimension(size(particles, 1))                       :: diff
        real(DP), dimension(size(particles, 1), size(particles, 1))   :: cov_inv
        
        ! 只需要这三个核心工作矩阵
        real(DP), dimension(size(particles, 2), n_clusters)           :: log_P
        real(DP), dimension(n_clusters, size(particles, 2))           :: omega
        real(DP), dimension(n_clusters)                               :: responsibility

        dim = size(particles, 1)
        n_particles = size(particles, 2)

        !===================================================================
        ! 步骤 1: E-step (计算对数概率密度)
        !===================================================================
        do i = 1, n_clusters
            ! ✨ 优化1：每个簇只求一次逆矩阵和行列式！
            call inverse_and_determinant(covariances(:,:,i), cov_inv, det_cov, info)
            
            if (info > 0 .or. det_cov <= 0.0_DP) then
                write(*,*) '[ERROR] 第', i, '个簇的协方差矩阵非正定或奇异！'
                return
            end if
            
            do j = 1, n_particles
                diff = particles(:, j) - means(:, i)
                mahalanobis_sq = dot_product(diff, matmul(cov_inv, diff))
                
                ! ✨ 优化2：直接将 weight 融入对数概率中，彻底抛弃线性 rho
                ! log( P(x_j | \mu_i, \Sigma_i) * weight_i )
                log_P(j, i) = log(weights(i)) &
                            - 0.5_DP * dim * log(2.0_DP * PI) &
                            - 0.5_DP * log(det_cov) &
                            - 0.5_DP * mahalanobis_sq
            end do
        end do

        !===================================================================
        ! 步骤 2: E-step (应用论文 Eq.19 计算责任度 omega)
        !===================================================================
        do j = 1, n_particles
            do i = 1, n_clusters
                sum_exp = 0.0_DP
                do h = 1, n_clusters
                    ! 完全在指数层级相减，完美避开 absolute 0 的下溢问题
                    sum_exp = sum_exp + exp(log_P(j, h) - log_P(j, i))
                end do
                ! 此时的 omega(i,j) 就是归一化后的确切概率，且绝对安全
                omega(i, j) = 1.0_DP / sum_exp
            end do
        end do

        !===================================================================
        ! 步骤 3: M-step (更新 GMM 参数)
        !===================================================================
        do i = 1, n_clusters
            ! 该簇的总责任度 (论文中的 \mathcal{W}_i)
            responsibility(i) = sum(omega(i, :))
            
            ! 1. 更新均值 
            means(:, i) = sum(particles * reshape(omega(i, :), [1, n_particles]), dim=2) / responsibility(i)
            
            ! 2. 更新协方差矩阵
            covariances(:, :, i) = 0.0_DP
            do j = 1, n_particles
                diff = particles(:, j) - means(:, i)
                
                ! ✨ 优化3：展开外积计算。避免使用 matmul 和 reshape，
                ! 这样可以阻止 Fortran 在循环内部分配临时内存块，极大提升速度
                do c = 1, dim
                    do r = 1, dim
                        covariances(r, c, i) = covariances(r, c, i) + omega(i, j) * diff(r) * diff(c)
                    end do
                end do
            end do
            covariances(:, :, i) = covariances(:, :, i) / responsibility(i)
            
            ! 3. 更新权重
            ! 按照 EM 算法理论，权重的分母就是所有簇的总责任度，等于粒子总数 n_particles
            weights(i) = responsibility(i) / real(n_particles, DP)
        end do

    end subroutine EM_step


end module pod_gmm_math_module