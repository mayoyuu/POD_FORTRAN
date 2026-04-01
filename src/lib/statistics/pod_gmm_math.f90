module pod_gmm_math_module
    use pod_global, only: DP
    use pod_basicmath_module
    use pod_uq_gmm_state_module, only: uq_gmm_state_type, gaussian_component
    
    implicit none   
    
    private :: kmeans_cluster
    public :: compute_gmm_moments, compute_gmm_correlation
contains

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
            ! 这里的计算完全独立，使用 OpenMP 多线程加速是极其高效的
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


end module pod_gmm_math_module