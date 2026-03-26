module pod_random_module
    use pod_global, only: DP
    implicit none
    private
    
    ! 暴露给外部调用的接口
    public :: init_random_seed
    public :: randn
    public :: generate_multivariate_normal

    real(DP), parameter :: PI = 3.14159265358979323846_DP

    ! ========================================================
    ! 显式声明 LAPACK 外部接口 (解决 implicit-interface 报错)
    ! ========================================================
    interface
        subroutine dpotrf(uplo, n, a, lda, info)
            import :: DP
            implicit none
            character(len=1), intent(in) :: uplo
            integer, intent(in)          :: n
            real(DP), intent(inout)      :: a(lda, *)
            integer, intent(in)          :: lda
            integer, intent(out)         :: info
        end subroutine dpotrf
    end interface

contains

    !> 初始化随机数种子 (支持 OpenMP 线程安全)
    subroutine init_random_seed(repeatable)
        logical, intent(in), optional :: repeatable
        logical :: is_repeatable
        
        if (present(repeatable)) then
            is_repeatable = repeatable
        else
            is_repeatable = .false. ! 默认每次运行生成不同的随机数
        end if
        
        ! Fortran 2018 标准：自动为所有 OpenMP 线程分配不同的种子
        ! arg1: repeatable (是否每次运行结果一致)
        ! arg2: image_distinct (多映像/多线程是否独立)
        call random_init(is_repeatable, .true.)
    end subroutine init_random_seed

    !> 生成单个标准正态分布随机数 N(0, 1) - Box-Muller 方法
    real(DP) function randn() result(z)
        real(DP) :: u1, u2
        
        ! 获取两个 (0, 1) 之间的均匀分布随机数
        call random_number(u1)
        call random_number(u2)
        
        ! 防止 u1 过小导致 log 报错
        if (u1 < epsilon(1.0_DP)) u1 = epsilon(1.0_DP)
        
        z = sqrt(-2.0_DP * log(u1)) * cos(2.0_DP * PI * u2)
    end function randn

    !> 基于给定的均值和协方差矩阵，生成 N 个多元正态分布粒子
    subroutine generate_multivariate_normal(mean, cov, samples)
        real(DP), intent(in)    :: mean(:)        ! 均值向量，大小 (dim)
        real(DP), intent(in)    :: cov(:,:)       ! 协方差矩阵，大小 (dim, dim)
        real(DP), intent(inout) :: samples(:,:)   ! 输出样本，大小 (dim, n_particles)
        
        integer :: dim, n_particles, i, j, info
        real(DP), allocatable :: L(:,:)
        real(DP), allocatable :: Z(:,:)
        
        dim = size(mean)
        n_particles = size(samples, 2)
        
        ! 1. 复制协方差矩阵 (因为 dpotrf 会覆盖输入矩阵)
        allocate(L(dim, dim))
        L = cov
        
        ! 2. 调用 LAPACK 进行 Cholesky 分解: C = L * L^T
        ! 'L' 表示计算并返回下三角矩阵
        call dpotrf('L', dim, L, dim, info)
        
        if (info /= 0) then
            write(*,*) '[ERROR] Cholesky 分解失败！协方差矩阵可能非正定。错误代码: ', info
            ! 可以在这里加入特征值分解(SVD)作为后备方案，以处理奇异矩阵
            return
        end if
        
        ! 清理上三角部分 (dpotrf 不会自动清零上三角)
        do i = 1, dim
            do j = i + 1, dim
                L(i, j) = 0.0_DP
            end do
        end do
        
        ! 3. 生成标准正态分布随机矩阵 Z (dim x n_particles)
        allocate(Z(dim, n_particles))
        
        !$omp parallel do default(none) shared(dim, n_particles, Z) private(i, j)
        do j = 1, n_particles
            do i = 1, dim
                Z(i, j) = randn()
            end do
        end do
        !$omp end parallel do
        
        ! 4. 线性映射 X = mean + L * Z
        ! 使用 BLAS 的 dgemm 或 Fortran 内置的 matmul (内置 matmul 通常已优化)
        samples = matmul(L, Z)
        
        ! 将均值加到每一个粒子上
        !$omp parallel do default(none) shared(dim, n_particles, samples, mean) private(i, j)
        do j = 1, n_particles
            do i = 1, dim
                samples(i, j) = samples(i, j) + mean(i)
            end do
        end do
        !$omp end parallel do
        
        deallocate(L)
        deallocate(Z)
    end subroutine generate_multivariate_normal

end module pod_random_module