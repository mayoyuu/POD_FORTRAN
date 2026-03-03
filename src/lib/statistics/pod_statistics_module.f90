module pod_statistics_module
    use pod_global, only: DP
    
    implicit none
    
contains

    subroutine least_squares_fit(A, b, x, residual)
        real(DP), dimension(:,:), intent(in) :: A
        real(DP), dimension(:), intent(in) :: b
        real(DP), dimension(:), intent(out) :: x
        real(DP), intent(out) :: residual
        
        real(DP), dimension(:,:), allocatable :: ATA, ATA_inv
        real(DP), dimension(:), allocatable :: ATb
        integer :: m, n
        
        m = size(A, 1)
        n = size(A, 2)
        
        ! 分配内存
        allocate(ATA(n, n), ATA_inv(n, n), ATb(n))
        
        ! 计算 A^T * A
        ATA = matmul(transpose(A), A)
        
        ! 计算 A^T * b
        ATb = matmul(transpose(A), b)
        
        ! 计算 (A^T * A)^(-1)
        call matrix_inverse(ATA, ATA_inv)
        
        ! 计算解: x = (A^T * A)^(-1) * A^T * b
        x = matmul(ATA_inv, ATb)
        
        ! 计算残差
        residual = sqrt(sum((matmul(A, x) - b)**2))
        
        deallocate(ATA, ATA_inv, ATb)
    end subroutine least_squares_fit
    
    subroutine weighted_least_squares_fit(A, b, W, x, residual)
        real(DP), dimension(:,:), intent(in) :: A
        real(DP), dimension(:), intent(in) :: b
        real(DP), dimension(:,:), intent(in) :: W
        real(DP), dimension(:), intent(out) :: x
        real(DP), intent(out) :: residual
        
        real(DP), dimension(:,:), allocatable :: ATWA, ATWA_inv
        real(DP), dimension(:), allocatable :: ATWb
        integer :: m, n
        
        m = size(A, 1)
        n = size(A, 2)
        
        ! 分配内存
        allocate(ATWA(n, n), ATWA_inv(n, n), ATWb(n))
        
        ! 计算 A^T * W * A
        ATWA = matmul(matmul(transpose(A), W), A)
        
        ! 计算 A^T * W * b
        ATWb = matmul(matmul(transpose(A), W), b)
        
        ! 计算 (A^T * W * A)^(-1)
        call matrix_inverse(ATWA, ATWA_inv)
        
        ! 计算解: x = (A^T * W * A)^(-1) * A^T * W * b
        x = matmul(ATWA_inv, ATWb)
        
        ! 计算加权残差
        residual = sqrt(sum((matmul(A, x) - b)**2))
        
        deallocate(ATWA, ATWA_inv, ATWb)
    end subroutine weighted_least_squares_fit
    
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
    
    subroutine compute_covariance_matrix(A, W, covariance)
        real(DP), dimension(:,:), intent(in) :: A
        real(DP), dimension(:,:), intent(in) :: W
        real(DP), dimension(:,:), intent(out) :: covariance
        
        real(DP), dimension(:,:), allocatable :: ATWA, ATWA_inv
        integer :: n
        
        n = size(A, 2)
        allocate(ATWA(n, n), ATWA_inv(n, n))
        
        ! 计算 A^T * W * A
        ATWA = matmul(matmul(transpose(A), W), A)
        
        ! 计算 (A^T * W * A)^(-1)
        call matrix_inverse(ATWA, ATWA_inv)
        
        ! 协方差矩阵就是 (A^T * W * A)^(-1)
        covariance = ATWA_inv
        
        deallocate(ATWA, ATWA_inv)
    end subroutine compute_covariance_matrix
    
    subroutine compute_chi_square(measurements, predicted, weights, chi_square, degrees_of_freedom)
        real(DP), dimension(:), intent(in) :: measurements, predicted, weights
        real(DP), intent(out) :: chi_square
        integer, intent(out) :: degrees_of_freedom
        
        real(DP), dimension(:), allocatable :: residuals
        
        allocate(residuals(size(measurements)))
        
        ! 计算残差
        residuals = measurements - predicted
        
        ! 计算加权卡方统计量
        chi_square = sum(weights * residuals**2)
        
        ! 自由度
        degrees_of_freedom = size(measurements) - 6  ! 假设有6个轨道参数
    end subroutine compute_chi_square
    
    subroutine compute_confidence_ellipse(covariance, confidence_level, semi_major, semi_minor, orientation)
        real(DP), dimension(6,6), intent(in) :: covariance
        real(DP), intent(in) :: confidence_level
        real(DP), intent(out) :: semi_major, semi_minor, orientation
        
        real(DP), dimension(6,6) :: position_covariance
        real(DP), dimension(3,3) :: pos_cov_3x3
        real(DP), dimension(3) :: eigenvalues
        real(DP), dimension(3,3) :: eigenvectors
        real(DP) :: chi_square_value
        
        ! 提取位置协方差矩阵 (3x3)
        pos_cov_3x3 = covariance(1:3, 1:3)
        
        ! 计算特征值和特征向量
        call compute_eigenvalues_eigenvectors(pos_cov_3x3, eigenvalues, eigenvectors)
        
        ! 计算置信度对应的卡方值
        call chi_square_inverse(confidence_level, 3, chi_square_value)
        
        ! 计算椭球半轴
        semi_major = sqrt(chi_square_value * eigenvalues(1))
        semi_minor = sqrt(chi_square_value * eigenvalues(2))
        
        ! 计算方向角
        orientation = atan2(eigenvectors(2,1), eigenvectors(1,1)) * 180.0_DP / (4.0_DP * atan(1.0_DP))
    end subroutine compute_confidence_ellipse
    
    subroutine compute_eigenvalues_eigenvectors(A, eigenvalues, eigenvectors)
        real(DP), dimension(:,:), intent(in) :: A
        real(DP), dimension(:), intent(out) :: eigenvalues
        real(DP), dimension(:,:), intent(out) :: eigenvectors
        
        integer :: n, info, i
        real(DP), allocatable, dimension(:) :: work
        real(DP), allocatable, dimension(:,:) :: A_copy
        
        n = size(A, 1)
        allocate(A_copy(n, n), work(3*n))
        
        A_copy = A
        eigenvectors = A_copy
        
        ! 简化实现，不使用LAPACK
        ! 这里应该实现特征值计算算法
        eigenvalues = 1.0_DP  ! 临时简化
        eigenvectors = 0.0_DP
        do i = 1, n
            eigenvectors(i, i) = 1.0_DP
        end do
        
        deallocate(A_copy, work)
    end subroutine compute_eigenvalues_eigenvectors
    
    subroutine chi_square_inverse(p, df, chi_square_value)
        real(DP), intent(in) :: p
        integer, intent(in) :: df
        real(DP), intent(out) :: chi_square_value
        
        ! 简化实现，实际应用中需要使用更精确的数值方法
        ! 这里使用近似公式
        real(DP) :: z, mu, sigma
        
        ! 使用正态分布近似
        call normal_inverse(p, z)
        
        mu = df
        sigma = sqrt(2.0_DP * df)
        
        chi_square_value = mu + sigma * z
    end subroutine chi_square_inverse
    
    subroutine normal_inverse(p, z)
        real(DP), intent(in) :: p
        real(DP), intent(out) :: z
        
        ! 简化实现，使用近似公式
        real(DP) :: t
        
        if (p < 0.5_DP) then
            t = sqrt(-2.0_DP * log(p))
        else
            t = sqrt(-2.0_DP * log(1.0_DP - p))
        end if
        
        z = t - (2.30753_DP + 0.27061_DP * t) / (1.0_DP + 0.99229_DP * t + 0.04481_DP * t**2)
        
        if (p < 0.5_DP) z = -z
    end subroutine normal_inverse
    
    subroutine compute_residual_statistics(residuals, mean, std, rms)
        real(DP), dimension(:), intent(in) :: residuals
        real(DP), intent(out) :: mean, std, rms
        integer :: n
        
        n = size(residuals)
        
        ! 计算均值
        mean = sum(residuals) / n
        
        ! 计算标准差
        std = sqrt(sum((residuals - mean)**2) / (n - 1))
        
        ! 计算均方根
        rms = sqrt(sum(residuals**2) / n)
    end subroutine compute_residual_statistics
    
    subroutine compute_correlation_matrix(covariance, correlation)
        real(DP), dimension(:,:), intent(in) :: covariance
        real(DP), dimension(:,:), intent(out) :: correlation
        
        integer :: n, i, j
        real(DP), dimension(:), allocatable :: std_dev
        
        n = size(covariance, 1)
        allocate(std_dev(n))
        
        ! 计算标准差
        do i = 1, n
            std_dev(i) = sqrt(covariance(i, i))
        end do
        
        ! 计算相关系数矩阵
        do i = 1, n
            do j = 1, n
                if (std_dev(i) > 0.0_DP .and. std_dev(j) > 0.0_DP) then
                    correlation(i, j) = covariance(i, j) / (std_dev(i) * std_dev(j))
                else
                    correlation(i, j) = 0.0_DP
                end if
            end do
        end do
        
        deallocate(std_dev)
    end subroutine compute_correlation_matrix

end module pod_statistics_module
