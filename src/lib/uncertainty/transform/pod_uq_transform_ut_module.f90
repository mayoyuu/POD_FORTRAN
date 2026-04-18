!-----------------------------------------------------------------------------------------
! 模块: pod_uq_transform_ut_module
! 功能: 纯数学工具 - 无迹变换 (Unscented Transform)
!       负责生成 Sigma 点，以及根据变换后的 Sigma 点重建统计矩
!-----------------------------------------------------------------------------------------
module pod_uq_transform_ut_module
    use pod_global, only: DP
    implicit none
    private

    ! 对外暴露的两个核心纯数学接口
    public :: generate_sigma_points
    public :: reconstruct_ut_moments

    ! 声明 LAPACK Cholesky 分解接口
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

    !> ======================================================================
    !> 步骤 1: 根据均值和协方差生成 2n+1 个 Sigma 点及对应权重
    !> 完美映射你 C++ 中的 pointut0 和 W 逻辑
    !> ======================================================================
    subroutine generate_sigma_points(mean_x, cov_x, kappa, sigmas_x, weights)
        real(DP), dimension(:), intent(in) :: mean_x            ! 输入均值 (n)
        real(DP), dimension(:,:), intent(in) :: cov_x           ! 输入协方差 (n, n)
        real(DP), intent(in) :: kappa                           ! UT参数 (你的代码中为 0.5)
        real(DP), dimension(:,:), allocatable, intent(out) :: sigmas_x  ! 输出Sigma点 (n, 2n+1)
        real(DP), dimension(:), allocatable, intent(out) :: weights     ! 输出权重 (2n+1)

        integer :: n, n_sigma, i, info
        real(DP) :: scale
        real(DP), allocatable :: L(:,:)

        n = size(mean_x)
        n_sigma = 2 * n + 1
        scale = sqrt(real(n, DP) + kappa)

        allocate(sigmas_x(n, n_sigma))
        allocate(weights(n_sigma))
        allocate(L(n, n))

        ! 1. 计算权重 W
        weights(1) = kappa / (real(n, DP) + kappa)
        do i = 2, n_sigma
            weights(i) = 1.0_DP / (2.0_DP * (real(n, DP) + kappa))
        end do

        ! 2. Cholesky 分解获取 cov 的平方根 L
        L = cov_x
        call dpotrf('L', n, L, n, info)
        if (info /= 0) then
            write(*,*) '[ERROR] UT Transform: Cholesky decomposition failed!'
            return
        end if
        
        ! 清理上三角部分
        do i = 1, n
            L(1:i-1, i) = 0.0_DP
        end do

        ! 3. 生成 2n+1 个 Sigma 点
        sigmas_x(:, 1) = mean_x
        do i = 1, n
            sigmas_x(:, i + 1)     = mean_x + scale * L(:, i)
            sigmas_x(:, i + 1 + n) = mean_x - scale * L(:, i)
        end do

        deallocate(L)
    end subroutine generate_sigma_points


    !> ======================================================================
    !> 步骤 2: 重建统计矩 (均值, 协方差, 交叉协方差)
    !> 对应你 C++ 中 mean_end 和 cov_out 的计算
    !> ======================================================================
    subroutine reconstruct_ut_moments(sigmas_x, sigmas_y, weights, mean_y, P_yy, P_xy, is_angle)
        real(DP), dimension(:,:), intent(in) :: sigmas_x  ! 原始点 (nx, 2n+1)
        real(DP), dimension(:,:), intent(in) :: sigmas_y  ! 经过非线性函数映射后的点 (ny, 2n+1)
        real(DP), dimension(:), intent(in)   :: weights   ! 权重 (2n+1)
        real(DP), dimension(:), allocatable, intent(out)   :: mean_y  ! 输出均值 (ny)
        real(DP), dimension(:,:), allocatable, intent(out) :: P_yy    ! 输出自协方差 (ny, ny)
        real(DP), dimension(:,:), allocatable, intent(out) :: P_xy    ! 输出交叉协方差 (nx, ny)
        ! 可选参数：标记输出变量中哪些维度是角度（例如 RA, DEC），用于修复 0-360 度跨界问题
        logical, dimension(:), optional, intent(in)        :: is_angle ! (ny) 

        integer :: nx, ny, n_sigma, i, j, r, c
        real(DP), allocatable :: diff_x(:), diff_y(:)
        real(DP) :: sum_sin, sum_cos
        real(DP) :: PI = 3.14159265358979323846_DP

        nx = size(sigmas_x, 1)
        ny = size(sigmas_y, 1)
        n_sigma = size(weights)

        allocate(mean_y(ny), P_yy(ny, ny), P_xy(nx, ny))
        allocate(diff_x(nx), diff_y(ny))

        ! 1. 计算均值
        mean_y = 0.0_DP
        do i = 1, n_sigma
            mean_y = mean_y + weights(i) * sigmas_y(:, i)
        end do

        ! ✨ 角度均值修正: 防止 359 度和 1 度平均出 180 度的致命错误
        if (present(is_angle)) then
            do j = 1, ny
                if (is_angle(j)) then
                    sum_sin = 0.0_DP; sum_cos = 0.0_DP
                    do i = 1, n_sigma
                        sum_sin = sum_sin + weights(i) * sin(sigmas_y(j, i))
                        sum_cos = sum_cos + weights(i) * cos(sigmas_y(j, i))
                    end do
                    mean_y(j) = atan2(sum_sin, sum_cos)
                end if
            end do
        end if

        ! 2. 计算自协方差 P_yy 与交叉协方差 P_xy
        P_yy = 0.0_DP
        P_xy = 0.0_DP

        do i = 1, n_sigma
            diff_x = sigmas_x(:, i) - sigmas_x(:, 1) ! sigmas_x(:, 1) 永远是未扰动的原始均值
            diff_y = sigmas_y(:, i) - mean_y

            ! 角度偏差修正 (强制映射回 -PI 到 PI 之间)
            if (present(is_angle)) then
                do j = 1, ny
                    if (is_angle(j)) then
                        if (diff_y(j) > PI)  diff_y(j) = diff_y(j) - 2.0_DP * PI
                        if (diff_y(j) < -PI) diff_y(j) = diff_y(j) + 2.0_DP * PI
                    end if
                end do
            end if

            ! 展开矩阵外积计算，对缓存命中率极高
            do c = 1, ny
                do r = 1, ny
                    P_yy(r, c) = P_yy(r, c) + weights(i) * diff_y(r) * diff_y(c)
                end do
                ! 卡尔曼更新的灵魂: X与Y的交叉协方差
                do r = 1, nx
                    P_xy(r, c) = P_xy(r, c) + weights(i) * diff_x(r) * diff_y(c)
                end do
            end do
        end do

        deallocate(diff_x, diff_y)
    end subroutine reconstruct_ut_moments

end module pod_uq_transform_ut_module