!-----------------------------------------------------------------------------------------
! 模块: pod_uq_transform_da_module
! 功能: 纯数学工具 - 基于微分代数 (DA) 的测量多项式极速求值与统计矩重建
! 对应: EMDAC-N 测量更新步的 DA-MC 快速映射
!-----------------------------------------------------------------------------------------
module pod_uq_transform_da_module
    use pod_global, only: DP
    use pod_dace_classes
    implicit none
    private

    ! 对外暴露的三个纯数学/代数工具接口 (与 UT 模块完美平行)
    public :: init_da_expansion
    public :: evaluate_da_polynomial
    public :: reconstruct_da_moments

contains

    !> ======================================================================
    !> 步骤 1: 初始化 DA 展开中心
    !> 对应 UT 的 generate_sigma_points
    !> ======================================================================
    subroutine init_da_expansion(mean_x, da_order, state_da_0)
        real(DP), dimension(:), intent(in)   :: mean_x     ! 当前高斯核均值 \hat{x}_i
        integer, intent(in)                  :: da_order   ! DA 展开阶数
        type(AlgebraicVector), intent(inout) :: state_da_0 ! 输出的 DA 初始状态

        integer :: nx, i
        nx = size(mean_x)

        ! 初始化 DACE 核心并分配变量
        call dace_initialize(da_order, nx)
        call state_da_0%init(nx)
        
        ! 注入中心值和独立偏差 dA (使用你暴露的 da_var 函数)
        do i = 1, nx
            state_da_0%elements(i) = mean_x(i) + da_var(i)
        end do
    end subroutine init_da_expansion


    !> ======================================================================
    !> 步骤 2: 将测量多项式对所有粒子进行极速求值
    !> 获取所有粒子的预测测量值 q_{j}
    !> ======================================================================
    subroutine evaluate_da_polynomial(meas_da_f, mean_x, particles_x, particles_y)
        type(AlgebraicVector), intent(in) :: meas_da_f     ! 用户物理函数输出的高阶测量多项式
        real(DP), dimension(:), intent(in):: mean_x        ! 展开中心 \hat{x}_i
        real(DP), dimension(:,:), intent(in) :: particles_x ! 输入粒子集 (nx, n_particles)
        real(DP), dimension(:,:), allocatable, intent(out) :: particles_y ! 输出求值结果 q_j (ny, n_particles)

        type(CompiledDA) :: compiled_meas
        integer :: nx, ny, n_particles, j
        real(DP), allocatable :: eval_inputs(:), eval_results(:)

        nx = size(particles_x, 1)
        n_particles = size(particles_x, 2)
        
        ! 🚀 修正点：调用 AlgebraicVector 的属性 %size，而不是函数
        ny = meas_da_f%size  

        allocate(particles_y(ny, n_particles))
        
        ! 编译多项式以获得极致性能
        compiled_meas = meas_da_f%compile()

        !$omp parallel do default(none) private(j, eval_inputs, eval_results) &
        !$omp shared(n_particles, nx, ny, particles_x, mean_x, compiled_meas, particles_y)
        do j = 1, n_particles
            allocate(eval_inputs(nx))
            ! 偏差 = 当前粒子坐标 - 展开中心
            eval_inputs = particles_x(:, j) - mean_x
            
            ! eval_results 会在调用 CompiledDA%eval 时被自动 allocate
            eval_results = compiled_meas%eval(eval_inputs)
            
            particles_y(:, j) = eval_results
            
            ! 循环末尾显式释放内存，防止 OMP 线程内的内存泄漏
            deallocate(eval_inputs, eval_results)
        end do
        !$omp end parallel do
        
        call compiled_meas%destroy()
    end subroutine evaluate_da_polynomial


    !> ======================================================================
    !> 步骤 3: 利用 EM 软权重重建交叉协方差与自协方差
    !> 对应 UT 的 reconstruct_ut_moments
    !> ======================================================================
    subroutine reconstruct_da_moments(particles_x, particles_y, em_weights, mean_x, &
                                      mean_z, P_zz, P_xz, is_angle)
        real(DP), dimension(:,:), intent(in) :: particles_x  ! \tilde{\delta}_j + \hat{x}_i
        real(DP), dimension(:,:), intent(in) :: particles_y  ! q_j
        real(DP), dimension(:), intent(in)   :: em_weights   ! EM 输出的归一化软权重 w_{i,j}
        real(DP), dimension(:), intent(in)   :: mean_x       ! \hat{x}_i
        
        real(DP), dimension(:), allocatable, intent(out)   :: mean_z
        real(DP), dimension(:,:), allocatable, intent(out) :: P_zz
        real(DP), dimension(:,:), allocatable, intent(out) :: P_xz
        logical, dimension(:), optional, intent(in)        :: is_angle

        integer :: nx, ny, n_particles, j, r, c
        real(DP), allocatable :: diff_x(:), diff_z(:)
        real(DP) :: sum_sin, sum_cos
        real(DP) :: PI = 3.14159265358979323846_DP

        nx = size(particles_x, 1)
        ny = size(particles_y, 1)
        n_particles = size(particles_x, 2)

        allocate(mean_z(ny), P_zz(ny, ny), P_xz(nx, ny))
        allocate(diff_x(nx), diff_z(ny))

        ! 1. 计算加权均值 \hat{z}_i
        mean_z = 0.0_DP
        do j = 1, n_particles
            mean_z = mean_z + em_weights(j) * particles_y(:, j)
        end do

        ! 角度均值跨界平滑修正
        if (present(is_angle)) then
            do r = 1, ny
                if (is_angle(r)) then
                    sum_sin = 0.0_DP; sum_cos = 0.0_DP
                    do j = 1, n_particles
                        sum_sin = sum_sin + em_weights(j) * sin(particles_y(r, j))
                        sum_cos = sum_cos + em_weights(j) * cos(particles_y(r, j))
                    end do
                    mean_z(r) = atan2(sum_sin, sum_cos)
                end if
            end do
        end if

        ! 2. 计算自协方差 P_zz 和 交叉协方差 P_xz
        P_zz = 0.0_DP
        P_xz = 0.0_DP

        do j = 1, n_particles
            diff_x = particles_x(:, j) - mean_x
            diff_z = particles_y(:, j) - mean_z

            ! 角度偏差边界截断 [-PI, PI]
            if (present(is_angle)) then
                do r = 1, ny
                    if (is_angle(r)) then
                        if (diff_z(r) > PI)  diff_z(r) = diff_z(r) - 2.0_DP * PI
                        if (diff_z(r) < -PI) diff_z(r) = diff_z(r) + 2.0_DP * PI
                    end if
                end do
            end if

            ! 展开矩阵外积计算，极大提升 Cache 命中率
            do c = 1, ny
                do r = 1, ny
                    P_zz(r, c) = P_zz(r, c) + em_weights(j) * diff_z(r) * diff_z(c)
                end do
                do r = 1, nx
                    P_xz(r, c) = P_xz(r, c) + em_weights(j) * diff_x(r) * diff_z(c)
                end do
            end do
        end do
        
        deallocate(diff_x, diff_z)
    end subroutine reconstruct_da_moments

end module pod_uq_transform_da_module