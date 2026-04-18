!-----------------------------------------------------------------------------------------
! 模块: pod_uq_gmm_state_module
! 功能: 定义高斯混合模型(GMM)的数据结构与基础操作
!-----------------------------------------------------------------------------------------
module pod_uq_gmm_state_module
    use pod_global, only: DP
    implicit none
    private

    ! 对外暴露的类型 (类似 C++ 的 public class)
    public :: gaussian_component
    public :: uq_gmm_state_type

    !===================================================================
    ! 1. 单个高斯分量定义 (类似 C++ 的 struct 或 class GaussianComponent)
    !===================================================================
    type :: gaussian_component
        real(DP) :: weight                     ! 权重 \mu_i
        real(DP), allocatable :: mean(:)       ! 均值向量 \hat{x}_i
        real(DP), allocatable :: cov(:,:)      ! 协方差矩阵 P_i
    contains
        ! 绑定过程 (类似 C++ 的成员函数 / Member functions)
        procedure :: init => init_gaussian_component
    end type gaussian_component

    !===================================================================
    ! 2. GMM 整体状态定义 (类似 C++ 的 class GMMState)
    !===================================================================
    type :: uq_gmm_state_type
        integer :: n_components                ! 高斯分量总数 N
        integer :: state_dim                   ! 状态维度 (例如 6 维)
        
        ! 核心数据: 存放分量的动态数组 (类似 C++ 的 std::vector<GaussianComponent>)
        type(gaussian_component), allocatable :: components(:)
    contains
        ! 绑定过程
        procedure :: init_from_single => init_gmm_from_single
        procedure :: allocate_components => allocate_gmm_components
    end type uq_gmm_state_type

contains

    !-------------------------------------------------------------------
    ! 实现: 初始化单个高斯分量
    !-------------------------------------------------------------------
    subroutine init_gaussian_component(this, w, m, c)
        class(gaussian_component), intent(inout) :: this   ! class 关键字类似 C++ 的 this 指针
        real(DP), intent(in) :: w
        real(DP), intent(in) :: m(:)
        real(DP), intent(in) :: c(:,:)
        
        this%weight = w
        
        ! 自动分配内存并赋值 (Modern Fortran 的 allocatable 赋值极其方便)
        this%mean = m
        this%cov = c
    end subroutine init_gaussian_component

    !-------------------------------------------------------------------
    ! 实现: 给 GMM 容器分配内存
    !-------------------------------------------------------------------
    subroutine allocate_gmm_components(this, n_comp, dim)
        class(uq_gmm_state_type), intent(inout) :: this
        integer, intent(in) :: n_comp
        integer, intent(in) :: dim
        
        this%n_components = n_comp
        this%state_dim = dim
        
        if (allocated(this%components)) deallocate(this%components)
        allocate(this%components(n_comp))
    end subroutine allocate_gmm_components

    !-------------------------------------------------------------------
    ! 实现: 论文中的初始化逻辑 (将 1 个初始高斯分布均分为 N 个)
    !-------------------------------------------------------------------
    subroutine init_gmm_from_single(this, n_comp, initial_mean, initial_cov)
        class(uq_gmm_state_type), intent(inout) :: this
        integer, intent(in) :: n_comp
        real(DP), intent(in) :: initial_mean(:)
        real(DP), intent(in) :: initial_cov(:,:)
        
        integer :: i
        real(DP) :: initial_weight
        
        ! 1. 分配容器内存
        call this%allocate_components(n_comp, size(initial_mean))
        
        ! 2. 初始权重为 1/N
        initial_weight = 1.0_DP / real(n_comp, DP)
        
        ! 3. 遍历赋值 (每个分量一开始共享相同的均值和协方差)
        do i = 1, n_comp
            call this%components(i)%init(initial_weight, initial_mean, initial_cov)
        end do
        
    end subroutine init_gmm_from_single

end module pod_uq_gmm_state_module