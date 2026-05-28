!> 数据结构定义：不确定性量化状态
module pod_uq_state_module
    use pod_global, only: DP
    implicit none
    private
    public :: uq_state_type

    type :: uq_state_type
        ! 注意：Fortran 中强烈建议 (状态维度, 粒子数) 的列优先布局
        real(DP), allocatable :: samples(:,:)  ! size: (dim, n_particles)
        real(DP), allocatable :: mean(:)       ! size: (dim)
        real(DP), allocatable :: cov(:,:)      ! size: (dim, dim)
    contains
        ! 重命名为 allocate_memory，明确它只做内存分配
        procedure :: allocate_memory => state_allocate_memory
        procedure :: deallocate_memory => state_deallocate_memory 
        procedure :: compute_moments => state_compute_moments
        procedure :: compute_higher_moments => state_compute_higher_moments
    end type uq_state_type

contains

    ! 只负责分配内存
    subroutine state_allocate_memory(this, dim, n_particles)
        class(uq_state_type), intent(inout) :: this
        integer, intent(in) :: dim, n_particles
        
        if (allocated(this%samples)) deallocate(this%samples)
        allocate(this%samples(dim, n_particles))
        this%samples = 0.0_DP
    end subroutine state_allocate_memory

    subroutine state_deallocate_memory(this)
        class(uq_state_type), intent(inout) :: this
        
        if (allocated(this%samples)) deallocate(this%samples)
        if (allocated(this%mean))    deallocate(this%mean)
        if (allocated(this%cov))     deallocate(this%cov)
    end subroutine state_deallocate_memory
    
    ! 计算均值和协方差
    subroutine state_compute_moments(this)
        class(uq_state_type), intent(inout) :: this
        integer :: dim, n_particles, i, k
        real(DP), allocatable :: diff(:,:)
        
        if (.not. allocated(this%samples)) then
            write(*,*) '[Error] UQ State: Samples not allocated.'
            return
        end if
        
        dim = size(this%samples, 1)
        n_particles = size(this%samples, 2)
        
        if (allocated(this%mean)) deallocate(this%mean)
        allocate(this%mean(dim))
        
        if (allocated(this%cov)) deallocate(this%cov)
        allocate(this%cov(dim, dim))
        
        ! 1. 计算均值
        do i = 1, dim
            this%mean(i) = sum(this%samples(i, :)) / real(n_particles, DP)
        end do
        
        ! 2. 计算协方差
        ! 开辟一个临时矩阵存储每个粒子与均值的偏差
        allocate(diff(dim, n_particles))
        
        ! 使用 Fortran 的 spread 函数将一维 mean 扩展为二维矩阵进行逐元素相减
        diff = this%samples - spread(this%mean, dim=2, ncopies=n_particles)
        
        ! 利用内置矩阵乘法 matmul 计算协方差矩阵 (dim x n_particles) * (n_particles x dim)
        ! 此操作如果链接了 OpenBLAS/MKL，底层会自动多线程加速
        this%cov = matmul(diff, transpose(diff)) / real(n_particles - 1, DP)
        
        deallocate(diff)
    end subroutine state_compute_moments

    subroutine state_compute_higher_moments(this, skewness, kurtosis)
        class(uq_state_type), intent(in) :: this
        real(DP), allocatable, intent(out) :: skewness(:)
        real(DP), allocatable, intent(out) :: kurtosis(:)

        integer :: dim, n_particles, i, j
        real(DP), allocatable :: diff(:), sigma(:)

        if (.not. allocated(this%samples)) then
            write(*,*) '[Error] UQ State: Samples not allocated for higher moments.'
            return
        end if
        if (.not. allocated(this%mean)) then
            write(*,*) '[Error] UQ State: Mean not computed. Call compute_moments first.'
            return
        end if

        dim = size(this%samples, 1)
        n_particles = size(this%samples, 2)

        allocate(skewness(dim), kurtosis(dim))
        allocate(diff(n_particles), sigma(dim))

        do i = 1, dim
            sigma(i) = sqrt(this%cov(i, i))
        end do

        do i = 1, dim
            do j = 1, n_particles
                diff(j) = this%samples(i, j) - this%mean(i)
            end do

            if (sigma(i) > 0.0_DP) then
                skewness(i) = sum(diff**3) / real(n_particles, DP) / (sigma(i)**3)
                kurtosis(i) = sum(diff**4) / real(n_particles, DP) / (sigma(i)**4)
            else
                skewness(i) = 0.0_DP
                kurtosis(i) = 0.0_DP
            end if
        end do

        deallocate(diff, sigma)
    end subroutine state_compute_higher_moments

end module pod_uq_state_module