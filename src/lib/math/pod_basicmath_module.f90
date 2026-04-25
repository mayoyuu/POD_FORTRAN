!--------------------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------------------
!> # POD Basic Mathematics Module
!> 
!> This module provides fundamental mathematical operations for the POD Fortran
!> space object monitoring system. It contains basic vector and matrix operations
!> that are commonly used in astrodynamics and numerical computations.
!> 
!> ## Features
!> 
!> - **Vector Operations**: Cross product, dot product, magnitude, normalization
!> - **Matrix Operations**: Basic matrix operations and transformations
!> - **Numerical Utilities**: Mathematical constants and utility functions
!> - **High Precision**: Uses double precision arithmetic for accuracy
!> 
!> ## Mathematical Functions
!> 
!> - **Cross Product**: 3D vector cross product (a × b)
!> - **Dot Product**: 3D vector dot product (a · b)
!> - **Vector Magnitude**: Euclidean norm of 3D vectors
!> - **Vector Normalization**: Unit vector calculation
!> - **Matrix Norm**: Frobenius norm of matrices
!> - **Vector Norm**: General vector norm calculation
!> 
!> ## Dependencies
!> 
!> - **Internal**: `pod_global` for precision constants
!> 
!> ## Author
!> 
!> Zhao Yuhui (PMO, zhaoyuhui@pmo.ac.cn)
!> 
!> ## Version
!> 
!> - **Created**: 2025-09-11
!> - **Updated**: 2025-09-11
!> - **Ref**: 1. Numerical Recipes in Fortran 90: The Art of Parallel Scientific Computing
!>            2. Matrix Computations by Gene H. Golub and Charles F. Van Loan
!>            3. Fundamentals of Astrodynamics and Applications by David A. Vallado
!> 
!> @note This module provides the fundamental mathematical building blocks for
!>       all higher-level computations in the POD Fortran system.
!> 
!> @warning All functions use double precision arithmetic. Ensure consistent
!>          precision throughout the calling code.
!> 
!> @todo Add more matrix operations (inverse, determinant, eigenvalues)
!> @todo Add quaternion operations for attitude representation
!> @todo Add interpolation and approximation functions
!--------------------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------------------

module pod_basicmath_module
    use pod_global, only: DP
    implicit none
    
    !========================================================
    ! LAPACK Interfaces for Matrix Inversion
    !========================================================
    interface

        ! 加入 Cholesky 分解接口
        subroutine dpotrf(uplo, n, a, lda, info)
            import :: DP
            implicit none
            character(len=1), intent(in) :: uplo
            integer, intent(in)          :: n
            real(DP), intent(inout)      :: a(lda, *)
            integer, intent(in)          :: lda
            integer, intent(out)         :: info
        end subroutine dpotrf

        subroutine dgetrf(m, n, a, lda, ipiv, info)
            import :: DP
            implicit none
            integer, intent(in)          :: m, n, lda
            real(DP), intent(inout)      :: a(lda, *)
            integer, intent(out)         :: ipiv(*)
            integer, intent(out)         :: info
        end subroutine dgetrf

        subroutine dgetri(n, a, lda, ipiv, work, lwork, info)
            import :: DP
            implicit none
            integer, intent(in)          :: n, lda, lwork
            real(DP), intent(inout)      :: a(lda, *)
            integer, intent(in)          :: ipiv(*)
            real(DP), intent(out)        :: work(*)
            integer, intent(out)         :: info
        end subroutine dgetri


        subroutine dgeev(jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr, ldvr, work, lwork, info)
            import :: DP
            implicit none
            character, intent(in)            :: jobvl, jobvr
            integer, intent(in)              :: n, lda, ldvl, ldvr, lwork
            real(DP), intent(inout)          :: a(lda, *)
            real(DP), intent(out)            :: wr(*), wi(*), vl(ldvl, *), vr(ldvr, *)
            real(DP), intent(out)            :: work(*)
            integer, intent(out)             :: info
        end subroutine dgeev
    end interface

    private
    public :: cross_product, dot_product_3d, vector_magnitude, normalize_vector
    public :: norm_vector, norm_matrix, inverse_matrix, matrix_determinant, inverse_and_determinant,PI
    public :: eigenvalue_decomposition,dpotrf, dgeev, dgetrf, dgetri

    real(DP), parameter :: PI = 3.14159265358979323846_DP
    
contains

    !> Calculate the cross product of two 3D vectors
    !> 
    !> This function computes the cross product (vector product) of two 3D vectors.
    !> The cross product is perpendicular to both input vectors and has magnitude
    !> equal to the area of the parallelogram formed by the vectors.
    !> 
    !> @param[in] a First vector [ax, ay, az]
    !> @param[in] b Second vector [bx, by, bz]
    !> @return Cross product vector a × b
    !> 
    !> @note The cross product is anti-commutative: a × b = -(b × a)
    !> 
    !> @example
    !> ```fortran
    !> real(DP), dimension(3) :: a, b, result
    !> a = [1.0_DP, 0.0_DP, 0.0_DP]
    !> b = [0.0_DP, 1.0_DP, 0.0_DP]
    !> result = cross_product(a, b)  ! result = [0.0, 0.0, 1.0]
    !> ```
    function cross_product(a, b)
        real(DP), dimension(3), intent(in) :: a, b
        real(DP), dimension(3) :: cross_product
        
        cross_product(1) = a(2) * b(3) - a(3) * b(2)
        cross_product(2) = a(3) * b(1) - a(1) * b(3)
        cross_product(3) = a(1) * b(2) - a(2) * b(1)
    end function cross_product
    
    !> Calculate the dot product of two 3D vectors
    !> 
    !> This function computes the dot product (scalar product) of two 3D vectors.
    !> The dot product is a scalar value equal to the product of the magnitudes
    !> of the vectors and the cosine of the angle between them.
    !> 
    !> @param[in] a First vector [ax, ay, az]
    !> @param[in] b Second vector [bx, by, bz]
    !> @return Dot product scalar a · b
    !> 
    !> @note The dot product is commutative: a · b = b · a
    !> 
    !> @example
    !> ```fortran
    !> real(DP), dimension(3) :: a, b
    !> real(DP) :: result
    !> a = [1.0_DP, 2.0_DP, 3.0_DP]
    !> b = [4.0_DP, 5.0_DP, 6.0_DP]
    !> result = dot_product_3d(a, b)  ! result = 32.0
    !> ```
    function dot_product_3d(a, b)
        real(DP), dimension(3), intent(in) :: a, b
        real(DP) :: dot_product_3d
        
        dot_product_3d = a(1) * b(1) + a(2) * b(2) + a(3) * b(3)
    end function dot_product_3d
    
    !> Calculate the magnitude (norm) of a 3D vector
    !> 
    !> This function computes the Euclidean norm (magnitude) of a 3D vector.
    !> The magnitude is the length of the vector in 3D space.
    !> 
    !> @param[in] v Input vector [vx, vy, vz]
    !> @return Vector magnitude ||v||
    !> 
    !> @note The magnitude is always non-negative
    !> 
    !> @example
    !> ```fortran
    !> real(DP), dimension(3) :: v
    !> real(DP) :: magnitude
    !> v = [3.0_DP, 4.0_DP, 0.0_DP]
    !> magnitude = vector_magnitude(v)  ! magnitude = 5.0
    !> ```
    function vector_magnitude(v)
        real(DP), dimension(3), intent(in) :: v
        real(DP) :: vector_magnitude
        
        vector_magnitude = sqrt(v(1)**2 + v(2)**2 + v(3)**2)
    end function vector_magnitude
    
    !> Normalize a 3D vector to unit length
    !> 
    !> This function normalizes a 3D vector to have unit magnitude.
    !> If the input vector has zero magnitude, returns a zero vector.
    !> 
    !> @param[in] v Input vector [vx, vy, vz]
    !> @return Normalized unit vector
    !> 
    !> @note The normalized vector has magnitude 1.0 (or 0.0 if input is zero)
    !> 
    !> @example
    !> ```fortran
    !> real(DP), dimension(3) :: v, unit_v
    !> v = [3.0_DP, 4.0_DP, 0.0_DP]
    !> unit_v = normalize_vector(v)  ! unit_v = [0.6, 0.8, 0.0]
    !> ```
    function normalize_vector(v)
        real(DP), dimension(3), intent(in) :: v
        real(DP), dimension(3) :: normalize_vector
        real(DP) :: mag
        
        mag = vector_magnitude(v)
        if (mag > epsilon(mag)) then
            normalize_vector = v / mag
        else
            normalize_vector = 0.0_DP
        end if
    end function normalize_vector
    
    !> Calculate the norm of a general vector
    !> 
    !> This function computes the Euclidean norm of a vector of any dimension.
    !> 
    !> @param[in] vector Input vector of any size
    !> @return Vector norm ||vector||
    !> 
    !> @note This is a general-purpose function for vectors of any dimension
    !> 
    !> @example
    !> ```fortran
    !> real(DP), dimension(5) :: v
    !> real(DP) :: norm_val
    !> v = [1.0_DP, 2.0_DP, 3.0_DP, 4.0_DP, 5.0_DP]
    !> norm_val = norm_vector(v)  ! norm_val = sqrt(55)
    !> ```
    function norm_vector(vector)
        real(DP), dimension(:), intent(in) :: vector
        real(DP) :: norm_vector
        norm_vector = sqrt(sum(vector**2))
    end function norm_vector
    
    !> Calculate the Frobenius norm of a matrix
    !> 
    !> This function computes the Frobenius norm (Euclidean norm) of a matrix.
    !> The Frobenius norm is the square root of the sum of squares of all elements.
    !> 
    !> @param[in] matrix Input matrix of any size
    !> @return Matrix Frobenius norm ||matrix||_F
    !> 
    !> @note The Frobenius norm is the matrix equivalent of the Euclidean vector norm
    !> 
    !> @example
    !> ```fortran
    !> real(DP), dimension(2,2) :: A
    !> real(DP) :: norm_val
    !> A = reshape([1.0_DP, 2.0_DP, 3.0_DP, 4.0_DP], [2,2])
    !> norm_val = norm_matrix(A)  ! norm_val = sqrt(30)
    !> ```
    function norm_matrix(matrix)
        real(DP), dimension(:,:), intent(in) :: matrix
        real(DP) :: norm_matrix
        norm_matrix = sqrt(sum(matrix**2))
    end function norm_matrix

    !> Calculate the inverse of a square matrix using LAPACK
    !> 
    !> @param[in]  A        Input square matrix
    !> @param[out] A_inv    Inverse of the matrix
    !> @param[out] info     Status code (0 = success, >0 = singular matrix)
    subroutine inverse_matrix(A, A_inv, info)
        real(DP), dimension(:,:), intent(in)  :: A
        real(DP), dimension(:,:), intent(out) :: A_inv
        integer, intent(out)                  :: info
        
        integer :: n, lwork
        integer, allocatable :: ipiv(:)
        real(DP), allocatable :: work(:)
        real(DP) :: work_query(1)
        
        n = size(A, 1)
        A_inv = A  ! LAPACK routines overwrite the input matrix
        
        allocate(ipiv(n))
        
        ! 1. LU Factorization
        call dgetrf(n, n, A_inv, n, ipiv, info)
        if (info /= 0) then
            deallocate(ipiv)
            return
        end if
        
        ! 2. Workspace query for optimal performance
        call dgetri(n, A_inv, n, ipiv, work_query, -1, info)
        lwork = int(work_query(1))
        allocate(work(lwork))
        
        ! 3. Compute Inverse
        call dgetri(n, A_inv, n, ipiv, work, lwork, info)
        
        deallocate(ipiv)
        deallocate(work)
    end subroutine inverse_matrix

    !> Calculate the determinant of a square matrix using LAPACK LU factorization
    !> 
    !> @param[in]  A        Input square matrix
    !> @param[out] det      Determinant of the matrix
    !> @param[out] info     Status code (0 = success, >0 = singular matrix)
    subroutine matrix_determinant(A, det, info)
        real(DP), dimension(:,:), intent(in)  :: A
        real(DP), intent(out)                 :: det
        integer, intent(out)                  :: info
        
        integer :: n, i
        real(DP), dimension(:,:), allocatable :: LU
        integer, dimension(:), allocatable    :: ipiv
        real(DP) :: det_sign
        
        n = size(A, 1)
        allocate(LU(n, n), ipiv(n))
        
        LU = A ! 复制输入矩阵，因为 dgetrf 会覆盖它
        
        ! 1. 进行 LU 分解
        call dgetrf(n, n, LU, n, ipiv, info)
        
        if (info > 0) then
            ! 矩阵是奇异的，行列式为 0
            det = 0.0_DP
            deallocate(LU, ipiv)
            return
        end if
        
        ! 2. 从 U 的对角线和行交换记录中计算行列式
        det = 1.0_DP
        det_sign = 1.0_DP
        
        do i = 1, n
            det = det * LU(i, i)
            ! LAPACK 的 ipiv 记录了置换操作，如果不等于当前索引 i，说明发生了行交换
            if (ipiv(i) /= i) then
                det_sign = -det_sign
            end if
        end do
        
        det = det * det_sign
        
        deallocate(LU, ipiv)
    end subroutine matrix_determinant

    !> Calculate both the inverse and determinant of a square matrix simultaneously
    !> Highly optimized for multivariate Gaussian PDF calculations where both are needed.
    !> 
    !> @param[in]  A        Input square matrix (e.g., Covariance matrix)
    !> @param[out] A_inv    Inverse of the matrix
    !> @param[out] det      Determinant of the matrix
    !> @param[out] info     Status code (0 = success, >0 = singular matrix)
    subroutine inverse_and_determinant(A, A_inv, det, info)
        real(DP), dimension(:,:), intent(in)  :: A
        real(DP), dimension(:,:), intent(out) :: A_inv
        real(DP), intent(out)                 :: det
        integer, intent(out)                  :: info
        
        integer :: n, i, lwork
        integer, allocatable :: ipiv(:)
        real(DP), allocatable :: work(:)
        real(DP) :: work_query(1)
        real(DP) :: det_sign
        
        n = size(A, 1)
        A_inv = A  ! 复制输入矩阵
        
        allocate(ipiv(n))
        
        ! 1. LU 分解 (只做一次!)
        call dgetrf(n, n, A_inv, n, ipiv, info)
        
        if (info > 0) then
            det = 0.0_DP
            deallocate(ipiv)
            return
        end if
        
        ! 2. 趁着 A_inv 里面还装着 U 矩阵，先把行列式算出来
        det = 1.0_DP
        det_sign = 1.0_DP
        do i = 1, n
            det = det * A_inv(i, i)
            if (ipiv(i) /= i) then
                det_sign = -det_sign
            end if
        end do
        det = det * det_sign
        
        ! 3. 查询求逆最优工作空间大小
        call dgetri(n, A_inv, n, ipiv, work_query, -1, info)
        lwork = int(work_query(1))
        allocate(work(lwork))
        
        ! 4. 原地计算逆矩阵
        call dgetri(n, A_inv, n, ipiv, work, lwork, info)
        
        deallocate(ipiv)
        deallocate(work)
    end subroutine inverse_and_determinant

    !> Calculate eigenvalues and eigenvectors of a general real square matrix
    !> 
    !> Uses LAPACK's DGEEV routine.
    !> 
    !> @param[in]  A       Input square matrix (N x N)
    !> @param[out] wr      Real parts of the eigenvalues
    !> @param[out] wi      Imaginary parts of the eigenvalues
    !> @param[out] vr      Right eigenvectors (optional)
    !> @param[out] info    Status code (0 = success, <0 = illegal value, >0 = convergence failed)
    subroutine eigenvalue_decomposition(A, wr, wi, vr, info)
        real(DP), dimension(:,:), intent(in)           :: A
        real(DP), dimension(:), intent(out)             :: wr, wi
        real(DP), dimension(:,:), intent(out), optional :: vr
        integer, intent(out)                            :: info
        
        integer :: n, lwork, ldvr
        character :: jobvr
        real(DP), allocatable :: work(:), A_tmp(:,:)
        real(DP), dimension(1,1) :: dummy_vl ! Left eigenvectors not requested
        real(DP) :: work_query(1)
        
        n = size(A, 1)
        allocate(A_tmp(n, n))
        A_tmp = A ! DGEEV overwrites the input matrix
        
        ! Determine if right eigenvectors are requested
        if (present(vr)) then
            jobvr = 'V'
            ldvr = n
        else
            jobvr = 'N'
            ldvr = 1
        end if
        
        ! 1. Workspace query
        call dgeev('N', jobvr, n, A_tmp, n, wr, wi, dummy_vl, 1, vr, ldvr, work_query, -1, info)
        lwork = int(work_query(1))
        allocate(work(lwork))
        
        ! 2. Compute Eigenvalues (and optionally Eigenvectors)
        if (present(vr)) then
            call dgeev('N', jobvr, n, A_tmp, n, wr, wi, dummy_vl, 1, vr, ldvr, work, lwork, info)
        else
            ! Use a dummy array for vr if not requested
            block
                real(DP) :: dummy_vr(1,1)
                call dgeev('N', jobvr, n, A_tmp, n, wr, wi, dummy_vl, 1, dummy_vr, 1, work, lwork, info)
            end block
        end if
        
        deallocate(work, A_tmp)
    end subroutine eigenvalue_decomposition

end module pod_basicmath_module
