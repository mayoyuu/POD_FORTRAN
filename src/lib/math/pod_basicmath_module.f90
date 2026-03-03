!--------------------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------------------
!> # CAT Basic Mathematics Module
!> 
!> This module provides fundamental mathematical operations for the CAT Fortran
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
!>       all higher-level computations in the CAT Fortran system.
!> 
!> @warning All functions use double precision arithmetic. Ensure consistent
!>          precision throughout the calling code.
!> 
!> @todo Add more matrix operations (inverse, determinant, eigenvalues)
!> @todo Add quaternion operations for attitude representation
!> @todo Add interpolation and approximation functions
!--------------------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------------------

module pod_basicmath
    use pod_global, only: DP
    implicit none
    
    private
    public :: cross_product, dot_product_3d, vector_magnitude, normalize_vector
    public :: norm_vector, norm_matrix
    
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

end module pod_basicmath
