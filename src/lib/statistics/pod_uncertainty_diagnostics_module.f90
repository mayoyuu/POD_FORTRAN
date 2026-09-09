!> @file pod_uncertainty_diagnostics_module.f90
!! @brief Deterministic probability quadrature and covariance diagnostics.
!!
!! The module deliberately accepts normalized probability weights.  It is
!! independent of any particular propagator and is reusable by SRP and other
!! low-dimensional uncertainty studies.
module pod_uncertainty_diagnostics_module
    use, intrinsic :: iso_fortran_env, only: int64
    use pod_global, only: DP
    use pod_basicmath_module, only: eigenvalue_decomposition
    implicit none
    private

    public :: gauss_legendre_probability_rule
    public :: gauss_normal_probability_rule
    public :: compute_weighted_moments
    public :: compute_covariance_axes
    public :: compute_effective_rank
    public :: generate_antithetic_standard_normal_samples

contains

    !> Generate a reproducible antithetic standard-normal sample ensemble.
    !!
    !! The first half uses midpoint-stratified standard-normal quantiles.  A
    !! local Park-Miller stream independently permutes those quantiles in every
    !! component, and the second half is the exact negative of the first.  This
    !! gives stable tail/fourth-moment coverage, exactly zero finite-ensemble
    !! mean, and does not change Fortran's global random-number state.  Each
    !! component is finally rescaled to unit second moment; the moment matching
    !! leaves every antithetic pair intact.
    subroutine generate_antithetic_standard_normal_samples(samples,seed,status)
        real(DP), intent(out) :: samples(:,:)
        integer, intent(in) :: seed
        integer, intent(out) :: status

        integer(int64), parameter :: LCG_A=48271_int64
        integer(int64), parameter :: LCG_M=2147483647_int64
        integer(int64) :: rng_state
        integer,allocatable :: permutation(:)
        real(DP),allocatable :: normal_grid(:)
        integer :: n_dim,n_samples,n_half,column,component,swap_column,temp_index
        real(DP) :: u1,radius

        status=0
        samples=0.0_DP
        n_dim=size(samples,1)
        n_samples=size(samples,2)
        if(n_dim<1 .or. n_samples<2 .or. mod(n_samples,2)/=0) then
            status=-1
            return
        end if
        if(seed<=0) then
            status=-2
            return
        end if

        n_half=n_samples/2
        rng_state=mod(int(seed,int64),LCG_M-1_int64)+1_int64
        allocate(permutation(n_half),normal_grid(n_half))
        do column=1,n_half
            normal_grid(column)=inverse_standard_normal( &
                (real(column,DP)-0.5_DP)/real(n_half,DP))
        end do
        do component=1,n_dim
            permutation=[(column,column=1,n_half)]
            do column=n_half,2,-1
                call next_uniform(rng_state,u1)
                swap_column=1+int(u1*real(column,DP))
                temp_index=permutation(column)
                permutation(column)=permutation(swap_column)
                permutation(swap_column)=temp_index
            end do
            samples(component,1:n_half)=normal_grid(permutation)
        end do
        samples(:,n_half+1:n_samples)=-samples(:,1:n_half)

        do component=1,n_dim
            radius=sqrt(sum(samples(component,:)**2)/real(n_samples,DP))
            if(radius<=tiny(1.0_DP)) then
                samples=0.0_DP
                status=-3
                return
            end if
            samples(component,:)=samples(component,:)/radius
        end do
        deallocate(permutation,normal_grid)

    contains

        subroutine next_uniform(state,value)
            integer(int64), intent(inout) :: state
            real(DP), intent(out) :: value

            state=mod(LCG_A*state,LCG_M)
            value=real(state,DP)/real(LCG_M,DP)
        end subroutine next_uniform

        function inverse_standard_normal(probability) result(value)
            real(DP),intent(in) :: probability
            real(DP) :: value,q,r,cdf,pdf
            real(DP),parameter :: P_LOW=0.02425_DP,P_HIGH=1.0_DP-P_LOW
            real(DP),parameter :: SQRT_TWO=sqrt(2.0_DP)
            real(DP),parameter :: SQRT_TWO_PI=sqrt(2.0_DP*acos(-1.0_DP))
            real(DP),parameter :: A(6)=[-3.969683028665376e1_DP, &
                2.209460984245205e2_DP,-2.759285104469687e2_DP, &
                1.383577518672690e2_DP,-3.066479806614716e1_DP, &
                2.506628277459239_DP]
            real(DP),parameter :: B(5)=[-5.447609879822406e1_DP, &
                1.615858368580409e2_DP,-1.556989798598866e2_DP, &
                6.680131188771972e1_DP,-1.328068155288572e1_DP]
            real(DP),parameter :: C(6)=[-7.784894002430293e-3_DP, &
                -3.223964580411365e-1_DP,-2.400758277161838_DP, &
                -2.549732539343734_DP,4.374664141464968_DP, &
                2.938163982698783_DP]
            real(DP),parameter :: D(4)=[7.784695709041462e-3_DP, &
                3.224671290700398e-1_DP,2.445134137142996_DP, &
                3.754408661907416_DP]

            if(probability<P_LOW) then
                q=sqrt(-2.0_DP*log(probability))
                value=(((((C(1)*q+C(2))*q+C(3))*q+C(4))*q+C(5))*q+C(6))/ &
                    ((((D(1)*q+D(2))*q+D(3))*q+D(4))*q+1.0_DP)
            else if(probability<=P_HIGH) then
                q=probability-0.5_DP
                r=q*q
                value=(((((A(1)*r+A(2))*r+A(3))*r+A(4))*r+A(5))*r+A(6))*q/ &
                    (((((B(1)*r+B(2))*r+B(3))*r+B(4))*r+B(5))*r+1.0_DP)
            else
                q=sqrt(-2.0_DP*log(1.0_DP-probability))
                value=-(((((C(1)*q+C(2))*q+C(3))*q+C(4))*q+C(5))*q+C(6))/ &
                    ((((D(1)*q+D(2))*q+D(3))*q+D(4))*q+1.0_DP)
            end if

            cdf=0.5_DP*(1.0_DP+erf(value/SQRT_TWO))
            pdf=exp(-0.5_DP*value*value)/SQRT_TWO_PI
            value=value-(cdf-probability)/pdf
        end function inverse_standard_normal

    end subroutine generate_antithetic_standard_normal_samples

    !> Construct an n-point Gauss-Legendre rule for probability U(-1,1).
    subroutine gauss_legendre_probability_rule(nodes, weights, status)
        real(DP), intent(out) :: nodes(:), weights(:)
        integer, intent(out) :: status
        real(DP), allocatable :: jacobi(:,:), wr(:), wi(:), vr(:,:)
        integer :: n, k, info

        status = 0
        nodes = 0.0_DP
        weights = 0.0_DP
        n = size(nodes)
        if (n < 1 .or. size(weights) /= n) then
            status = -1
            return
        end if

        allocate(jacobi(n,n), wr(n), wi(n), vr(n,n))
        jacobi = 0.0_DP
        do k = 1, n-1
            jacobi(k,k+1) = real(k,DP)/sqrt(4.0_DP*real(k,DP)**2-1.0_DP)
            jacobi(k+1,k) = jacobi(k,k+1)
        end do
        call eigenvalue_decomposition(jacobi, wr, wi, vr, info)
        if (info /= 0 .or. .not. eigenvalues_are_real(wr,wi)) then
            status = -2
        else
            nodes = wr
            weights = vr(1,:)**2
            call sort_rule_ascending(nodes, weights)
            weights = weights/sum(weights)
        end if
        deallocate(jacobi,wr,wi,vr)
    end subroutine gauss_legendre_probability_rule

    !> Construct an n-point Gauss-Hermite rule for N(0,sigma^2).
    !!
    !! The Jacobi matrix uses the probabilists' Hermite recurrence, so its
    !! eigenvalues are standard-normal nodes and first-row eigenvector squares
    !! are already normalized probability weights.
    subroutine gauss_normal_probability_rule(sigma, nodes, weights, status)
        real(DP), intent(in) :: sigma
        real(DP), intent(out) :: nodes(:), weights(:)
        integer, intent(out) :: status
        real(DP), allocatable :: jacobi(:,:), wr(:), wi(:), vr(:,:)
        integer :: n, k, info

        status = 0
        nodes = 0.0_DP
        weights = 0.0_DP
        n = size(nodes)
        if (n < 1 .or. size(weights) /= n .or. sigma < 0.0_DP) then
            status = -1
            return
        end if

        allocate(jacobi(n,n), wr(n), wi(n), vr(n,n))
        jacobi = 0.0_DP
        do k = 1, n-1
            jacobi(k,k+1) = sqrt(real(k,DP))
            jacobi(k+1,k) = jacobi(k,k+1)
        end do
        call eigenvalue_decomposition(jacobi, wr, wi, vr, info)
        if (info /= 0 .or. .not. eigenvalues_are_real(wr,wi)) then
            status = -2
        else
            nodes = sigma*wr
            weights = vr(1,:)**2
            call sort_rule_ascending(nodes, weights)
            weights = weights/sum(weights)
        end if
        deallocate(jacobi,wr,wi,vr)
    end subroutine gauss_normal_probability_rule

    !> Compute weighted mean, covariance, skewness and excess kurtosis.
    subroutine compute_weighted_moments(samples, weights, mean_value, covariance, &
                                        skewness, excess_kurtosis, status)
        real(DP), intent(in) :: samples(:,:), weights(:)
        real(DP), intent(out) :: mean_value(:), covariance(:,:)
        real(DP), intent(out) :: skewness(:), excess_kurtosis(:)
        integer, intent(out) :: status

        real(DP), allocatable :: normalized_weights(:), deviation(:)
        real(DP) :: weight_sum, variance, sigma, scale
        integer :: n_state, n_nodes, i, j

        status = 0
        mean_value = 0.0_DP
        covariance = 0.0_DP
        skewness = 0.0_DP
        excess_kurtosis = 0.0_DP
        n_state = size(samples,1)
        n_nodes = size(samples,2)
        if (n_state < 1 .or. n_nodes < 1 .or. size(weights) /= n_nodes .or. &
            size(mean_value) /= n_state .or. size(skewness) /= n_state .or. &
            size(excess_kurtosis) /= n_state .or. &
            size(covariance,1) /= n_state .or. size(covariance,2) /= n_state) then
            status = -1
            return
        end if
        if (any(weights < 0.0_DP)) then
            status = -2
            return
        end if
        weight_sum = sum(weights)
        if (weight_sum <= 0.0_DP) then
            status = -3
            return
        end if

        allocate(normalized_weights(n_nodes), deviation(n_state))
        normalized_weights = weights/weight_sum
        mean_value = matmul(samples,normalized_weights)
        do j = 1, n_nodes
            deviation = samples(:,j)-mean_value
            do i = 1, n_state
                covariance(:,i) = covariance(:,i) + &
                    normalized_weights(j)*deviation*deviation(i)
            end do
        end do
        covariance = 0.5_DP*(covariance+transpose(covariance))

        do i = 1, n_state
            variance = max(covariance(i,i),0.0_DP)
            sigma = sqrt(variance)
            scale = max(1.0_DP,abs(mean_value(i)),maxval(abs(samples(i,:))))
            if (sigma > 100.0_DP*epsilon(1.0_DP)*scale) then
                skewness(i) = sum(normalized_weights* &
                    ((samples(i,:)-mean_value(i))/sigma)**3)
                excess_kurtosis(i) = sum(normalized_weights* &
                    ((samples(i,:)-mean_value(i))/sigma)**4)-3.0_DP
            end if
        end do
        deallocate(normalized_weights,deviation)
    end subroutine compute_weighted_moments

    !> Return ordered eigenvalues/eigenvectors and one-sigma semi-axis scales.
    subroutine compute_covariance_axes(covariance, eigenvalues, eigenvectors, &
                                       axis_sigma, status)
        real(DP), intent(in) :: covariance(3,3)
        real(DP), intent(out) :: eigenvalues(3), eigenvectors(3,3), axis_sigma(3)
        integer, intent(out) :: status

        real(DP) :: symmetric(3,3), wr(3), wi(3), vr(3,3)
        real(DP) :: negative_tolerance, vector_norm
        integer :: info, i

        status = 0
        symmetric = 0.5_DP*(covariance+transpose(covariance))
        call eigenvalue_decomposition(symmetric,wr,wi,vr,info)
        if (info /= 0 .or. .not. eigenvalues_are_real(wr,wi)) then
            status = -1
            eigenvalues = 0.0_DP
            eigenvectors = 0.0_DP
            axis_sigma = 0.0_DP
            return
        end if

        call sort_eigenpairs_descending(wr,vr)
        negative_tolerance = 1000.0_DP*epsilon(1.0_DP)* &
            max(maxval(abs(wr)),tiny(1.0_DP))
        if (minval(wr) < -negative_tolerance) then
            status = -2
            eigenvalues = wr
            eigenvectors = vr
            axis_sigma = 0.0_DP
            return
        end if

        eigenvalues = max(wr,0.0_DP)
        eigenvectors = vr
        do i = 1, 3
            vector_norm = sqrt(dot_product(eigenvectors(:,i),eigenvectors(:,i)))
            if (vector_norm <= tiny(1.0_DP)) then
                status = -3
                return
            end if
            eigenvectors(:,i) = eigenvectors(:,i)/vector_norm
        end do
        axis_sigma = sqrt(eigenvalues)
    end subroutine compute_covariance_axes

    !> Count covariance eigenvalues significant relative to the largest one.
    subroutine compute_effective_rank(covariance, relative_tolerance, rank_value, status)
        real(DP), intent(in) :: covariance(:,:), relative_tolerance
        integer, intent(out) :: rank_value, status

        real(DP), allocatable :: symmetric(:,:), wr(:), wi(:)
        real(DP) :: largest
        integer :: n, info

        status = 0
        rank_value = 0
        n = size(covariance,1)
        if (n < 1 .or. size(covariance,2) /= n .or. relative_tolerance < 0.0_DP) then
            status = -1
            return
        end if
        allocate(symmetric(n,n),wr(n),wi(n))
        symmetric = 0.5_DP*(covariance+transpose(covariance))
        call eigenvalue_decomposition(symmetric,wr,wi,info=info)
        if (info /= 0 .or. .not. eigenvalues_are_real(wr,wi)) then
            status = -2
        else
            largest = max(maxval(wr),0.0_DP)
            if (largest > 0.0_DP) then
                rank_value = count(wr > relative_tolerance*largest)
            end if
        end if
        deallocate(symmetric,wr,wi)
    end subroutine compute_effective_rank

    pure logical function eigenvalues_are_real(wr,wi)
        real(DP), intent(in) :: wr(:), wi(:)
        real(DP) :: tolerance
        tolerance = 1000.0_DP*epsilon(1.0_DP)* &
                    max(maxval(abs(wr)),tiny(1.0_DP))
        eigenvalues_are_real = maxval(abs(wi)) <= tolerance
    end function eigenvalues_are_real

    pure subroutine sort_rule_ascending(nodes,weights)
        real(DP), intent(inout) :: nodes(:), weights(:)
        real(DP) :: node_value, weight_value
        integer :: i, j
        do i = 2, size(nodes)
            node_value = nodes(i)
            weight_value = weights(i)
            j = i-1
            do while (j >= 1)
                if (nodes(j) <= node_value) exit
                nodes(j+1) = nodes(j)
                weights(j+1) = weights(j)
                j = j-1
            end do
            nodes(j+1) = node_value
            weights(j+1) = weight_value
        end do
    end subroutine sort_rule_ascending

    pure subroutine sort_eigenpairs_descending(values,vectors)
        real(DP), intent(inout) :: values(:), vectors(:,:)
        real(DP) :: value, vector(size(vectors,1))
        integer :: i, j
        do i = 2, size(values)
            value = values(i)
            vector = vectors(:,i)
            j = i-1
            do while (j >= 1)
                if (values(j) >= value) exit
                values(j+1) = values(j)
                vectors(:,j+1) = vectors(:,j)
                j = j-1
            end do
            values(j+1) = value
            vectors(:,j+1) = vector
        end do
    end subroutine sort_eigenpairs_descending

end module pod_uncertainty_diagnostics_module
