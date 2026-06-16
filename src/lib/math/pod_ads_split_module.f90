!> ADS (Automatic Domain Splitting) core module
!> Provides SplittingHistory, Patch, Manifold types and operations
!> Translated from C++ reference: ADS_cpp_Core_file/
module pod_ads_split_module
    use pod_global, only: DP
    use pod_dace_classes, only: AlgebraicVector, DA, da_var, da_estim_norm
    implicit none
    private

    ! =========================================================================
    ! Types
    ! =========================================================================
    type :: splitting_history_type
        integer, allocatable :: entries(:)
    end type

    type :: patch_type
        type(AlgebraicVector) :: da_vec
        type(splitting_history_type) :: history
    end type

    type :: manifold_type
        type(patch_type), allocatable :: patches(:)
        integer :: n_patches
    end type

    public :: splitting_history_type, patch_type, manifold_type
    public :: sh_push, sh_pop, sh_count, sh_replay, sh_center, sh_width, sh_contain, sh_map_point

contains

    ! =========================================================================
    ! SplittingHistory: sh_push
    ! =========================================================================
    subroutine sh_push(history, dir)
        type(splitting_history_type), intent(inout) :: history
        integer, intent(in) :: dir
        integer, allocatable :: tmp(:)
        integer :: n
        if (.not. allocated(history%entries)) then
            allocate(history%entries(1))
            history%entries(1) = dir
        else
            n = size(history%entries)
            allocate(tmp(n+1))
            tmp(1:n) = history%entries
            tmp(n+1) = dir
            call move_alloc(tmp, history%entries)
        end if
    end subroutine sh_push

    ! =========================================================================
    ! SplittingHistory: sh_pop
    ! =========================================================================
    subroutine sh_pop(history)
        type(splitting_history_type), intent(inout) :: history
        integer, allocatable :: tmp(:)
        integer :: n
        if (.not. allocated(history%entries)) return
        n = size(history%entries)
        if (n <= 1) then
            deallocate(history%entries)
        else
            allocate(tmp(n-1))
            tmp = history%entries(1:n-1)
            call move_alloc(tmp, history%entries)
        end if
    end subroutine sh_pop

    ! =========================================================================
    ! SplittingHistory: sh_count
    ! =========================================================================
    integer function sh_count(history, n) result(c)
        type(splitting_history_type), intent(in) :: history
        integer, intent(in) :: n
        integer :: i
        c = 0
        if (.not. allocated(history%entries)) return
        if (n == 0) then
            c = size(history%entries)
        else
            do i = 1, size(history%entries)
                if (abs(history%entries(i)) == n) c = c + 1
            end do
        end if
    end function sh_count

    ! =========================================================================
    ! SplittingHistory: sh_center
    ! =========================================================================
    function sh_center(history) result(c)
        type(splitting_history_type), intent(in) :: history
        real(DP) :: c(6)
        real(DP) :: w(6)
        integer :: i, n, sgn
        real(DP) :: half_w
        w = 2.0_DP
        c = 0.0_DP
        if (.not. allocated(history%entries)) return
        do i = 1, size(history%entries)
            n = abs(history%entries(i)) - 1
            sgn = history%entries(i) / abs(history%entries(i))
            w(n+1) = 0.5_DP * w(n+1)
            half_w = abs(w(n+1))
            c(n+1) = c(n+1) + 0.5_DP * real(sgn, DP) * half_w
        end do
    end function sh_center

    ! =========================================================================
    ! SplittingHistory: sh_width
    ! =========================================================================
    function sh_width(history) result(w)
        type(splitting_history_type), intent(in) :: history
        real(DP) :: w(6)
        integer :: i, n
        w = 2.0_DP
        if (.not. allocated(history%entries)) return
        do i = 1, size(history%entries)
            n = abs(history%entries(i)) - 1
            w(n+1) = 0.5_DP * abs(w(n+1))
        end do
    end function sh_width

    ! =========================================================================
    ! SplittingHistory: sh_contain
    ! =========================================================================
    logical function sh_contain(history, pt) result(ok)
        type(splitting_history_type), intent(in) :: history
        real(DP), intent(in) :: pt(6)
        real(DP) :: c(6), w(6)
        integer :: i
        c = sh_center(history)
        w = sh_width(history)
        ok = .true.
        do i = 1, 6
            if (abs(pt(i) - c(i)) > 0.5_DP * w(i)) then
                ok = .false.
                return
            end if
        end do
    end function sh_contain

    ! =========================================================================
    ! SplittingHistory: sh_map_point
    ! =========================================================================
    subroutine sh_map_point(history, pt)
        type(splitting_history_type), intent(in) :: history
        real(DP), intent(inout) :: pt(6)
        integer :: i, n
        if (.not. allocated(history%entries)) return
        do i = 1, size(history%entries)
            n = abs(history%entries(i))
            if (history%entries(i) > 0) then
                pt(n) = 2.0_DP * pt(n) - 1.0_DP   ! right split
            else
                pt(n) = 2.0_DP * pt(n) + 1.0_DP   ! left split
            end if
        end do
    end subroutine sh_map_point

    ! =========================================================================
    ! SplittingHistory: sh_replay
    !
    ! NOTE: The canonical C++ implementation composes the DA map with the split
    ! transformations by performing DA-to-DA substitution (evaluating a DA
    ! polynomial at DA-valued coordinates). The current pod_dace_classes API
    ! does not expose this operation: AlgebraicVector%eval only accepts
    ! (integer, real) for single-variable real substitution or (real(:)) for
    ! full real-vector evaluation. There is no DA-valued eval / composition.
    !
    ! This implementation builds the identity vector x(i)=da_var(i) and applies
    ! the affine split transformations to x, but the actual composition step
    ! (obj = obj.eval(x) in the C++ reference) is not available.
    ! =========================================================================
    subroutine sh_replay(history, obj)
        type(splitting_history_type), intent(in) :: history
        type(AlgebraicVector), intent(inout) :: obj
        type(AlgebraicVector) :: x
        type(DA) :: tmp_da
        integer :: i, n, nvars, sgn
        real(DP) :: sign_val

        nvars = 6

        ! Build identity DA vector x(i) = da_var(i)
        call x%init(nvars)
        do i = 1, nvars
            x%elements(i) = da_var(i)
        end do

        if (.not. allocated(history%entries)) return

        do i = 1, size(history%entries)
            n = abs(history%entries(i)) - 1
            sgn = history%entries(i) / abs(history%entries(i))
            sign_val = 0.5_DP * real(sgn, DP)

            ! x(n) = 0.5*sign + 0.5*da_var(n+1)
            tmp_da = 0.5_DP * da_var(n+1)
            x%elements(n+1) = sign_val + tmp_da

            ! BLOCKED: DA-to-DA substitution is not available in the current
            ! pod_dace_classes API. The following call from the spec cannot
            ! compile because AlgebraicVector%eval does not accept an
            ! AlgebraicVector argument. Available signatures:
            !   eval(integer, real(DP))  -> AlgebraicVector  (real substitution)
            !   eval(real(DP)(:))        -> real(DP), allocatable (full eval)
            ! Required: obj = obj.eval(x)  -- substitutes variables in each
            ! element of obj with the DA expressions in x (polynomial composition)
            !
            ! Once DA-to-DA substitution is added to the API, uncomment:
            !   obj = obj%eval(x)

            ! x(n) = da_var(n+1)  (restore)
            x%elements(n+1) = da_var(n+1)
        end do

        call x%destroy()
    end subroutine sh_replay

end module pod_ads_split_module
