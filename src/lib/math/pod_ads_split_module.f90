!> ADS (Automatic Domain Splitting) core module
!> Provides SplittingHistory, Patch, Manifold types and operations
!> Translated from C++ reference: ADS_cpp_Core_file/
module pod_ads_split_module
    use pod_global, only: DP
    use pod_dace_classes, only: AlgebraicVector, DA, da_var, da_estim_norm, &
        da_translate_variable, operator(+), operator(*), assignment(=)
    use iso_c_binding, only: c_int
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
    public :: patch_init, patch_destroy, patch_get_trunc_err, patch_get_split_dir, patch_split
    public :: mf_init, mf_destroy, mf_push, mf_pop_front

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

            ! DA-to-DA composition: substitute variables in obj with DA expressions from x
            obj = obj%eval(x)

            ! x(n) = da_var(n+1)  (restore)
            x%elements(n+1) = da_var(n+1)
        end do

        call x%destroy()
    end subroutine sh_replay

    ! =========================================================================
    ! Patch: patch_init
    !
    ! Moves DA handle ownership from da_vec into p. After the call, da_vec
    ! elements have handle = -1 (safe to destroy as a no-op).
    ! To preserve the source, deep-copy into a temporary first.
    ! =========================================================================
    subroutine patch_init(p, da_vec, history)
        type(patch_type), intent(out) :: p
        type(AlgebraicVector), intent(inout) :: da_vec
        type(splitting_history_type), intent(in), optional :: history
        type(splitting_history_type) :: hist_copy
        integer :: i

        if (present(history)) hist_copy = history

        call p%da_vec%destroy()
        call p%da_vec%init(6)
        do i = 1, 6
            p%da_vec%elements(i)%handle = da_vec%elements(i)%handle
            da_vec%elements(i)%handle = -1
        end do
        if (present(history)) then
            p%history = hist_copy
        end if
    end subroutine patch_init

    ! =========================================================================
    ! Patch: patch_destroy
    ! =========================================================================
    subroutine patch_destroy(p)
        type(patch_type), intent(inout) :: p
        call p%da_vec%destroy()
        if (allocated(p%history%entries)) deallocate(p%history%entries)
    end subroutine patch_destroy

    ! =========================================================================
    ! Patch: patch_get_trunc_err
    ! =========================================================================
    subroutine patch_get_trunc_err(p, order, errors)
        type(patch_type), intent(in) :: p
        integer, intent(in) :: order
        real(DP), intent(out) :: errors(6)
        integer :: i
        do i = 1, 6
            call da_estim_norm(p%da_vec%elements(i)%handle, 0, order, errors(i))
        end do
    end subroutine patch_get_trunc_err

    ! =========================================================================
    ! Patch: patch_get_split_dir
    ! =========================================================================
    integer function patch_get_split_dir(p, comp, order) result(dir)
        type(patch_type), intent(in) :: p
        integer, intent(in) :: comp, order
        real(DP) :: err_m, top_norm
        integer :: i
        dir = 1
        err_m = 0.0_DP
        do i = 1, 6
            call da_estim_norm(p%da_vec%elements(comp)%handle, i, order, top_norm)
            if (top_norm > err_m) then
                err_m = top_norm
                dir = i
            end if
        end do
    end function patch_get_split_dir

    ! =========================================================================
    ! Patch: patch_split
    !
    ! Splits patch p along direction dir using DACE translateVariable to
    ! perform the affine transformation v_dir -> 0.5*v_dir ± 0.5 on each
    ! component of the DA vector. This avoids the generic DA-to-DA eval path.
    !
    ! DIAGNOSTIC MODE: Set DIAGNOSTIC_SPLIT = .true. to bypass the affine
    ! transformation and simply deep-copy the original DA into both children.
    ! This isolates whether the corruption is in translateVariable or downstream.
    ! =========================================================================
    subroutine patch_split(p, dir, left, right)
        type(patch_type), intent(inout) :: p
        integer, intent(in) :: dir
        type(patch_type), intent(out) :: left, right
        type(AlgebraicVector) :: temp_vec
        type(splitting_history_type) :: saved_hist
        integer(c_int) :: new_handle
        integer :: i
        logical, parameter :: DIAGNOSTIC_SPLIT = .false.

        if (DIAGNOSTIC_SPLIT) then
            ! === DIAGNOSTIC: deep-copy original DA, no affine transform ===
            write(*,'(A,I0,A)') '[DIAG] patch_split called, dir=', dir, &
                ' -- BYPASSING translateVariable, deep-copying instead'

            left%history = p%history
            call sh_push(left%history, -dir)
            call temp_vec%init(6)
            do i = 1, 6
                temp_vec%elements(i) = p%da_vec%elements(i)
            end do
            saved_hist = left%history  ! break aliasing before patch_init
            call patch_init(left, temp_vec, saved_hist)
            call temp_vec%destroy()

            right%history = p%history
            call sh_push(right%history, dir)
            call temp_vec%init(6)
            do i = 1, 6
                temp_vec%elements(i) = p%da_vec%elements(i)
            end do
            saved_hist = right%history  ! break aliasing before patch_init
            call patch_init(right, temp_vec, saved_hist)
            call temp_vec%destroy()
            return
        end if

        ! ---- Left half: v_dir -> 0.5*v_dir - 0.5 ----
        left%history = p%history
        call sh_push(left%history, -dir)

        call temp_vec%init(6)
        do i = 1, 6
            call da_translate_variable( &
                p%da_vec%elements(i)%handle, dir, 0.5_DP, -0.5_DP, new_handle)
            temp_vec%elements(i)%handle = new_handle
        end do
        saved_hist = left%history  ! break aliasing: left%history undefined when patch_init intent(out) fires
        call patch_init(left, temp_vec, saved_hist)
        call temp_vec%destroy()

        ! ---- Right half: v_dir -> 0.5*v_dir + 0.5 ----
        right%history = p%history
        call sh_push(right%history, dir)

        call temp_vec%init(6)
        do i = 1, 6
            call da_translate_variable( &
                p%da_vec%elements(i)%handle, dir, 0.5_DP, 0.5_DP, new_handle)
            temp_vec%elements(i)%handle = new_handle
        end do
        saved_hist = right%history  ! break aliasing before patch_init
        call patch_init(right, temp_vec, saved_hist)
        call temp_vec%destroy()
    end subroutine patch_split

    ! =========================================================================
    ! Manifold: mf_init
    ! =========================================================================
    subroutine mf_init(m)
        type(manifold_type), intent(out) :: m
        allocate(m%patches(0))
        m%n_patches = 0
    end subroutine mf_init

    ! =========================================================================
    ! Manifold: mf_destroy
    ! =========================================================================
    subroutine mf_destroy(m)
        type(manifold_type), intent(inout) :: m
        integer :: i
        do i = 1, m%n_patches
            call patch_destroy(m%patches(i))
        end do
        if (allocated(m%patches)) deallocate(m%patches)
        m%n_patches = 0
    end subroutine mf_destroy

    ! =========================================================================
    ! Manifold: mf_push
    ! =========================================================================
    subroutine mf_push(m, p)
        type(manifold_type), intent(inout) :: m
        type(patch_type), intent(in) :: p
        type(patch_type), allocatable :: tmp(:)
        integer :: n
        if (.not. allocated(m%patches)) then
            allocate(m%patches(1))
            m%n_patches = 1
            m%patches(1) = p
        else
            n = m%n_patches
            allocate(tmp(n+1))
            tmp(1:n) = m%patches(1:n)
            tmp(n+1) = p
            call move_alloc(tmp, m%patches)
            m%n_patches = n + 1
        end if
    end subroutine mf_push

    ! =========================================================================
    ! Manifold: mf_pop_front — removes and returns the FIRST patch (FIFO queue)
    ! =========================================================================
    subroutine mf_pop_front(m, p)
        type(manifold_type), intent(inout) :: m
        type(patch_type), intent(out) :: p
        type(patch_type), allocatable :: tmp(:)
        integer :: n, i
        if (m%n_patches == 0) return
        p = m%patches(1)
        n = m%n_patches
        if (n == 1) then
            deallocate(m%patches)
            m%n_patches = 0
        else
            allocate(tmp(n-1))
            do i = 2, n
                tmp(i-1) = m%patches(i)
            end do
            call move_alloc(tmp, m%patches)
            m%n_patches = n - 1
        end if
    end subroutine mf_pop_front

end module pod_ads_split_module
