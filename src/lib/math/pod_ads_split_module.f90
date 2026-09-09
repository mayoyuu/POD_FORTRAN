!> ADS (Automatic Domain Splitting) core module
!> Provides SplittingHistory, Patch, Manifold types and operations
!> Translated from C++ reference: ADS_cpp_Core_file/
module pod_ads_split_module
    use, intrinsic :: ieee_arithmetic, only: ieee_quiet_nan, ieee_value
    use pod_global, only: DP
    use pod_dace_classes, only: AlgebraicVector, CompiledDA, DA, da_var, da_estim_norm, &
        dace_max_variables, da_add, da_mul, vector_eval_da_vec_sub, &
        operator(+), operator(*), assignment(=)
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
    public :: mf_find_patch, mf_evaluate_point, mf_evaluate_points

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
        real(DP), allocatable :: c(:)
        real(DP), allocatable :: w(:)
        integer :: i, n, sgn, nvars
        real(DP) :: half_w
        nvars = dace_max_variables()
        allocate(c(nvars), w(nvars))
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
        real(DP), allocatable :: w(:)
        integer :: i, n, nvars
        nvars = dace_max_variables()
        allocate(w(nvars))
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
        real(DP), intent(in) :: pt(:)
        real(DP), allocatable :: c(:), w(:)
        integer :: i, nvars
        nvars = dace_max_variables()
        if (size(pt) /= nvars) then
            ok = .false.
            return
        end if
        c = sh_center(history)
        w = sh_width(history)
        ok = .true.
        do i = 1, nvars
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
        real(DP), intent(inout) :: pt(:)
        integer :: i, n
        if (size(pt) /= dace_max_variables()) return
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

        nvars = dace_max_variables()

        ! Build identity DA vector x(i) = da_var(i)
        call x%init(nvars)
        do i = 1, nvars
            x%elements(i) = da_var(i)
        end do

        if (.not. allocated(history%entries)) then
            call x%destroy()
            return
        end if

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
    ! Transfers DA storage ownership from da_vec into p without allocating or
    ! overwriting handles. After the call, da_vec has no allocated elements.
    ! To preserve the source, deep-copy into a temporary first.
    ! =========================================================================
    subroutine patch_init(p, da_vec, history)
        type(patch_type), intent(out) :: p
        type(AlgebraicVector), intent(inout) :: da_vec
        type(splitting_history_type), intent(in), optional :: history
        type(splitting_history_type) :: hist_copy
        integer :: n_components

        if (present(history)) hist_copy = history

        if (allocated(da_vec%elements)) then
            n_components = size(da_vec%elements)
        else
            n_components = 0
        end if
        ! The allocated container is authoritative, matching C++ vector::size().
        call p%da_vec%destroy()
        if (allocated(p%da_vec%h_list)) deallocate(p%da_vec%h_list)

        if (allocated(da_vec%elements)) then
            call move_alloc(da_vec%elements, p%da_vec%elements)
        end if
        if (allocated(da_vec%h_list)) then
            call move_alloc(da_vec%h_list, p%da_vec%h_list)
        end if

        p%da_vec%size = n_components
        da_vec%size = 0
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
        real(DP), intent(out) :: errors(:)
        integer :: i
        errors = 0.0_DP
        if (size(errors) < p%da_vec%size) then
            error stop 'patch_get_trunc_err: errors array is smaller than Patch output size'
        end if
        do i = 1, p%da_vec%size
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
        integer :: i, nvars
        dir = 1
        err_m = 0.0_DP
        nvars = dace_max_variables()
        do i = 1, nvars
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
    ! Matches C++ Patch::split: construct an identity DA map, replace the split
    ! coordinate with 0.5*x +/- 0.5, then compose every Patch output with it.
    ! =========================================================================
    subroutine patch_split(p, dir, left, right)
        type(patch_type), intent(inout) :: p
        integer, intent(in) :: dir
        type(patch_type), intent(out) :: left, right
        type(AlgebraicVector) :: map_vec
        type(DA) :: affine_da
        integer :: i, nvars

        if (.not. allocated(p%da_vec%elements)) &
            error stop 'patch_split: parent Patch has no DA vector'
        nvars = dace_max_variables()
        if (dir < 1 .or. dir > nvars) &
            error stop 'patch_split: split direction is outside the active DACE variables'

        call map_vec%init(nvars)
        do i = 1, nvars
            call map_vec%elements(i)%destroy()
            call map_vec%elements(i)%init_var(i)
        end do
        call affine_da%init()

        ! ---- Left half: v_dir -> 0.5*v_dir - 0.5 ----
        left%history = p%history
        call sh_push(left%history, -dir)
        call da_mul(map_vec%elements(dir), 0.5_DP, affine_da)
        call da_add(affine_da, -0.5_DP, map_vec%elements(dir))
        call vector_eval_da_vec_sub(p%da_vec, map_vec, left%da_vec)

        call map_vec%elements(dir)%destroy()
        call map_vec%elements(dir)%init_var(dir)

        ! ---- Right half: v_dir -> 0.5*v_dir + 0.5 ----
        right%history = p%history
        call sh_push(right%history, dir)
        call da_mul(map_vec%elements(dir), 0.5_DP, affine_da)
        call da_add(affine_da, 0.5_DP, map_vec%elements(dir))
        call vector_eval_da_vec_sub(p%da_vec, map_vec, right%da_vec)

        call affine_da%destroy()
        call map_vec%destroy()
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
        integer :: i, n
        if (.not. allocated(m%patches)) then
            allocate(m%patches(1))
            m%n_patches = 1
            m%patches(1) = p
        else
            n = m%n_patches
            allocate(tmp(n+1))
            tmp(1:n) = m%patches(1:n)
            tmp(n+1) = p
            do i = 1, n
                call patch_destroy(m%patches(i))
            end do
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
            call patch_destroy(m%patches(1))
            deallocate(m%patches)
            m%n_patches = 0
        else
            allocate(tmp(n-1))
            do i = 2, n
                tmp(i-1) = m%patches(i)
            end do
            do i = 1, n
                call patch_destroy(m%patches(i))
            end do
            call move_alloc(tmp, m%patches)
            m%n_patches = n - 1
        end if
    end subroutine mf_pop_front

    !> Return the first Patch containing a point in the global ADS unit box.
    !! Zero means the point is outside every Patch (or has the wrong dimension).
    subroutine mf_find_patch(m, point, patch_index)
        type(manifold_type), intent(in) :: m
        real(DP), intent(in) :: point(:)
        integer, intent(out) :: patch_index
        integer :: i

        patch_index = 0
        if (size(point) /= dace_max_variables()) return
        do i = 1, m%n_patches
            if (sh_contain(m%patches(i)%history, point)) then
                patch_index = i
                return
            end if
        end do
    end subroutine mf_find_patch

    !> Evaluate one global unit-box point using its owning Patch.
    !! value is set to NaN and found is false when no Patch contains the point.
    subroutine mf_evaluate_point(m, point, value, found, status)
        type(manifold_type), intent(in) :: m
        real(DP), intent(in) :: point(:)
        real(DP), intent(out) :: value(:)
        logical, intent(out) :: found
        integer, intent(out) :: status
        type(CompiledDA) :: compiled
        real(DP), allocatable :: local_point(:)
        integer :: patch_index

        status = 0
        found = .false.
        value = ieee_value(0.0_DP, ieee_quiet_nan)
        if (size(point) /= dace_max_variables()) then
            status = -1
            return
        end if
        if (m%n_patches < 1) return
        if (size(value) /= m%patches(1)%da_vec%size) then
            status = -2
            return
        end if

        call mf_find_patch(m, point, patch_index)
        if (patch_index == 0) return
        local_point = point
        call sh_map_point(m%patches(patch_index)%history, local_point)
        compiled = m%patches(patch_index)%da_vec%compile()
        call compiled%eval_into(local_point, value, status)
        call compiled%destroy()
        found = status == 0
    end subroutine mf_evaluate_point

    !> Evaluate sample columns across a split manifold in Patch-sized batches.
    !! Each Patch is compiled once, then all samples belonging to it are mapped
    !! to local coordinates and evaluated by one batch call. Outside points stay
    !! NaN with found=false; they are never silently replaced by a nominal state.
    subroutine mf_evaluate_points(m, points, values, found, status)
        type(manifold_type), intent(in) :: m
        real(DP), intent(in) :: points(:,:)
        real(DP), intent(out) :: values(:,:)
        logical, intent(out) :: found(:)
        integer, intent(out) :: status
        type(CompiledDA) :: compiled
        real(DP), allocatable :: local_points(:,:), local_values(:,:)
        integer, allocatable :: owners(:), sample_indices(:)
        integer :: i, j, n_assigned, n_samples, n_variables, n_outputs

        status = 0
        values = ieee_value(0.0_DP, ieee_quiet_nan)
        found = .false.
        n_variables = dace_max_variables()
        n_samples = size(points,2)
        if (size(points,1) /= n_variables .or. size(found) /= n_samples .or. &
            size(values,2) /= n_samples) then
            status = -1
            return
        end if
        if (m%n_patches < 1) then
            if (size(values,1) /= 0) status = -2
            return
        end if
        n_outputs = m%patches(1)%da_vec%size
        if (size(values,1) /= n_outputs) then
            status = -2
            return
        end if
        do i = 2, m%n_patches
            if (m%patches(i)%da_vec%size /= n_outputs) then
                status = -3
                return
            end if
        end do
        if (n_samples == 0) return

        allocate(owners(n_samples))
        owners = 0
        do j = 1, n_samples
            call mf_find_patch(m, points(:,j), owners(j))
        end do

        do i = 1, m%n_patches
            n_assigned = count(owners == i)
            if (n_assigned == 0) cycle
            allocate(local_points(n_variables,n_assigned))
            allocate(local_values(n_outputs,n_assigned))
            allocate(sample_indices(n_assigned))
            n_assigned = 0
            do j = 1, n_samples
                if (owners(j) /= i) cycle
                n_assigned = n_assigned + 1
                sample_indices(n_assigned) = j
                local_points(:,n_assigned) = points(:,j)
                call sh_map_point(m%patches(i)%history, local_points(:,n_assigned))
            end do
            compiled = m%patches(i)%da_vec%compile()
            call compiled%eval_batch_into(local_points, local_values, status)
            call compiled%destroy()
            if (status /= 0) return
            do j = 1, n_assigned
                values(:,sample_indices(j)) = local_values(:,j)
                found(sample_indices(j)) = .true.
            end do
            deallocate(local_points, local_values, sample_indices)
        end do
    end subroutine mf_evaluate_points

end module pod_ads_split_module
