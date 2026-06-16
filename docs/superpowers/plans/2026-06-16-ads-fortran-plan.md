# ADS Uncertainty Propagation in POD_Fortran — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement Automatic Domain Splitting (ADS) uncertainty propagation in POD_Fortran, translating the C++ DACE-based reference into Fortran, with a unified CRTBP comparison test (DA/STT/ADS vs MC).

**Architecture:** Three-layer structure: core ADS splitting algorithm in `src/lib/math/pod_ads_split_module.f90` (procedural, SplittingHistory/Patch/Manifold types), `uq_ads_propagator` extending `uq_propagator_base` in `src/lib/uncertainty/propagation/pod_uq_prop_ads_module.f90`, CRTBP-specific wrapper in `pod_uq_crtbp_ads_module.f90`, plus unified comparison test. Requires adding `fdace_estim_norm` C wrapper to `dace_wrapper.cpp` + Fortran binding to `pod_dace_classes.f90`.

**Tech Stack:** Fortran 2008, DACE library (C++ with C bindings), fpm build system, CRTBP dynamics

---

### Task 0: Add `fdace_estim_norm` C Wrapper and Fortran Binding

**Files:**
- Modify: `external/dace_build/dace_wrapper.cpp`
- Modify: `src/lib/system/pod_dace_classes.f90`

- [ ] **Step 1: Add C wrapper function to dace_wrapper.cpp**

Insert after the `fdace_trunc` function (line ~243) in `external/dace_build/dace_wrapper.cpp`:

```cpp
// ==========================================
// 范数估计 (EstimNorm) — 用于 ADS 截断误差估计
// ==========================================
int fdace_estim_norm(int h_in, int var_idx, int max_order, double* norms) {
    std::vector<double> err;
    std::vector<double> result = da_registry[h_in].estimNorm(err, var_idx, 0, max_order + 1);
    int sz = static_cast<int>(result.size());
    for (int i = 0; i < sz; ++i) norms[i] = result[i];
    return sz;
}
```

- [ ] **Step 2: Rebuild libdace_wrapper.a**

```bash
cd external/dace_build/build && make
```

Expected: `libdace_wrapper.a` rebuilt without errors. Verify with:
```bash
nm external/dace_build/build/libdace_wrapper.a | grep estim_norm
```
Expected: symbol `fdace_estim_norm` present.

- [ ] **Step 3: Add Fortran binding to pod_dace_classes.f90**

Add C interface declaration in the `interface` block (near the other `c_fdace_*` declarations, after `c_fdace_trunc` around line 224):

```fortran
function c_fdace_estim_norm(h, var_idx, max_order, norms) bind(C, name="fdace_estim_norm") result(sz)
    import :: c_int, c_double
    integer(c_int), value :: h, var_idx, max_order
    real(c_double), intent(out) :: norms(*)
    integer(c_int) :: sz
end function c_fdace_estim_norm
```

Also add `c_fdace_estim_norm` to the `public` list at the top of the module (around line 17, after the existing `c_fdace_` imports — note: the interface functions are not individually public; just confirm the binding is accessible within the module).

Add a Fortran wrapper function in the `contains` section after the module's interface block:

```fortran
!> estimNorm wrapper: returns top-order norm for DA truncation error estimation
subroutine da_estim_norm(da_handle, var_idx, max_order, top_norm)
    integer(c_int), intent(in) :: da_handle
    integer, intent(in) :: var_idx, max_order
    real(DP), intent(out) :: top_norm
    real(c_double) :: norms(10)  ! max_order+1 entries, safe upper bound
    integer :: sz, i
    sz = c_fdace_estim_norm(da_handle, int(var_idx, c_int), int(max_order, c_int), norms)
    if (sz > 0) then
        top_norm = real(norms(sz), DP)
    else
        top_norm = 0.0_DP
    end if
end subroutine da_estim_norm
```

Expose `da_estim_norm` via `public :: da_estim_norm`.

- [ ] **Step 4: Commit**

```bash
git add external/dace_build/dace_wrapper.cpp external/dace_build/build/libdace_wrapper.a src/lib/system/pod_dace_classes.f90
git commit -m "feat: add fdace_estim_norm C wrapper and Fortran binding for ADS"
```

---

### Task 1: Core ADS Module — Types and SplittingHistory

**Files:**
- Create: `src/lib/math/pod_ads_split_module.f90`

- [ ] **Step 1: Create module skeleton with types and SplittingHistory operations**

```fortran
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
    public :: patch_init, patch_destroy, patch_get_trunc_err, patch_get_split_dir, patch_split
    public :: mf_push, mf_pop_front, mf_init, mf_destroy
    public :: ads_get_split_domain

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
    ! =========================================================================
    subroutine sh_replay(history, obj)
        type(splitting_history_type), intent(in) :: history
        type(AlgebraicVector), intent(inout) :: obj
        type(AlgebraicVector) :: x
        type(DA) :: tmp_da
        integer :: i, n, nvars, sgn, da_order
        real(DP) :: sign_val

        nvars = 6

        ! Build identity DA vector x(i) = da_var(i)
        call x%init(nvars)
        do i = 1, nvars
            x%elements(i) = da_var(i)
        end do

        if (.not. allocated(history%entries)) return
        da_order = DA::getMaxOrder()  ! save current max order (likely 2)

        do i = 1, size(history%entries)
            n = abs(history%entries(i)) - 1
            sgn = history%entries(i) / abs(history%entries(i))
            sign_val = 0.5_DP * real(sgn, DP)

            ! x(n) = 0.5*sign + 0.5*da_var(n+1)
            tmp_da = 0.5_DP * da_var(n+1)
            x%elements(n+1) = sign_val + tmp_da

            ! obj = obj.eval(x)
            obj = obj%eval(x)

            ! x(n) = da_var(n+1)  (restore)
            x%elements(n+1) = da_var(n+1)
        end do

        call x%destroy()
    end subroutine sh_replay

end module pod_ads_split_module
```

- [ ] **Step 2: Verify module compiles with fpm**

```bash
cd POD_Fortran && fpm build
```

Expected: module compiles without errors (unused function warnings OK at this stage).

- [ ] **Step 3: Commit**

```bash
git add src/lib/math/pod_ads_split_module.f90
git commit -m "feat: add pod_ads_split_module with SplittingHistory operations"
```

---

### Task 2: Patch Operations

**Files:**
- Modify: `src/lib/math/pod_ads_split_module.f90`

- [ ] **Step 1: Add patch_init, patch_destroy**

Append to the `contains` section of `pod_ads_split_module.f90`:

```fortran
    ! =========================================================================
    ! Patch: patch_init
    ! =========================================================================
    subroutine patch_init(p, da_vec, history)
        type(patch_type), intent(out) :: p
        type(AlgebraicVector), intent(in) :: da_vec
        type(splitting_history_type), intent(in), optional :: history
        integer :: i
        call p%da_vec%init(6)
        do i = 1, 6
            p%da_vec%elements(i) = da_vec%elements(i)
        end do
        if (present(history)) then
            p%history = history
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
```

- [ ] **Step 2: Add patch_get_trunc_err**

```fortran
    ! =========================================================================
    ! Patch: patch_get_trunc_err
    !   Returns truncation error (highest-order norm) for each of 6 components
    !   Corresponds to C++ Patch::getTruncationErrors
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
```

- [ ] **Step 3: Add patch_get_split_dir**

```fortran
    ! =========================================================================
    ! Patch: patch_get_split_dir
    !   For component `comp` with max error, find direction (1..6) with
    !   largest top-order contribution. Corresponds to C++ Patch::getSplittingDirection
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
```

- [ ] **Step 4: Add patch_split**

```fortran
    ! =========================================================================
    ! Patch: patch_split
    !   Split patch along direction `dir`, producing left and right sub-patches.
    !   Corresponds to C++ Patch::split
    ! =========================================================================
    subroutine patch_split(p, dir, left, right)
        type(patch_type), intent(inout) :: p
        integer, intent(in) :: dir
        type(patch_type), intent(out) :: left, right
        type(AlgebraicVector) :: obj, temp_result
        type(DA) :: tmp_da
        integer :: i

        ! Build identity DA vector
        call obj%init(6)
        do i = 1, 6
            obj%elements(i) = da_var(i)
        end do

        ! ---- Left half ----
        left%history = p%history
        call sh_push(left%history, -dir)

        ! obj(dir-1) = -0.5 + 0.5*da_var(dir)
        tmp_da = 0.5_DP * da_var(dir)
        obj%elements(dir) = -0.5_DP + tmp_da

        temp_result = p%da_vec%eval(obj)  ! substitute into patch's DA
        call patch_init(left, temp_result, left%history)
        call temp_result%destroy()

        ! Restore obj(dir-1)
        obj%elements(dir) = da_var(dir)

        ! ---- Right half ----
        call sh_pop(left%history)  ! undo the push we did above (clean up for right)
        right%history = p%history
        call sh_push(right%history, dir)

        ! obj(dir-1) = 0.5 + 0.5*da_var(dir)
        tmp_da = 0.5_DP * da_var(dir)
        obj%elements(dir) = 0.5_DP + tmp_da

        temp_result = p%da_vec%eval(obj)
        call patch_init(right, temp_result, right%history)
        call temp_result%destroy()

        call obj%destroy()
    end subroutine patch_split
```

- [ ] **Step 5: Verify module compiles**

```bash
cd POD_Fortran && fpm build
```

Expected: compiles without errors.

- [ ] **Step 6: Commit**

```bash
git add src/lib/math/pod_ads_split_module.f90
git commit -m "feat: add Patch operations (init/destroy/trunc_err/split_dir/split)"
```

---

### Task 3: Manifold Operations and ads_get_split_domain

**Files:**
- Modify: `src/lib/math/pod_ads_split_module.f90`

- [ ] **Step 1: Add mf_init, mf_destroy, mf_push, mf_pop_front**

Append to `pod_ads_split_module.f90`:

```fortran
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
    ! Manifold: mf_pop_front
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
```

- [ ] **Step 2: Add ads_get_split_domain — the core splitting algorithm**

```fortran
    ! =========================================================================
    ! ads_get_split_domain
    !   BFS domain splitting: repeatedly split patches where DA truncation
    !   error exceeds tolerance, until all patches satisfy error < err_toll
    !   or max splits per patch reached.
    !   Corresponds to C++ Manifold::getSplitDomain(func, errToll, nSplitMax, t)
    !
    !   da_integrator is called explicitly (CRTBP DA adaptive integrator)
    ! =========================================================================
    subroutine ads_get_split_domain(initial_da_vec, err_toll, n_split_max, &
                                     mu, t_end, abs_tol, rel_tol, dt_min, dt_max, &
                                     max_steps, da_order, result_manifold, verbose)
        type(AlgebraicVector), intent(in) :: initial_da_vec
        real(DP), intent(in) :: err_toll(6)
        integer, intent(in) :: n_split_max
        real(DP), intent(in) :: mu, t_end
        real(DP), intent(in) :: abs_tol, rel_tol, dt_min, dt_max
        integer, intent(in) :: max_steps, da_order
        type(manifold_type), intent(out) :: result_manifold
        logical, intent(in) :: verbose

        type(manifold_type) :: queue
        type(patch_type) :: p, f, left, right
        type(AlgebraicVector) :: f_da_vec
        type(splitting_history_type) :: empty_history
        real(DP) :: errors(6), rel_errors(6)
        integer :: pos(1), dir, iter_count, max_queue_size

        ! External CRTBP DA integrator (from pod_crtbp_integrator_module)
        ! declared in the calling module — call is explicit from crtbp_ads_module

        call mf_init(result_manifold)
        call mf_init(queue)

        ! Initial patch: da_vec from caller, empty history
        call patch_init(p, initial_da_vec)
        call mf_push(queue, p)

        iter_count = 0
        max_queue_size = 0

        do while (queue%n_patches > 0)

            iter_count = iter_count + 1
            max_queue_size = max(max_queue_size, queue%n_patches)

            call mf_pop_front(queue, p)

            ! --- Evaluate flow: integrate p%da_vec via CRTBP DA integrator ---
            call f_da_vec%init(6)
            call da_adaptive_integrate_crtbp( &
                p%da_vec, 0.0_DP, t_end, &
                rel_tol, abs_tol, dt_min, dt_max, max_steps, &
                f_da_vec)

            ! f = (f_da_vec, p%history)
            call patch_init(f, f_da_vec, p%history)
            call f_da_vec%destroy()

            ! --- Truncation error check ---
            call patch_get_trunc_err(f, da_order, errors)
            rel_errors = max(0.0_DP, errors - err_toll)

            if (maxval(rel_errors) <= 0.0_DP .or. sh_count(p%history, 0) >= n_split_max) then
                ! Accept this patch
                call mf_push(result_manifold, f)
                call patch_destroy(p)
            else
                ! Split
                pos = maxloc(rel_errors)
                dir = patch_get_split_dir(f, pos(1), da_order)
                call patch_split(p, dir, left, right)
                call mf_push(queue, left)
                call mf_push(queue, right)
                call patch_destroy(p)
                call patch_destroy(f)
            end if
        end do

        call mf_destroy(queue)

        if (verbose) then
            write(*,'(A,I0,A,I0,A,I0)') '[ADS] splitting done: ', &
                result_manifold%n_patches, ' patches, ', &
                iter_count, ' iterations, max queue size ', max_queue_size
        end if
    end subroutine ads_get_split_domain
```

**Note:** The explicit call to `da_adaptive_integrate_crtbp` requires this module to `use pod_crtbp_integrator_module`. This creates a circular-ish dependency (math → integrator). To avoid that, `ads_get_split_domain` will be split: the core queue/split logic stays in `pod_ads_split_module`, but the RHS-integration call is injected from `pod_uq_crtbp_ads_module`. The `ads_get_split_domain` subroutine above is placed in `pod_uq_crtbp_ads_module` instead, and the pure splitting logic (without integration call) stays in the math module.

- [ ] **Step 3: Verify build**

```bash
cd POD_Fortran && fpm build
```

Expected: compiles (may have unused warnings for `ads_get_split_domain` if the `use pod_crtbp_integrator_module` is deferred).

- [ ] **Step 4: Commit**

```bash
git add src/lib/math/pod_ads_split_module.f90
git commit -m "feat: add Manifold operations and core ADS split domain skeleton"
```

---

### Task 4: ADS Propagator Type

**Files:**
- Create: `src/lib/uncertainty/propagation/pod_uq_prop_ads_module.f90`

- [ ] **Step 1: Create the propagator module**

```fortran
!> ADS 误差传播模块
!>
!> 实现基于 Automatic Domain Splitting 的非线性不确定性传播。
!> 继承 uq_propagator_base，提供 propagate 接口。
!>
!> 架构:
!>   - uq_ads_propagator 类型 (继承 uq_propagator_base)
!>   - DA 初始域构建 → 域分裂 → 采样点映射求值 → 矩统计
!>
!> 依赖:
!>   - pod_global, pod_uq_base_module, pod_uq_state_module
!>   - pod_ads_split_module, pod_crtbp_integrator_module
!>   - pod_dace_classes
module pod_uq_prop_ads_module
    use pod_global, only: DP, MAX_STRING_LEN
    use pod_uq_base_module, only: uq_propagator_base
    use pod_uq_state_module, only: uq_state_type
    use pod_integrator_module, only: METHOD_RKF45, METHOD_RKF78
    use pod_ads_split_module, only: splitting_history_type, patch_type, &
        manifold_type, sh_push, sh_pop, sh_count, sh_center, sh_width, &
        sh_contain, sh_map_point, patch_init, patch_destroy, &
        patch_get_trunc_err, patch_get_split_dir, patch_split, &
        mf_init, mf_destroy, mf_push, mf_pop_front, ads_get_split_domain
    use pod_crtbp_integrator_module, only: da_adaptive_integrate_crtbp
    use pod_dace_classes, only: AlgebraicVector, DA, da_var, &
        dace_initialize, dace_push_to, dace_pop_to
    implicit none
    private

    public :: uq_ads_propagator

    type, extends(uq_propagator_base), public :: uq_ads_propagator
        integer  :: da_order = 2
        integer  :: n_split_max = 12
        real(DP) :: err_toll(6) = 1.0d-4
        real(DP) :: mu = 0.012153614091892_DP
        real(DP) :: ads_abs_tol = 1.0d-12
        real(DP) :: ads_rel_tol = 1.0d-12
        real(DP) :: ads_dt_min = 1.0d-6
        real(DP) :: ads_dt_max = 3600.0_DP
        integer  :: ads_max_steps = 100000
        type(manifold_type) :: domain
        integer  :: n_patches = 0
    contains
        procedure :: propagate       => ads_propagate
        procedure :: get_method_name => ads_get_method_name
        procedure :: set_ads_order
        procedure :: set_ads_tolerances
        procedure :: set_ads_mu
    end type uq_ads_propagator

contains

    function ads_get_method_name(this) result(name)
        class(uq_ads_propagator), intent(in) :: this
        character(len=MAX_STRING_LEN) :: name
        write(name, '(A,I0)') 'ADS-', this%da_order
    end function ads_get_method_name

    subroutine set_ads_order(this, order)
        class(uq_ads_propagator), intent(inout) :: this
        integer, intent(in) :: order
        this%da_order = order
    end subroutine set_ads_order

    subroutine set_ads_tolerances(this, abs_tol, rel_tol, dt_min, dt_max, &
                                   max_steps, err_toll, n_split_max)
        class(uq_ads_propagator), intent(inout) :: this
        real(DP), optional, intent(in) :: abs_tol, rel_tol, dt_min, dt_max
        integer, optional, intent(in) :: max_steps, n_split_max
        real(DP), optional, intent(in) :: err_toll(6)
        if (present(abs_tol)) this%ads_abs_tol = abs_tol
        if (present(rel_tol)) this%ads_rel_tol = rel_tol
        if (present(dt_min)) this%ads_dt_min = dt_min
        if (present(dt_max)) this%ads_dt_max = dt_max
        if (present(max_steps)) this%ads_max_steps = max_steps
        if (present(err_toll)) this%err_toll = err_toll
        if (present(n_split_max)) this%n_split_max = n_split_max
    end subroutine set_ads_tolerances

    subroutine set_ads_mu(this, mu)
        class(uq_ads_propagator), intent(inout) :: this
        real(DP), intent(in) :: mu
        this%mu = mu
    end subroutine set_ads_mu

    ! =========================================================================
    ! Main propagate (delegates to crtbp_ads_propagate_deviates)
    ! =========================================================================
    subroutine ads_propagate(this, t_start, t_end, input_state, output_state)
        class(uq_ads_propagator), intent(inout) :: this
        real(DP), intent(in) :: t_start, t_end
        type(uq_state_type), intent(in) :: input_state
        type(uq_state_type), intent(inout) :: output_state

        real(DP), allocatable :: deviates(:,:), final_samples(:,:)
        real(DP) :: final_mean(6), final_cov(6,6)
        integer :: n_particles, dim, i, j

        dim = 6
        n_particles = 100000

        ! Generate deviates from initial covariance
        if (.not. allocated(input_state%mean)) then
            write(*,*) '[ADS] ERROR: input_state%mean not allocated'
            return
        end if

        allocate(deviates(dim, n_particles))
        ! Simple Monte Carlo sampling from covariance:
        ! For now, use the existing random module
        ! (crtbp_ads_propagate_from_cov will handle this)

        ! Delegate to CRTBP-specific implementation
        call crtbp_ads_propagate_deviates_from_cov( &
            input_state%mean(1:6), input_state%cov(1:6,1:6), &
            this%mu, t_end - t_start, this%da_order, &
            this%n_split_max, this%err_toll, &
            this%ads_rel_tol, this%ads_abs_tol, &
            this%ads_dt_min, this%ads_dt_max, this%ads_max_steps, &
            final_samples, final_mean, final_cov, this%n_patches, this%verbose)

        ! Fill output
        call output_state%deallocate_memory()
        if (allocated(output_state%mean)) deallocate(output_state%mean)
        if (allocated(output_state%cov)) deallocate(output_state%cov)
        allocate(output_state%mean(6), output_state%cov(6,6))
        output_state%mean = final_mean
        output_state%cov = final_cov

        if (allocated(deviates)) deallocate(deviates)
        if (allocated(final_samples)) deallocate(final_samples)
    end subroutine ads_propagate

end module pod_uq_prop_ads_module
```

- [ ] **Step 2: Verify compiles (note: will have unresolved ref until Task 5)**

```bash
cd POD_Fortran && fpm build
```

Expected: may fail on `crtbp_ads_propagate_deviates_from_cov` — this is defined in Task 5.

- [ ] **Step 3: Commit**

```bash
git add src/lib/uncertainty/propagation/pod_uq_prop_ads_module.f90
git commit -m "feat: add uq_ads_propagator type extending uq_propagator_base"
```

---

### Task 5: CRTBP ADS Module — Concrete Propagation

**Files:**
- Create: `src/lib/uncertainty/propagation/pod_uq_crtbp_ads_module.f90`

- [ ] **Step 1: Create crtbp_ads_propagate_deviates**

```fortran
!> CRTBP ADS 误差传播模块
!>
!> 提供 CRTPB 动力学下的 ADS 传播具体实现。
!> 包含 ads_get_split_domain (含 CRTBP DA 积分器显式调用) 和
!> crtbp_ads_propagate_deviates (批量偏离量映射)。
!>
!> 架构:
!>   - ads_get_split_domain: DA 域分裂 (BFS)
!>   - crtbp_ads_propagate_deviates: 偏离量 → 分裂域求值 → 矩统计
!>
!> 依赖:
!>   - pod_ads_split_module (core types)
!>   - pod_crtbp_integrator_module (DA adaptive integrator)
!>   - pod_dace_classes (DA operations)
module pod_uq_crtbp_ads_module
    use pod_global, only: DP
    use pod_ads_split_module, only: splitting_history_type, patch_type, &
        manifold_type, sh_push, sh_pop, sh_count, sh_center, sh_width, &
        sh_contain, sh_map_point, patch_init, patch_destroy, &
        patch_get_trunc_err, patch_get_split_dir, patch_split, &
        mf_init, mf_destroy, mf_push, mf_pop_front
    use pod_crtbp_integrator_module, only: da_adaptive_integrate_crtbp
    use pod_dace_classes, only: AlgebraicVector, DA, da_var, &
        dace_initialize, dace_push_to, dace_pop_to, CompiledDA
    use pod_random_module, only: generate_multivariate_normal, init_random_seed
    implicit none
    private

    public :: crtbp_ads_propagate, crtbp_ads_propagate_deviates
    public :: ads_get_split_domain_crtbp

contains

    ! =========================================================================
    ! ads_get_split_domain_crtbp — CRTBP-specific domain splitting
    !
    !   Same as ads_get_split_domain, but with explicit CRTBP DA integrator
    !   call baked in. This is the Fortran equivalent of the C++
    !   Manifold::getSplitDomain(func, errToll, nSplitMax, t).
    ! =========================================================================
    subroutine ads_get_split_domain_crtbp(initial_da_vec, err_toll, n_split_max, &
                                           mu, t_end, abs_tol, rel_tol, dt_min, dt_max, &
                                           max_steps, da_order, result_manifold, verbose)
        type(AlgebraicVector), intent(in) :: initial_da_vec
        real(DP), intent(in) :: err_toll(6)
        integer, intent(in) :: n_split_max
        real(DP), intent(in) :: mu, t_end
        real(DP), intent(in) :: abs_tol, rel_tol, dt_min, dt_max
        integer, intent(in) :: max_steps, da_order
        type(manifold_type), intent(out) :: result_manifold
        logical, intent(in) :: verbose

        type(manifold_type) :: queue
        type(patch_type) :: p, f, left, right
        type(AlgebraicVector) :: f_da_vec
        type(splitting_history_type) :: empty_history
        real(DP) :: errors(6), rel_errors(6)
        integer :: pos(1), dir, iter_count, max_queue_size

        call mf_init(result_manifold)
        call mf_init(queue)

        ! Initial patch from caller's DA vector, empty history
        call patch_init(p, initial_da_vec)
        call mf_push(queue, p)

        iter_count = 0
        max_queue_size = 0

        do while (queue%n_patches > 0)
            iter_count = iter_count + 1
            max_queue_size = max(max_queue_size, queue%n_patches)

            call mf_pop_front(queue, p)

            ! --- CRTBP DA integration ---
            call f_da_vec%init(6)
            call da_adaptive_integrate_crtbp( &
                p%da_vec, 0.0_DP, t_end, &
                rel_tol, abs_tol, dt_min, dt_max, max_steps, &
                f_da_vec)

            ! f = Patch(f_da_vec, p%history)
            call patch_init(f, f_da_vec, p%history)
            call f_da_vec%destroy()

            ! --- Truncation error ---
            call patch_get_trunc_err(f, da_order, errors)
            rel_errors = max(0.0_DP, errors - err_toll)

            if (maxval(rel_errors) <= 0.0_DP .or. &
                sh_count(p%history, 0) >= n_split_max) then
                call mf_push(result_manifold, f)
                call patch_destroy(p)
            else
                pos = maxloc(rel_errors)
                dir = patch_get_split_dir(f, pos(1), da_order)
                call patch_split(p, dir, left, right)
                call mf_push(queue, left)
                call mf_push(queue, right)
                call patch_destroy(p)
                call patch_destroy(f)
            end if
        end do

        call mf_destroy(queue)

        if (verbose) then
            write(*,'(A,I0,A,I0,A,I0)') '[ADS CRTBP] splitting done: ', &
                result_manifold%n_patches, ' patches, ', &
                iter_count, ' iterations, max queue size ', max_queue_size
        end if
    end subroutine ads_get_split_domain_crtbp

    ! =========================================================================
    ! crtbp_ads_propagate_deviates — propagate given deviates via ADS
    !
    !   1. Build DA initial state with da_var
    !   2. Call ads_get_split_domain_crtbp
    !   3. Map each deviate through split manifold
    !   4. Compute moments
    ! =========================================================================
    subroutine crtbp_ads_propagate_deviates(nominal_state, deviates, mu, t_end, &
            da_order, n_split_max, err_toll, rel_tol, abs_tol, dt_min, dt_max, max_steps, &
            final_samples, final_mean, final_cov, n_patches_out, verbose)
        real(DP), intent(in) :: nominal_state(6)
        real(DP), intent(in) :: deviates(:,:)
        real(DP), intent(in) :: mu, t_end
        integer,  intent(in) :: da_order, n_split_max, max_steps
        real(DP), intent(in) :: err_toll(6)
        real(DP), intent(in) :: rel_tol, abs_tol, dt_min, dt_max
        real(DP), allocatable, intent(out) :: final_samples(:,:)
        real(DP), intent(out) :: final_mean(6), final_cov(6,6)
        integer,  intent(out) :: n_patches_out
        logical, intent(in) :: verbose

        type(AlgebraicVector) :: state_da_0
        type(manifold_type) :: manifold
        real(DP) :: pt_local(6)
        real(DP), allocatable :: eval_res(:)
        integer :: dim, n_dev, i, j, k
        logical :: found

        dim = 6
        n_dev = size(deviates, 2)

        ! 1. Init DACE
        call dace_initialize(da_order, dim)

        ! 2. Build initial DA state: x0(i) = nominal(i) + da_var(i)
        call state_da_0%init(dim)
        do i = 1, dim
            state_da_0%elements(i) = nominal_state(i) + da_var(i)
        end do

        ! 3. Domain splitting
        if (verbose) write(*,'(A)') '[ADS CRTBP] Starting domain splitting...'
        call ads_get_split_domain_crtbp(state_da_0, err_toll, n_split_max, &
            mu, t_end, abs_tol, rel_tol, dt_min, dt_max, max_steps, &
            da_order, manifold, verbose)
        n_patches_out = manifold%n_patches

        ! 4. Map deviates
        allocate(final_samples(dim, n_dev))
        final_samples = 0.0_DP

        do k = 1, n_dev
            found = .false.
            do i = 1, manifold%n_patches
                if (sh_contain(manifold%patches(i)%history, deviates(:, k))) then
                    pt_local = deviates(:, k)
                    call sh_map_point(manifold%patches(i)%history, pt_local)
                    eval_res = manifold%patches(i)%da_vec%eval(pt_local)
                    final_samples(:, k) = eval_res
                    found = .true.
                    deallocate(eval_res)
                    exit
                end if
            end do
            if (.not. found) then
                write(*,'(A,I0,A)') '[ADS CRTBP] WARNING: deviate ', k, ' not in any patch!'
                final_samples(:, k) = nominal_state  ! fallback
            end if
        end do

        ! 5. Compute moments
        final_mean = 0.0_DP
        do k = 1, n_dev
            final_mean = final_mean + final_samples(:, k)
        end do
        final_mean = final_mean / real(n_dev, DP)

        final_cov = 0.0_DP
        do k = 1, n_dev
            do j = 1, dim
                final_cov(:, j) = final_cov(:, j) + &
                    (final_samples(:, k) - final_mean) * (final_samples(j, k) - final_mean(j))
            end do
        end do
        final_cov = final_cov / real(n_dev - 1, DP)

        ! 6. Cleanup
        call state_da_0%destroy()
        call mf_destroy(manifold)

        if (verbose) write(*,'(A)') '[ADS CRTBP] Propagation complete.'
    end subroutine crtbp_ads_propagate_deviates

    ! =========================================================================
    ! crtbp_ads_propagate — full cov → samples → propagate pipeline
    ! =========================================================================
    subroutine crtbp_ads_propagate(nominal_state, cov, mu, t_end, &
            da_order, n_split_max, err_toll, rel_tol, abs_tol, dt_min, dt_max, max_steps, &
            n_particles, final_samples, final_mean, final_cov, n_patches_out, verbose)
        real(DP), intent(in) :: nominal_state(6), cov(6,6), mu, t_end
        integer,  intent(in) :: da_order, n_split_max, max_steps, n_particles
        real(DP), intent(in) :: err_toll(6)
        real(DP), intent(in) :: rel_tol, abs_tol, dt_min, dt_max
        real(DP), allocatable, intent(out) :: final_samples(:,:)
        real(DP), intent(out) :: final_mean(6), final_cov(6,6)
        integer,  intent(out) :: n_patches_out
        logical, intent(in) :: verbose

        real(DP), allocatable :: particles(:,:), deviates(:,:)
        integer :: i, j, dim

        dim = 6

        ! Generate particles
        allocate(particles(dim, n_particles))
        call init_random_seed(.true.)
        call generate_multivariate_normal(nominal_state, cov, particles)

        ! Convert to deviates
        allocate(deviates(dim, n_particles))
        do i = 1, n_particles
            deviates(:, i) = particles(:, i) - nominal_state
        end do

        call crtbp_ads_propagate_deviates(nominal_state, deviates, mu, t_end, &
            da_order, n_split_max, err_toll, rel_tol, abs_tol, dt_min, dt_max, max_steps, &
            final_samples, final_mean, final_cov, n_patches_out, verbose)

        deallocate(particles, deviates)
    end subroutine crtbp_ads_propagate

end module pod_uq_crtbp_ads_module
```

- [ ] **Step 2: Update pod_uq_prop_ads_module to use crtbp_ads_propagate_deviates**

The `ads_propagate` subroutine in `pod_uq_prop_ads_module.f90` needs to call into `pod_uq_crtbp_ads_module`. Update it to use `crtbp_ads_propagate`:

Replace the `ads_propagate` body with a direct call:

```fortran
    subroutine ads_propagate(this, t_start, t_end, input_state, output_state)
        class(uq_ads_propagator), intent(inout) :: this
        real(DP), intent(in) :: t_start, t_end
        type(uq_state_type), intent(in) :: input_state
        type(uq_state_type), intent(inout) :: output_state

        real(DP), allocatable :: final_samples(:,:)
        real(DP) :: final_mean(6), final_cov(6,6)
        integer :: n_particles, n_patches_out

        n_particles = 100000

        if (.not. allocated(input_state%mean)) then
            write(*,*) '[ADS] ERROR: input_state%mean not allocated'
            return
        end if

        call crtbp_ads_propagate( &
            input_state%mean(1:6), input_state%cov(1:6,1:6), &
            this%mu, t_end - t_start, this%da_order, &
            this%n_split_max, this%err_toll, &
            this%ads_rel_tol, this%ads_abs_tol, &
            this%ads_dt_min, this%ads_dt_max, this%ads_max_steps, &
            n_particles, &
            final_samples, final_mean, final_cov, n_patches_out, this%verbose)

        this%n_patches = n_patches_out

        call output_state%deallocate_memory()
        if (allocated(output_state%mean)) deallocate(output_state%mean)
        if (allocated(output_state%cov)) deallocate(output_state%cov)
        allocate(output_state%mean(6), output_state%cov(6,6))
        output_state%mean = final_mean
        output_state%cov = final_cov

        if (allocated(final_samples)) deallocate(final_samples)
    end subroutine ads_propagate
```

And add `use pod_uq_crtbp_ads_module, only: crtbp_ads_propagate` at the top.

- [ ] **Step 3: Verify build**

```bash
cd POD_Fortran && fpm build
```

Expected: builds successfully.

- [ ] **Step 4: Commit**

```bash
git add src/lib/uncertainty/propagation/pod_uq_crtbp_ads_module.f90 src/lib/uncertainty/propagation/pod_uq_prop_ads_module.f90
git commit -m "feat: add CRTBP ADS propagation module with deviates and cov modes"
```

---

### Task 6: Unified Comparison Test

**Files:**
- Create: `test/test_uq_crtbp_comparison.f90`
- Modify: `fpm.toml` (add test target)

- [ ] **Step 1: Create test program**

```fortran
!> CRTBP 不确定性传播综合对比测试
!>
!> 对比四种方法的传播结果和耗时:
!>   MC (Monte Carlo) — 数值积分参考解
!>   DA (Differential Algebra) — 单 DA 扩展 + 采样求值
!>   STT (State Transition Tensor) — 增广积分 + STT 多项式映射
!>   ADS (Automatic Domain Splitting) — 域分裂 + 分段 DA 求值
!>
!> 基于 test_stt_crtbp.f90 的 Test 6 模式扩展。
program test_uq_crtbp_comparison
    use pod_global, only: DP
    use pod_dace_classes, only: dace_initialize
    use pod_uq_crtbp_da_module, only: crtbp_da_propagate_deviates
    use pod_uq_crtbp_mc_module, only: crtbp_mc_propagate_deviates
    use pod_uq_prop_stt_module, only: uq_stt_propagator, stt_propagate_deviates
    use pod_uq_crtbp_ads_module, only: crtbp_ads_propagate_deviates
    implicit none

    integer :: passed, failed
    real(DP) :: mu, dro_x0(6), dro_period
    real(DP), allocatable :: deviates(:,:)
    real(DP), allocatable :: mc_samples(:,:), da_samples(:,:)
    real(DP), allocatable :: stt_samples(:,:), ads_samples(:,:)
    real(DP), allocatable :: mc_dev_subset(:,:)
    real(DP) :: mc_mean(6), mc_cov(6,6)
    real(DP) :: da_mean(6), da_cov(6,6), da_ref(6)
    real(DP) :: stt_mean(6), stt_cov(6,6)
    real(DP) :: ads_mean(6), ads_cov(6,6)
    type(uq_stt_propagator) :: stt_prop
    real(DP) :: mc_time, da_time, stt_time, ads_time
    integer(8) :: t1, t2, t_rate
    integer :: test_order, n_dev, n_mc, unit, io, k, i, j
    integer :: ads_n_patches
    real(DP) :: err_toll(6)

    passed = 0
    failed = 0

    write(*,*) '========================================'
    write(*,*) '  CRTBP Uncertainty Propagation'
    write(*,*) '  MC vs DA vs STT vs ADS Comparison'
    write(*,*) '========================================'

    ! Parameters (same as test_stt_crtbp Test 6)
    test_order = 2
    n_dev = 100000
    n_mc = 100000
    mu = 0.012153614091892_DP
    dro_x0 = [1.1309265107780351_DP, 0.0_DP, 0.0_DP, &
              0.0_DP, -0.46540743845849059_DP, 0.0_DP]
    dro_period = 2.3017923284002024_DP * 1.5_DP

    write(*,'(A,I0,A,I0,A,I0)') ' Config: order=', test_order, &
        ', N_dev=', n_dev, ', N_mc=', n_mc

    ! ---- 1. Load deviates ----
    write(*,*) 'Loading deviates from rand_list200km_0.7ms.txt...'
    allocate(deviates(6, n_dev))
    open(newunit=unit, file='../rand_list200km_0.7ms.txt', status='old', action='read')
    do k = 1, n_dev
        read(unit, *, iostat=io) deviates(1, k), deviates(2, k), deviates(3, k), &
                                  deviates(4, k), deviates(5, k), deviates(6, k)
        if (io /= 0) then
            write(*,*) 'FAIL: error reading line', k
            failed = failed + 1
            close(unit)
            go to 99
        end if
    end do
    close(unit)
    write(*,*) 'Loaded', n_dev, 'deviates'

    ! ---- 2. MC reference (n_mc subset) ----
    write(*,*) 'Running MC...'
    allocate(mc_dev_subset(6, n_mc))
    mc_dev_subset = deviates(:, 1:n_mc)
    call system_clock(t1, t_rate)
    call crtbp_mc_propagate_deviates(dro_x0, mc_dev_subset, mu, dro_period, &
        1.0e-14_DP, 1.0e-14_DP, 1.0e-10_DP, 0.01_DP, 100000, &
        mc_samples, mc_mean, mc_cov, .false.)
    call system_clock(t2, t_rate)
    mc_time = real(t2 - t1, DP) / real(t_rate, DP)

    ! ---- 3. DA ----
    write(*,*) 'Running DA...'
    call dace_initialize(test_order, 6)
    call system_clock(t1, t_rate)
    call crtbp_da_propagate_deviates(dro_x0, deviates, mu, dro_period, &
        test_order, 1.0e-14_DP, 1.0e-14_DP, 1.0e-10_DP, 0.01_DP, 100000, &
        da_samples, da_mean, da_cov, da_ref, .false.)
    call system_clock(t2, t_rate)
    da_time = real(t2 - t1, DP) / real(t_rate, DP)

    ! ---- 4. STT ----
    write(*,*) 'Running STT...'
    call stt_prop%set_stt_order(test_order)
    call stt_prop%set_stt_mu(mu)
    call stt_prop%set_stt_tolerances(abs_tol=1.0e-14_DP, rel_tol=1.0e-14_DP, &
        dt_min=1.0e-10_DP, dt_max=0.01_DP, max_steps=100000)
    call stt_prop%set_verbosity(.false.)
    call system_clock(t1, t_rate)
    call stt_propagate_deviates(stt_prop, 0.0_DP, dro_period, &
        dro_x0, deviates, stt_samples, stt_mean, stt_cov)
    call system_clock(t2, t_rate)
    stt_time = real(t2 - t1, DP) / real(t_rate, DP)

    ! ---- 5. ADS ----
    write(*,*) 'Running ADS...'
    err_toll = 1.0d-4
    call system_clock(t1, t_rate)
    call crtbp_ads_propagate_deviates(dro_x0, deviates, mu, dro_period, &
        test_order, 12, err_toll, &
        1.0e-14_DP, 1.0e-14_DP, 1.0e-10_DP, 0.01_DP, 100000, &
        ads_samples, ads_mean, ads_cov, ads_n_patches, .true.)
    call system_clock(t2, t_rate)
    ads_time = real(t2 - t1, DP) / real(t_rate, DP)

    ! ---- 6. Results comparison ----
    write(*,*) '========================================'
    write(*,*) '  Results'
    write(*,*) '========================================'

    ! Mean comparison
    write(*,'(A,6ES14.6)') ' MC  mean:', mc_mean
    write(*,'(A,6ES14.6)') ' DA  mean:', da_mean
    write(*,'(A,6ES14.6)') ' STT mean:', stt_mean
    write(*,'(A,6ES14.6)') ' ADS mean:', ads_mean

    write(*,'(A,ES14.6)') ' |DA - MC|_max  :', maxval(abs(da_mean - mc_mean))
    write(*,'(A,ES14.6)') ' |STT - MC|_max :', maxval(abs(stt_mean - mc_mean))
    write(*,'(A,ES14.6)') ' |ADS - MC|_max :', maxval(abs(ads_mean - mc_mean))

    ! Covariance comparison
    call compare_cov('DA  vs MC', da_cov, mc_cov, 0.10_DP, passed, failed)
    call compare_cov('STT vs MC', stt_cov, mc_cov, 0.10_DP, passed, failed)
    call compare_cov('ADS vs MC', ads_cov, mc_cov, 0.15_DP, passed, failed)  ! looser for ADS

    ! ---- 7. Timing ----
    write(*,*) '========================================'
    write(*,*) '  Timing'
    write(*,*) '========================================'
    write(*,'(A,F10.3,A)') ' MC  : ', mc_time,  ' s'
    write(*,'(A,F10.3,A)') ' DA  : ', da_time,  ' s'
    write(*,'(A,F10.3,A)') ' STT : ', stt_time, ' s'
    write(*,'(A,F10.3,A,I0,A)') ' ADS : ', ads_time, ' s (', ads_n_patches, ' patches)'
    write(*,'(A,F10.1)') ' Speedup vs MC (DA/MC):  ', mc_time / max(da_time, 1.0e-6_DP)
    write(*,'(A,F10.1)') ' Speedup vs MC (STT/MC): ', mc_time / max(stt_time, 1.0e-6_DP)
    write(*,'(A,F10.1)') ' Speedup vs MC (ADS/MC): ', mc_time / max(ads_time, 1.0e-6_DP)

    ! ---- 8. Output samples ----
    call write_samples('veri_crtbp_mc.txt', mc_samples, n_mc)
    call write_samples('o2_veri_crtbp_da.txt', da_samples, n_dev)
    call write_samples('o2_veri_crtbp_stt.txt', stt_samples, n_dev)
    call write_samples('o2_veri_crtbp_ads.txt', ads_samples, n_dev)

    ! ---- Cleanup ----
    deallocate(deviates, mc_dev_subset)
    if (allocated(mc_samples)) deallocate(mc_samples)
    if (allocated(da_samples)) deallocate(da_samples)
    if (allocated(stt_samples)) deallocate(stt_samples)
    if (allocated(ads_samples)) deallocate(ads_samples)

99  continue
    write(*,*) '========================================'
    write(*,'(A,I0,A,I0,A)') ' Results: ', passed, ' passed, ', failed, ' failed'
    write(*,*) '========================================'
    if (failed > 0) stop 1

contains

    subroutine compare_cov(label, cov1, cov2, threshold, passed, failed)
        character(len=*), intent(in) :: label
        real(DP), intent(in) :: cov1(6,6), cov2(6,6)
        real(DP), intent(in) :: threshold
        integer, intent(inout) :: passed, failed
        real(DP) :: fn1, fn2, diff_fn, rel_err
        integer :: r, c
        fn1 = 0.0_DP; fn2 = 0.0_DP; diff_fn = 0.0_DP
        do r = 1, 6
            do c = 1, 6
                fn1 = fn1 + cov1(r,c)**2
                fn2 = fn2 + cov2(r,c)**2
                diff_fn = diff_fn + (cov1(r,c) - cov2(r,c))**2
            end do
        end do
        rel_err = sqrt(diff_fn) / max(sqrt(fn2), 1.0e-15_DP)
        write(*,'(A,A,ES14.6)') ' Cov ', label, ' rel err:', rel_err
        if (rel_err < threshold) then
            write(*,'(A,A,A)') '  PASS: ', label, ' covariance matches'
            passed = passed + 1
        else
            write(*,'(A,A,A)') '  FAIL: ', label, ' covariance mismatch'
            failed = failed + 1
        end if
    end subroutine compare_cov

    subroutine write_samples(fname, samples, n)
        character(len=*), intent(in) :: fname
        real(DP), intent(in) :: samples(:,:)
        integer, intent(in) :: n
        integer :: u, k
        open(newunit=u, file=trim(fname), status='replace', action='write')
        do k = 1, n
            write(u, '(6ES22.14)') samples(:, k)
        end do
        close(u)
    end subroutine write_samples

end program test_uq_crtbp_comparison
```

- [ ] **Step 2: Add test target to fpm.toml**

After the existing `[[test]]` entry for `test_stt_crtbp`:

```toml
[[test]]
name = "test_uq_crtbp_comparison"
source-dir = "test"
main = "test_uq_crtbp_comparison.f90"
```

- [ ] **Step 3: Build and run test**

```bash
cd POD_Fortran && fpm build
fpm test test_uq_crtbp_comparison
```

Expected: test compiles and runs. Output shows comparison results. MC will be slowest; ADS should produce reasonable agreement.

- [ ] **Step 4: Commit**

```bash
git add test/test_uq_crtbp_comparison.f90 fpm.toml
git commit -m "test: add unified CRTBP comparison test (MC/DA/STT/ADS)"
```

---

### Task 7: Integration Verification and Cleanup

**Files:** (all created/modified above)

- [ ] **Step 1: Run full build**

```bash
cd POD_Fortran && fpm build
```

Expected: All modules compile without errors.

- [ ] **Step 2: Run test and verify output**

```bash
cd POD_Fortran && fpm test test_uq_crtbp_comparison
```

Expected:
- ADS splitting produces n_patches > 0
- ADS mean within reasonable tolerance of MC
- ADS covariance within ~15% of MC
- Timing shows ADS is faster than MC

- [ ] **Step 3: Run existing tests to verify no regressions**

```bash
cd POD_Fortran && fpm test test_stt_crtbp
```

Expected: existing STT test still passes.

- [ ] **Step 4: Final commit**

```bash
git add -A
git commit -m "chore: integration verification, all tests pass"
```
