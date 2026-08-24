# ADS Eval-Based Split Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the Fortran ADS `translateVariable` split path with the same DA-vector composition used by the C++ reference, without modifying DACE core.

**Architecture:** `patch_split` builds a full identity DA map, changes the selected coordinate to `0.5*x +/- 0.5`, and composes every parent output into preallocated child handles through the existing `fdace_eval_da_vec` wrapper. `patch_init` transfers DA storage with `move_alloc`, so no allocated handle is overwritten or leaked.

**Tech Stack:** Fortran 2008, ISO C binding, DACE C++ wrapper, fpm tests.

## Global Constraints

- Do not modify or rebuild DACE core.
- Keep `patch_split(parent, direction, left, right)` unchanged.
- Keep Patch output count separate from the active DACE variable count.
- Do not change ADS normalization or HFEM force-model scaling.
- Codex edits source and tests but does not compile or run; the user performs verification.

---

### Task 1: Strengthen the dynamic-dimension regression

**Files:**
- Modify: `test/test_ads_core_dynamic_dimensions.f90`

**Interfaces:**
- Consumes: `active_da_count()`, `patch_split`, `AlgebraicVector%eval(real(:))`.
- Produces: a regression that distinguishes left/right DA composition and detects leaked active handles.

- [ ] **Step 1: Record the handle baseline**

Import `active_da_count`, declare `active_before`, and record it immediately after DACE initialization:

```fortran
use pod_dace_classes, only: AlgebraicVector, DA, dace_initialize, &
    da_exp_sub, active_da_count, assignment(=)
integer :: n_fail, direction, active_before

call dace_initialize(2, 7)
active_before = active_da_count()
```

- [ ] **Step 2: Correct and extend child-value assertions**

The left global point is `x7=-0.75`, so the order-two expansion of `exp(x7)` has a negative linear term. Add the corresponding right-child evaluation:

```fortran
values = left%da_vec%eval(point)
call assert_close(values(1), 1.0_DP - 0.75_DP + 0.5_DP * 0.75_DP**2, 1.0e-13_DP, &
    'translated left polynomial retains the global value', n_fail)

point = 0.0_DP
point(7) = 0.75_DP
call assert_true(sh_contain(right%history, point), &
    'right child contains a global unit-domain point', n_fail)
call sh_map_point(right%history, point)
call assert_close(point(7), 0.5_DP, 1.0e-14_DP, &
    'global point maps to the right local coordinate', n_fail)
values = right%da_vec%eval(point)
call assert_close(values(1), 1.0_DP + 0.75_DP + 0.5_DP * 0.75_DP**2, 1.0e-13_DP, &
    'translated right polynomial retains the global value', n_fail)
```

- [ ] **Step 3: Assert cleanup returns to the baseline**

After destroying all DA-bearing objects, add:

```fortran
call assert_equal_int(active_da_count(), active_before, &
    'ADS split releases every temporary and Patch DA handle', n_fail)
```

- [ ] **Step 4: Preserve the observed RED evidence**

Do not run this command as Codex. The user already observed the current failure:

```text
fpm test test_ads_core_dynamic_dimensions
```

Expected pre-fix result: abort in `da_translate_variable`/`patch_split` with `free(): invalid next size`.

---

### Task 2: Add no-temporary DA-vector composition

**Files:**
- Modify: `src/lib/system/pod_dace_classes.f90`
- Test: `test/test_ads_core_dynamic_dimensions.f90`

**Interfaces:**
- Consumes: existing `c_fdace_eval_da_vec(hi, h_map, n_vars, ho)`.
- Produces: `vector_eval_da_vec_sub(this, map_vec, res)` with two input `AlgebraicVector` objects and one preallocated/reusable output `AlgebraicVector`.

- [ ] **Step 1: Export the new subroutine**

Add it to the module public list:

```fortran
public :: vector_eval_da_vec_sub
```

- [ ] **Step 2: Implement direct composition into existing handles**

Place this beside `vector_eval_da_vec`:

```fortran
subroutine vector_eval_da_vec_sub(this, map_vec, res)
    class(AlgebraicVector), intent(in) :: this, map_vec
    type(AlgebraicVector), intent(inout) :: res
    integer(c_int), allocatable :: h_map(:)
    integer :: i

    if (.not. allocated(this%elements)) &
        error stop 'vector_eval_da_vec_sub: source vector is not allocated'
    if (.not. allocated(map_vec%elements)) &
        error stop 'vector_eval_da_vec_sub: DA map is not allocated'

    allocate(h_map(map_vec%size))
    do i = 1, map_vec%size
        h_map(i) = map_vec%elements(i)%handle
    end do

    if (.not. allocated(res%elements) .or. res%size /= this%size) then
        call res%destroy()
        call res%init(this%size)
    end if

    do i = 1, this%size
        call c_fdace_eval_da_vec(this%elements(i)%handle, h_map, &
            int(map_vec%size, c_int), res%elements(i)%handle)
    end do
end subroutine vector_eval_da_vec_sub
```

This must not use the existing function result or overloaded assignment; both would create finalizable DA temporaries.

- [ ] **Step 3: Static interface verification**

Confirm the wrapper argument order remains:

```fortran
c_fdace_eval_da_vec(input_handle, map_handles, map_size, output_handle)
```

No C++ wrapper or DACE source change is part of this task.

---

### Task 3: Replace ADS translation with C++-style composition

**Files:**
- Modify: `src/lib/math/pod_ads_split_module.f90`
- Test: `test/test_ads_core_dynamic_dimensions.f90`

**Interfaces:**
- Consumes: `dace_max_variables`, `da_mul`, `da_add`, `vector_eval_da_vec_sub`.
- Produces: unchanged `patch_split(p, dir, left, right)` behavior implemented through DA composition.

- [ ] **Step 1: Replace ADS imports**

Remove the ADS import of `da_translate_variable` and import:

```fortran
use pod_dace_classes, only: DA, AlgebraicVector, dace_max_variables, &
    da_mul, da_add, vector_eval_da_vec_sub
```

Retain other existing imports used elsewhere in the module.

- [ ] **Step 2: Build an identity map without function temporaries**

Inside `patch_split`, use:

```fortran
type(AlgebraicVector) :: map_vec
type(DA) :: affine_da
integer :: i, n_components, nvars

n_components = p%da_vec%size
nvars = dace_max_variables()
if (dir < 1 .or. dir > nvars) &
    error stop 'patch_split: split direction is outside the active DACE variables'
if (.not. allocated(p%da_vec%elements)) &
    error stop 'patch_split: parent Patch has no DA vector'

call map_vec%init(nvars)
do i = 1, nvars
    call map_vec%elements(i)%destroy()
    call map_vec%elements(i)%init_var(i)
end do
call affine_da%init()
```

- [ ] **Step 3: Compose the left child**

Replace `translateVariable` with:

```fortran
left%history = p%history
call sh_push(left%history, -dir)
call da_mul(map_vec%elements(dir), 0.5_DP, affine_da)
call da_add(affine_da, -0.5_DP, map_vec%elements(dir))
call vector_eval_da_vec_sub(p%da_vec, map_vec, left%da_vec)
```

- [ ] **Step 4: Restore the selected identity coordinate and compose the right child**

```fortran
call map_vec%elements(dir)%destroy()
call map_vec%elements(dir)%init_var(dir)

right%history = p%history
call sh_push(right%history, dir)
call da_mul(map_vec%elements(dir), 0.5_DP, affine_da)
call da_add(affine_da, 0.5_DP, map_vec%elements(dir))
call vector_eval_da_vec_sub(p%da_vec, map_vec, right%da_vec)

call affine_da%destroy()
call map_vec%destroy()
```

- [ ] **Step 5: Delete the obsolete split implementation**

Remove `temp_vec`, `saved_hist`, `new_handle`, `DIAGNOSTIC_SPLIT`, both raw-handle loops, and every call to `da_translate_variable` from `patch_split`.

---

### Task 4: Remove Patch-construction handle leaks

**Files:**
- Modify: `src/lib/math/pod_ads_split_module.f90`
- Test: `test/test_ads_core_dynamic_dimensions.f90`

**Interfaces:**
- Consumes: allocatable `AlgebraicVector%elements` and `%h_list`.
- Produces: the existing move-style `patch_init(p, da_vec, history)` without allocating and overwriting destination handles.

- [ ] **Step 1: Transfer the DA storage with `move_alloc`**

Replace the destination initialization and raw handle loop with:

```fortran
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
```

Remove the obsolete component loop and its loop index.

- [ ] **Step 2: Preserve history without aliasing**

Keep the existing local `hist_copy` captured before the `intent(out)` destination is changed, then assign it after the storage transfer:

```fortran
if (present(history)) p%history = hist_copy
```

- [ ] **Step 3: Review ownership invariants**

After `patch_init`, the source has no allocated DA element array and size zero. The Patch exclusively owns every transferred handle. `source%destroy()` remains a safe no-op in callers.

---

### Task 5: Static verification and user-run handoff

**Files:**
- Review: `src/lib/system/pod_dace_classes.f90`
- Review: `src/lib/math/pod_ads_split_module.f90`
- Review: `test/test_ads_core_dynamic_dimensions.f90`

**Interfaces:**
- Consumes: all changes from Tasks 1-4.
- Produces: a source-only handoff for user compilation.

- [ ] **Step 1: Run source-level checks only**

Codex may run read-only checks:

```text
git diff --check
rg -n "da_translate_variable" src/lib/math/pod_ads_split_module.f90
rg -n "vector_eval_da_vec_sub" src/lib/system/pod_dace_classes.f90 src/lib/math/pod_ads_split_module.f90
```

Expected: no whitespace errors; no ADS call to `da_translate_variable`; the new subroutine has one definition and one ADS call.

- [ ] **Step 2: User verifies the focused regression**

The user runs:

```text
fpm test test_ads_core_dynamic_dimensions
```

Expected final line:

```text
PASS: ADS Core separates 7 DA variables from 3 Patch components
```

- [ ] **Step 3: User verifies dependent ADS/HFEM tests**

The user runs:

```text
fpm test test_uq_crtbp_comparison
fpm test test_hfem_ads_propagation
```

Expected: both compile and execute without entering `da_translate_variable`; numerical/timing diagnostics remain test-controlled.
