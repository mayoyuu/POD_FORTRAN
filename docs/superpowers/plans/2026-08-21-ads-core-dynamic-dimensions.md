# ADS Core Dynamic Dimensions Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Restore the Fortran ADS Core's C++ behavior by separating the runtime DA independent-variable count from the runtime Patch output-component count.

**Architecture:** DACE remains the authority for the global independent-variable count, exposed through a small public Fortran wrapper. Splitting-history operations query that count, while Patch operations use `da_vec%size`. Existing six-variable CRTBP callers remain source-compatible.

**Tech Stack:** Fortran 2008, DACE C bindings, fpm.

## Global Constraints

- Do not compile or run; the user will perform all RED/GREEN verification.
- Write the regression test before changing production modules.
- Preserve the existing six-variable CRTBP API and behavior.
- Do not touch the user's modified `rename_to_pod.sh` or `setup_env.sh`.
- Do not create a Git commit unless the user explicitly asks.

---

### Task 1: Add a dynamic-dimension ADS regression test

**Files:**
- Create: `test/test_ads_core_dynamic_dimensions.f90`
- Modify: `fpm.toml`

**Interfaces:**
- Consumes: `dace_initialize`, `da_var`, `patch_init`, `patch_split`, `patch_get_split_dir`, `sh_center`, `sh_width`, `sh_contain`, and `sh_map_point`.
- Produces: a test requiring seven DA variables with a three-component Patch.

- [ ] **Step 1: Write the failing test**

The test initializes DACE with seven variables, constructs a three-component
map whose first component is `DA(7)**2`, verifies that the output size remains
three, requires split-direction selection to return seven, splits along variable
seven, and checks the child histories and local point mapping.

```fortran
call dace_initialize(2, 7)
call source%init(3)
source%elements(1) = da_var(7) * da_var(7)
source%elements(2) = 2.0_DP + da_var(2)
source%elements(3) = -3.0_DP
call patch_init(parent, source)
call assert_equal_int(parent%da_vec%size, 3, 'Patch output size')
call assert_equal_int(patch_get_split_dir(parent, 1, 2), 7, 'split direction')
call patch_split(parent, 7, left, right)
center = sh_center(left%history)
width = sh_width(left%history)
call assert_equal_int(size(center), 7, 'history center size')
call assert_close(center(7), -0.5_DP, 1.0e-14_DP, 'left center')
call assert_close(width(7), 1.0_DP, 1.0e-14_DP, 'left width')
point = 0.0_DP
point(7) = -0.75_DP
call assert_true(sh_contain(left%history, point), 'left contains point')
call sh_map_point(left%history, point)
call assert_close(point(7), -0.5_DP, 1.0e-14_DP, 'left local coordinate')
```

- [ ] **Step 2: Register the test**

```toml
[[test]]
name = "test_ads_core_dynamic_dimensions"
source-dir = "test"
main = "test_ads_core_dynamic_dimensions.f90"
```

- [ ] **Step 3: User verifies RED**

Run: `fpm test test_ads_core_dynamic_dimensions`

Expected before production changes: compilation failure caused by fixed-shape
six-element history APIs or a runtime assertion that direction seven is not
selected.

---

### Task 2: Expose the active DACE independent-variable count

**Files:**
- Modify: `src/lib/system/pod_dace_classes.f90`
- Test: `test/test_ads_core_dynamic_dimensions.f90`

**Interfaces:**
- Produces: `integer function dace_max_variables()`.

- [ ] **Step 1: Add the public API**

Add `dace_max_variables` to the module public list and implement:

```fortran
integer function dace_max_variables() result(nvars)
    nvars = int(c_fdace_get_max_variables())
end function dace_max_variables
```

- [ ] **Step 2: Keep the existing C binding private**

The ADS module imports only `dace_max_variables`; it does not import
`c_fdace_get_max_variables` directly.

---

### Task 3: Generalize SplittingHistory to the active DACE dimension

**Files:**
- Modify: `src/lib/math/pod_ads_split_module.f90`
- Test: `test/test_ads_core_dynamic_dimensions.f90`

**Interfaces:**
- `sh_center(history)` returns `real(DP), allocatable :: c(:)`.
- `sh_width(history)` returns `real(DP), allocatable :: w(:)`.
- `sh_contain(history, pt)` accepts `real(DP) :: pt(:)`.
- `sh_map_point(history, pt)` accepts `real(DP) :: pt(:)`.

- [ ] **Step 1: Query DACE instead of using literal six**

Import `dace_max_variables`. Allocate center, width, and replay identity arrays
with `dace_max_variables()`.

```fortran
nvars = dace_max_variables()
allocate(c(nvars), w(nvars))
c = 0.0_DP
w = 2.0_DP
```

- [ ] **Step 2: Make point arguments assumed-shape**

Containment returns `.false.` when `size(pt) /= dace_max_variables()`.
Mapping leaves an invalid-size point unchanged.

- [ ] **Step 3: Remove the empty-history replay leak**

Destroy the temporary identity vector before returning for an empty history.

---

### Task 4: Generalize Patch operations independently of DA dimension

**Files:**
- Modify: `src/lib/math/pod_ads_split_module.f90`
- Test: `test/test_ads_core_dynamic_dimensions.f90`

**Interfaces:**
- Patch output count is always `p%da_vec%size`.
- Split-direction search count is always `dace_max_variables()`.
- `patch_get_trunc_err` accepts `errors(:)`.

- [ ] **Step 1: Preserve the source vector size in `patch_init`**

```fortran
n_components = da_vec%size
call p%da_vec%init(n_components)
do i = 1, n_components
    p%da_vec%elements(i)%handle = da_vec%elements(i)%handle
    da_vec%elements(i)%handle = -1
end do
```

- [ ] **Step 2: Loop over output components for truncation errors**

Validate `size(errors) >= p%da_vec%size`, initialize the caller array to zero,
then estimate each actual output component.

- [ ] **Step 3: Loop over independent variables for direction selection**

```fortran
do i = 1, dace_max_variables()
    call da_estim_norm(p%da_vec%elements(comp)%handle, i, order, top_norm)
end do
```

- [ ] **Step 4: Preserve output size during splitting**

Both diagnostic-copy and translated children allocate `temp_vec` with
`p%da_vec%size` and translate every output component along the requested DA
variable, including direction seven.

---

### Task 5: Manual verification checkpoint

- [ ] **Step 1: User verifies GREEN for the new test**

Run: `fpm test test_ads_core_dynamic_dimensions`

Expected: `PASS: ADS Core separates 7 DA variables from 3 Patch components`.

- [ ] **Step 2: User verifies six-variable regression tests**

Run:

```text
fpm test test_ads_domain_scale
fpm test test_uq_crtbp_comparison
```

Expected: existing tests pass with unchanged six-variable behavior.

