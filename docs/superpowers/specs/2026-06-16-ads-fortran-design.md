# ADS Uncertainty Propagation in POD_Fortran — Design Spec

## Overview

Implement Automatic Domain Splitting (ADS) uncertainty propagation in POD_Fortran, translating the C++ reference implementation (`ADS_cpp_Core_file/`) from DACE-based C++ to DACE-based Fortran. First target: CRTBP dynamics.

## Files Changed

**New (4):**
```
src/lib/math/
  pod_ads_split_module.f90          ← Core ADS: SplittingHistory, Patch, Manifold

src/lib/uncertainty/propagation/
  pod_uq_prop_ads_module.f90        ← uq_ads_propagator (extends uq_propagator_base)
  pod_uq_crtbp_ads_module.f90       ← CRTBP-specific: crtbp_ads_propagate(_deviates)

test/
  test_uq_crtbp_comparison.f90      ← Unified comparison: DA / STT / ADS vs MC
```

**Modified (2):**
```
external/dace_build/dace_wrapper.cpp        ← Add fdace_estim_norm C wrapper
src/lib/system/pod_dace_classes.f90         ← Add Fortran binding for estimNorm
```

## Dependency Graph

```
pod_dace_classes (DACE C bindings)  ←  add estimNorm binding
       ↓
pod_ads_split_module (core S/H/P/M)
       ↓
pod_uq_crtbp_ads_module ← pod_uq_prop_ads_module ← pod_uq_prop_base_module
       ↓
test_uq_crtbp_comparison (+ existing STT/DA/MC imports)
```

## Step 0: Add `estimNorm` to DACE Wrapper

### C wrapper (dace_wrapper.cpp)

The DACE C++ class has `DA::estimNorm(int varIdx, int minOrder, int maxOrder)` which returns `std::vector<double>` — the norm of coefficients grouped by order, restricted to variable `varIdx` (0 = all variables).

Add:
```cpp
int fdace_estim_norm(int h_in, int var_idx, int max_order, double* norms) {
    std::vector<double> result = da_registry[h_in].estimNorm(var_idx, 0, max_order + 1);
    int sz = result.size();
    for (int i = 0; i < sz; ++i) norms[i] = result[i];
    return sz;  // number of norm values written
}
```

### Fortran binding (pod_dace_classes.f90)

```fortran
function c_fdace_estim_norm(h, var_idx, max_order, norms) bind(C, name="fdace_estim_norm") result(sz)
    import :: c_int, c_double
    integer(c_int), value :: h, var_idx, max_order
    real(c_double), intent(out) :: norms(*)
    integer(c_int) :: sz
end function
```

Usage: `sz = c_fdace_estim_norm(handle, 0, max_order, norms_array)` returns number of entries; `norms_array(sz)` = top-order norm.

## Module 1: `pod_ads_split_module` — Core ADS

Style: procedural (no type-bound procedures), matching POD_Fortran conventions.

### Types

```fortran
type :: splitting_history_type
    integer, allocatable :: entries(:)   ! ±dir, dir ∈ 1..6
end type

type :: patch_type
    type(AlgebraicVector) :: da_vec       ! DA expansion (6 components)
    type(splitting_history_type) :: history
end type

type :: manifold_type
    type(patch_type), allocatable :: patches(:)
    integer :: n_patches
end type
```

### Subroutines/Functions

**SplittingHistory:**
- `sh_push(history, dir)` — append ±dir to entries
- `sh_pop(history)` — remove last entry
- `sh_count(history, n)` → integer — total splits (n=0) or count of dir |n|
- `sh_replay(history, obj)` — rebuild sub-domain DA vector from initial identity via successive eval (inout obj)
- `sh_center(history)` → real(6) — sub-box center in [-1,1]^6
- `sh_width(history)` → real(6) — sub-box width
- `sh_contain(history, pt)` → logical — pt ∈ this sub-box?
- `sh_map_point(history, pt)` — transform pt in-place: global[-1,1]^6 → local[-1,1]^6

**Key algorithm for `sh_replay`:**
```
x = identity_DA_vector(6)            ! x(i) = da_var(i)
for each entry s in history:
    n = abs(s) - 1, sign = s/abs(s)
    x(n) = 0.5*sign + 0.5*da_var(n+1)
    obj = obj.eval(x)                 ! substitute x into obj
    x(n) = da_var(n+1)                ! restore
```
Corresponds exactly to C++ `SplittingHistory::replay`.

**Patch:**
- `patch_init(p, da_vec, history)` / `patch_destroy(p)`
- `patch_get_trunc_err(p, max_order, errors)` — errors(6): for each component, call `estimNorm(0, 0, max_order+1)` → take last element (= highest-order norm)
- `patch_get_split_dir(p, comp, max_order)` → integer direction (1..6): for each var i, call `estimNorm(i, 0, max_order+1)` → last element; direction with max value wins
- `patch_split(p, dir, left, right)` — split along `dir`:
  1. left: push -dir to history, substitute `obj(dir)=-0.5+0.5*da_var(dir)`, eval
  2. right: push +dir to history, substitute `obj(dir)=0.5+0.5*da_var(dir)`, eval

**Manifold:**
- `mf_push(m, p)` / `mf_pop_front(m, p)`
- `ads_get_split_domain(initial_da_vec, err_toll, n_split_max, mu, t_end, max_order, result_manifold)` — BFS domain splitting

### Key Algorithm: `ads_get_split_domain`

```
queue ← [initial_patch]          ! initial_patch has empty history
while queue not empty:
    p ← pop_front(queue)
    f_da_vec ← da_adaptive_integrate_crtbp(p.da_vec, 0.0, t_end, ...)  ! explicit CRTBP call
    f.history ← p.history
    errors ← get_truncation_errors(f, max_order)
    relative_errors ← max(0.0, errors - err_toll)   ! per-component, scalar if err_toll is scalar
    if maxval(relative_errors) == 0 or sh_count(p.history) >= n_split_max:
        results ← push(f)
    else:
        pos ← maxloc(relative_errors)
        dir ← get_splitting_direction(f, pos, max_order)
        (left, right) ← split(p, dir)
        queue ← push(left, right)
return results
```

## Module 2: `pod_uq_prop_ads_module` — Propagator Type

Extends `uq_propagator_base` (same pattern as `uq_stt_propagator`).

```fortran
type, extends(uq_propagator_base) :: uq_ads_propagator
    integer  :: da_order = 2
    integer  :: n_split_max = 12
    real(DP) :: err_toll(6) = 1.0d-4       ! per-component tolerance
    real(DP) :: mu = 0.012153614091892_DP
    real(DP) :: ads_abs_tol = 1.0d-12
    real(DP) :: ads_rel_tol = 1.0d-12
    real(DP) :: ads_dt_min = 1.0d-6
    real(DP) :: ads_dt_max = 3600.0_DP
    integer  :: ads_max_steps = 100000
    type(manifold_type) :: domain           ! split result
    integer  :: n_patches = 0               ! after splitting
contains
    procedure :: propagate       => ads_propagate
    procedure :: get_method_name => ads_get_method_name
    procedure :: set_ads_order
    procedure :: set_ads_tolerances
    procedure :: set_ads_mu
end type
```

- `set_ads_tolerances` allows setting abs_tol, rel_tol, dt_min, dt_max, max_steps, err_toll, n_split_max (all optional)
- `get_method_name` → "ADS-N" where N = da_order

## Module 3: `pod_uq_crtbp_ads_module` — CRTBP Binding

Public API mirroring existing `pod_uq_crtbp_da_module` / `pod_uq_crtbp_mc_module`:

```fortran
public :: crtbp_ads_propagate, crtbp_ads_propagate_deviates
```

### `crtbp_ads_propagate_deviates` flow:

1. DACE init: `dace_initialize(da_order, 6)`
2. Build initial DA state: `x0%elements(i) = nominal(i) + da_var(i)` → initial `AlgebraicVector`
3. Call `ads_get_split_domain(x0, err_toll, n_split_max, mu, t_end, da_order, manifold)`
4. For each deviates(:, k):
   - Find patch: test `sh_contain(patch%history, deviates(:, k))` across manifold
   - Map point: `sh_map_point(patch%history, pt_local)` where pt_local starts as deviates(:,k)
   - Evaluate: `final_samples(:, k) = patch%da_vec%eval(pt_local)`
5. Compute mean and covariance
6. Cleanup: `dace_pop_to()`

Signature mirrors existing modules:
```fortran
subroutine crtbp_ads_propagate_deviates(nominal_state, deviates, mu, t_end, &
    da_order, n_split_max, err_toll, rel_tol, abs_tol, dt_min, dt_max, max_steps, &
    final_samples, final_mean, final_cov, n_patches_out, verbose)
```

Same for `crtbp_ads_propagate`: generates particles from covariance, then same pipeline.

## Module 4: `test_uq_crtbp_comparison.f90` — Unified Test

Extends the Test 6 pattern from `test_stt_crtbp.f90` with ADS added.

Test parameters (same as existing test):
- DRO orbit: mu=0.012153614091892, x0=[1.1309..., 0, 0, 0, -0.4654..., 0]
- Propagation time: 1.5 × DRO period
- Deviates: `rand_list200km_0.7ms.txt`
- n_dev: 100k (DA/STT/ADS), n_mc: 100k (MC reference)
- DA/ADS order: 2, STT order: 2

Test items:
1. ADS splitting — verify n_patches > 0
2. Mean comparison — ADS vs DA vs MC (max absolute diff)
3. Covariance comparison — ADS vs DA vs MC (relative Frobenius error)
4. Timing — wall-clock time for each of MC/DA/STT/ADS, speedup vs MC
5. Output files — sample results for each method (mirror existing pattern)

## Non-Goals (this iteration)

- No 2D envelope export
- No JSON config parsing
- No multi-threaded sample evaluation (sequential first)
- Only CRTBP dynamics
