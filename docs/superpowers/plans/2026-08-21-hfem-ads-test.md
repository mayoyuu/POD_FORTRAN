# HFEM ADS Test Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add one manually runnable HFEM ADS test covering both six orbital variables and six orbital variables plus the SRP scale parameter `eta_srp`.

**Architecture:** A test-local HFEM BFS driver reuses the now-generic ADS Core and the production HFEM DA RKF78 integrator. Orbital uncertainty is embedded in the initial DA state; the non-state SRP parameter is reconstructed from each patch history before each integration.

**Tech Stack:** Fortran 2008, DACE, HFEM force model, RKF78, fpm.

## Global Constraints

- Begin only after the dynamic ADS Core changes are written.
- Do not compile or run; the user performs verification.
- Default case is `L1Halo-1`, seven days, 100 km and 0.3 m/s 1-sigma orbital uncertainty.
- The orbital ADS box is plus or minus 3 sigma.
- The SRP multiplier is `1 + eta_srp`, with `eta_srp` in `[-0.1, 0.1]`.
- Do not create a Git commit unless the user explicitly asks.

---

### Task 1: Add patch-local SRP DA scaling

**Files:**
- Modify: `src/lib/forcemodel/pod_da_force_model_module.f90`
- Test: `test/test_hfem_ads_propagation.f90`

**Interfaces:**
- `set_srp_scale_uncertainty(var_index, nominal_scale, da_span)` where the last two arguments are optional.
- Existing two-argument callers remain valid and retain `da_span = 1`.

- [ ] **Step 1: Add an optional span**

Store a module value `srp_scale_da_span = 1.0_DP`. Reset it to one in both the
setter and clearer, then override it when `da_span` is present.

- [ ] **Step 2: Apply the span to the formal variable**

Construct the force multiplier as:

```fortran
srp_multiplier_da = 1.0_DP + srp_scale_nominal + &
                    srp_scale_da_span * srp_scale_da
```

Then multiply the nominal SRP factor by `srp_multiplier_da`.

---

### Task 2: Write the HFEM ADS test before completing its driver

**Files:**
- Create: `test/test_hfem_ads_propagation.f90`
- Modify: `fpm.toml`

**Interfaces:**
- Consumes: generic ADS Core, `load_initial_opm`, `da_adaptive_step_integrate`, HFEM epoch and SRP setters.
- Produces: a single test executable running `orbit-6D` and `orbit+SRP-7D`.

- [ ] **Step 1: Define the default case and assertions**

Use DA order four, seven-day propagation, RKF78, maximum split depth eight,
position/velocity domain half-widths of 300 km and 0.0009 km/s, and SRP half-
width 0.1. Require accepted patches, complete deterministic-point coverage,
finite results, correct zero-duration physical mapping, and nonzero final SRP
sensitivity.

- [ ] **Step 2: Register the test**

```toml
[[test]]
name = "test_hfem_ads_propagation"
source-dir = "test"
main = "test_hfem_ads_propagation.f90"
```

---

### Task 3: Implement the test-local HFEM ADS BFS driver

**Files:**
- Modify: `test/test_hfem_ads_propagation.f90`

**Interfaces:**
- `build_hfem_ads_domain(..., n_variables, use_srp_uncertainty, manifold, stats)`.
- `evaluate_hfem_ads_point(manifold, point_unit, state_out, found)`.

- [ ] **Step 1: Build the normalized orbital initial map**

```fortran
state_da_0%elements(i) = (state0(i) + domain_half_width(i) * da_var(i)) / config%LU
state_da_0%elements(i+3) = (state0(i+3) + domain_half_width(i+3) * da_var(i+3)) / config%VU
```

- [ ] **Step 2: Reconstruct the SRP map per queued patch**

For each patch history, compute `center = sh_center(history)` and
`width = sh_width(history)`, then call:

```fortran
call set_srp_scale_uncertainty(7, &
    ETA_MEAN + ETA_HALF_WIDTH * center(7), &
    ETA_HALF_WIDTH * width(7) / 2.0_DP)
```

The 6D scenario calls `clear_srp_scale_uncertainty()`.

- [ ] **Step 3: Integrate and convert accepted maps to physical units**

Call `da_adaptive_step_integrate` with nondimensional time and RKF78. Multiply
the final position DA components by `config%LU` and velocity components by
`config%VU` before truncation-error estimation and acceptance.

- [ ] **Step 4: Split the input patch, not the propagated patch**

Use the propagated physical map only for error estimation and accepted output.
When error exceeds tolerance, obtain the direction from that propagated map but
call `patch_split` on the queued input map so each child is re-integrated.

- [ ] **Step 5: Evaluate deterministic unit-domain points**

Use `sh_contain`, `sh_map_point`, and compiled accepted Patch maps. Evaluate the
center, orbital interior points, and `eta_srp` endpoints `-1` and `+1`.

---

### Task 4: Manual verification checkpoint

- [ ] **Step 1: User builds**

Run: `fpm build`

Expected: successful compilation of the generic ADS Core and HFEM ADS test.

- [ ] **Step 2: User runs the Core test first**

Run: `fpm test test_ads_core_dynamic_dimensions`

Expected: PASS before starting the expensive HFEM propagation.

- [ ] **Step 3: User runs the HFEM test**

Run: `fpm test test_hfem_ads_propagation`

Expected: both 6D and 7D scenarios accept at least one patch, all deterministic
points are covered and finite, and the 7D case reports nonzero SRP sensitivity.

