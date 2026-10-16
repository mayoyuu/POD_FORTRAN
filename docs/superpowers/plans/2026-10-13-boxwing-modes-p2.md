# Box-Wing Modes P2 Implementation Plan

> **For agentic workers:** Use superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox syntax for tracking.

**Goal:** Implement and verify the mode-table lifecycle, Real/DA surrogate SRP and eight-variable physical coordinate map using the audited P1 dataset.

**Architecture:** A Real-only table module validates finalized metadata and interpolates absolute TDB modes. Real/DA force dispatchers share the proxy formula and use separate explicit coefficient state. The coordinate map resolves NONE/ETA/MODES and exposes patch-local s/u centers/spans. P3 will connect the complete HFEM runner, study-point and snapshot pipeline.

**Tech Stack:** Fortran/FPM, SPICE, existing DACE classes and ADS splitter; no new external dependency.

## Global Constraints

- P1 mode data remain immutable; full-sun/J2000/EARTH/TDB and km/km/s/km/s2 conventions remain fixed.
- Fortran validates finalized markers, schema, units, row/grid/vector identities, epoch and common Cannonball physical parameters. Python remains responsible for SHA-256 verification before formal runs.
- Missing table, unset coefficients, extrapolation, invalid indices/domain or mixed eta/mode uncertainty fail explicitly.
- Preserve legacy 6D/7D API wrappers and the existing force behavior.
- Root s/u = 0.5*(1+q7/q8); patch center = 0.5*(1+c7/c8), span = 0.25*w7/w8.
- Node endpoints and Real/DA constants must agree to relative 1e-12; patch mapping must retain global physical values and release all owned DA handles.
- Core force functions and maps are P2; formal 8D propagation runner, snapshots and long-arc convergence are P3/P4.

## Task 1: Typed metadata and mode table

**Files:** Create `src/lib/data/pod_json_object_module.f90`, `src/lib/forcemodel/pod_srp_mode_module.f90`, `test/test_srp_mode_table.f90`.

**Interfaces:** `load_srp_mode_table(csv_file,meta_file,expected_epoch0,status,message)`, `clear_srp_mode_table()`, `srp_mode_table_is_loaded()`, `interpolate_srp_modes(epoch,eE,eM,status,message)`.

- [x] Write fixture tests for node/midpoint/endpoints, missing/fake metadata, schema/unit/time/finite/vector mismatches, duplicate JSON fields, physical parameter and epoch mismatch, and invalidated reloads.
- [x] Run `source setup_env.sh && fpm test test_srp_mode_table`; confirm missing modules fail compilation.
- [x] Implement a typed JSON object reader and Real-only table module. Load clears prior storage first; interpolation clamps only <=1e-6 s endpoint roundoff and never extrapolates.
- [x] Run the table tests and confirm no DA handle dependency.

## Task 2: Explicit eight-variable coordinates

**Files:** Modify `src/lib/uncertainty/propagation/pod_uq_ads_coordinates_module.f90`; create `test/test_srp_mode_parameter_map.f90`.

**Interfaces:** Add ADS_PARAM_AUTO/NONE/SRP_ETA/SRP_MODES, `ads_physical_parameters_type`, optional parameterization on map construction, generic physical/unit conversions and `ads_patch_mode_parameters(map,center,width,s0,ds,u0,du,status)`.

- [x] Write tests for AUTO legacy resolution, 8D corner/center/random round trips, triangle inequalities, invalid domains, and actual ADS left/right splits in variables seven/eight.
- [x] Run `fpm test test_srp_mode_parameter_map`; confirm missing interfaces fail compilation.
- [x] Implement explicit parameterization and retain eta-scalar wrappers. Reject ambiguous mixed modes/eta; clip physical boundaries only within 1000*epsilon.
- [x] Re-run new and existing `test_ads_coordinate_mapping` and dynamic ADS core tests.

## Task 3: Configuration, lifecycle and surrogate force

**Files:** Modify `pod_config_module.f90`, `pod_engine_module.f90`, `pod_force_model_module.f90`, `pod_da_force_model_module.f90`; create `test/test_srp_modes_force_consistency.f90` and Python negative-path driver.

**Interfaces:** Add config paths and model selector; Real `set/clear_srp_mode_coefficients_real`; DA `set/clear_srp_mode_uncertainty`; internal Cannonball full-sun helpers independent of model dispatch.

- [x] Write tests for config semantics, init clearing, absent table/unset coefficient/invalid setter failures, P1's 361 three-endpoint checks, Real/DA constants and local coefficient sensitivity.
- [x] Run focused tests and confirm the new APIs are missing.
- [x] Add model config and path validation; engine clears table and coefficient state before reading configuration, without loading a table before the OPM epoch is known.
- [x] Extract Cannonball helpers without changing classic eta/ballistic override behavior. Mode branch rejects mixed legacy uncertainty/overrides and evaluates CB + norm(CB)*(s*eE+(1-s)*u*eM).
- [x] Verify all 361 reference nodes, representative off-reference states, short 8D force evaluations, split-local coefficients and repeated-call handle counts.

## Task 4: Regression and result

- [x] Run focused existing Real/DA, configuration, 6D/7D coordinate and ADS history tests, plus P1 Python suites and mode-table re-audits.
- [x] Inspect exact diff and `git diff --check`; record peak endpoint/constant errors and handle counts in a P2 report.
- [x] Update phase status in the design document; deliver the core implementation with its verified scope. Do not claim the unimplemented P3 runner or long-arc accuracy passes.

The existing isolated checkout is reused. Keep P1 additions and artifacts intact; leave commits/merge for a separately authorized integration step.
