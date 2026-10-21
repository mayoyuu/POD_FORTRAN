# Box-Wing Modes P3 Implementation Plan

> **For agentic workers:** Use superpowers:executing-plans task-by-task, with test-driven-development and verification-before-completion. Preserve the current uncommitted P1/P2 work; no commit, merge, push or campaign mutation in this phase.

**Goal:** Connect the eight-variable C2 propagation pipeline and pass a one-hour ADS versus independent Real smoke, including reconstructible snapshots and audited publication.

**Architecture:** Extend the existing HFEM domain using explicit NONE/ETA/MODES dispatch. Dimension-aware study points and Real mode probes feed a separate C2 runner. New v2 snapshots are opt-in for C2 while existing callers retain v1; a Python wrapper audits identities, outputs and snapshots before atomic directory publication.

**Tech Stack:** Existing Fortran/FPM, SPICE/DACE, Python standard library; no integrator or ADS-core rewrite.

## Global Constraints

- `s,u in [0,1]`; root coefficients `0.5*(1+q7/q8)`, patch nominal `0.5*(1+center7/8)`, span `0.25*width7/8`.
- Preserve existing 6D/7D APIs, old study points, runner CLI/schema, and v1 snapshots. Any old event serialization must still emit exactly seven score/center/width entries.
- No probability measure for s/u. Envelope diagnostics are sampled domain diagnostics, not claimed rigorous global extrema.
- Full-sun, EARTH/J2000/TDB, km/km/s/km/s2; use the immutable finalized P1 table and matching common Cr/area/mass/pressure.
- C2 uses eight active variables, DA order 4, depth 8, domain_sigma 3; tolerances 0.1 km and 1e-6 km/s, RKF78 1e-12/1e-12. Default duration 360 h, save 1 h; P3 actual execution is 1 h.
- Zero-mode C2 versus C0: position <=1e-6 km, velocity <=1e-9 km/s. ADS versus same-proxy Real: position <=0.1 km, velocity <=1e-6 km/s.
- All owned DA handles and eta/mode coefficient state are cleared on normal/status-error exits. No fallback to Cannonball on a modes error.
- Formal wrapper re-audits mode provenance, binds runtime OPM/config/mode CSV/meta/executable/source hashes, stages all files and atomically publishes only structurally valid successful results; never overwrite an existing output directory.
- P4's 15-day experiment, physical Box-Wing surrogate error and convergence study are not run or declared passed in P3.

## Task 1: HFEM 8D integration

**Files:** Modify `src/lib/uncertainty/propagation/pod_uq_hfem_ads_module.f90`; minimally slice legacy event writers in `app/run_srp_ads_history.f90`, `app/run_srp_ads_mc_pair.f90`, `app/run_srp_ads_mc_initial_only.f90`; create `test/test_hfem_ads_modes_8d.f90`.

**Interfaces:** `hfem_ads_options_type%parameterization=ADS_PARAM_AUTO`; stats record resolved parameterization. Existing propagate/history calls remain unchanged. Modes samples have physical rows `[state6,s,u]`; history evaluate still consumes unit points and returns state6.

- [x] Write tests: explicit modes init/evaluate, 8D one-shot output retains both parameters, short actual propagation versus Real, actual split directions 7/8 and 8-entry event records, invalid model/parameterization, lifecycle handle baseline.
- [x] Run `source setup_env.sh && fpm test test_hfem_ads_modes_8d`; verify missing parameterization fails compilation.
- [x] Pass options parameterization to both map constructors; allow da_storage_variables up to 8; expand shared event arrays to 8 and pass only active storage slices to split scoring. Preserve old serializers with `(1:7)` slices.

```fortran
select case(map%parameterization)
case(ADS_PARAM_NONE)
    call clear_srp_scale_uncertainty()
    call clear_srp_mode_uncertainty()
case(ADS_PARAM_SRP_ETA)
    call clear_srp_mode_uncertainty()
    call set_srp_scale_uncertainty(7,map%domain_sigma*map%srp_sigma*center(7), &
        map%domain_sigma*map%srp_sigma*width(7)/2)
case(ADS_PARAM_SRP_MODES)
    call clear_srp_scale_uncertainty()
    call ads_patch_mode_parameters(map,center(1:8),width(1:8),s0,ds,u0,du,status)
    call set_srp_mode_uncertainty(7,s0,ds,8,u0,du)
end select
```

- [x] Add model/table/time preflight and end-time checks without changing integrators; reject unsupported combinations before owning DA handles. Test new path plus existing ADS history/zero-duration/embedded-dimension targets.

## Task 2: Study points, Real probes and v2 snapshots

**Files:** Modify `pod_srp_ads_history_module.f90`, `pod_srp_ads_snapshot_module.f90`, `tools/reconstruct_srp_ads_snapshot.py`; add `test/test_srp_modes_study.f90`, `test/test_srp_modes_snapshot.f90`, `test/test_srp_modes_reconstruct.py`. Existing history runner passes explicit 7 to study builder.

**Interfaces:** `build_study_points(points,validation_ids,n_variables)` defaults to legacy 7 if omitted; supports 6/7/8. `advance_real_mode_probes(states,s_values,u_values,epoch0,t_start,t_end,options,status,message)`. Snapshot writer adds optional final `schema_version,mode_csv_sha256,mode_metadata_sha256` arguments; default v1. v2 map contains schema_version=2, n_variables, da_variables, parameterization (NONE/SRP_ETA/SRP_MODES), variable_role rows, mode CSV/meta SHA. v2 header uses actual storage dimensions, v1 remains seven padded columns.

- [x] Test exact old grid preservation and all 256 modes corners, 16 axes, center, 8-base Halton points. Validation includes those 273 boundary/center points plus 64 interior points (337 total), shape remains 512.
- [x] Test Real mode probe finite/domain/shape/time/coverage errors and clearing; write short actual-propagation checks, then implement explicit per-point setters and integration completion checks.

```fortran
call clear_srp_scale_error()
call clear_srp_mode_coefficients_real()
do j=1,size(s_values)
    call set_srp_mode_coefficients_real(s_values(j),u_values(j))
    ! Existing adaptive_step_integrate with nondimensional state/time.
end do
call clear_srp_mode_coefficients_real()
```

- [x] Test 8D power8 polynomial and split-history round trips; confirm new APIs/header fail before implementation. Generate dynamic CSV headers with internal integer formatting and explicit v2 identity rows; Python loads v1 inferred dimensions and v2 explicit dimensions, checks roles/parameterization and 64-hex modes hashes.
- [x] Run new snapshot tests and existing `test_srp_ads_snapshot_io` / `test_srp_ads_reconstruct.py`; no old-schema regression.

## Task 3: Independent C2 runner and publication wrapper

**Files:** Create `app/run_srp_modes_ads.f90`, `tools/run_srp_modes_ads.py`, `test/test_srp_modes_ads_runner.py`, `test/test_run_srp_modes_ads.py`; new complete `config/config_srp_modes_DRO-2.txt` matches P1 non-SRP settings and physical values, with mode paths.

**Interfaces:** Required CLI `-opm INPUT -cfg CONFIG -o PREFIX`; optional `--hours`, `--save-hours`, `--da-order`, `--max-depth`, `--pos-tol-km`, `--vel-tol-kms`. Default 360 h/1 h/4/8/0.1/1e-6. Wrapper arguments `--opm --config --output-dir` and same numeric options; executable selected by optional `--executable`, otherwise verified FPM target. Raw file prefix is `c2` inside staging/publication directory.

- [x] Write tests for help, finite CLI, wrong model, missing/mismatched tables, no overwrite, 1-hour output structure and zero-mode/ADS tolerances; missing target/output is a failure even if FPM returns 0.
- [x] Runner initializes engine, reads validated OPM, loads finalized table at OPM epoch, creates MODES history, maps 512 study points and advances 337 independent Real probes. Store all states and s/u, sampled envelope and 8 split counts/events at each checkpoint. Compare zero-mode nominal-state probe with C0 under identical forces/tolerances; do not mutate mode table or geometry for control.
- [x] Emit raw JSON `schema_version=srp_modes_ads_v1`, `formal_result=false`, parameterization/roles, input paths, OPM epoch, settings, numerical maxima, zero-mode maxima, depth-limited/precision flags and snapshot prefixes. Save v2 snapshots at initial/final and event checkpoints, including previous domain when a split occurs.
- [x] Wrapper re-audits P1, compares runtime OPM hash and common physical/non-SRP parameters, validates each CSV record and snapshot, records executable and source hashes, rewrites only its owned staging report and atomically publishes a new directory. Failures preserve staging.

```python
audit = audit_table(csv_file, meta_file)
assert audit['audit_complete'] and audit['gate_passed']
assert file_sha256(opm_file) == audit['opm_sha256']
# Run the bound executable into a new tempfile.mkdtemp beside output-dir.
# Audit reports, 8D points/events, snapshots, and threshold/precision flags.
# Recheck input/executable hashes, add provenance, then staging.rename(output_dir).
```

- [x] Run native one-hour C2 smoke via wrapper and independently reconstruct initial/final snapshots. Verify numerical and zero-mode limits and handle cleanup; preserve failed staging rather than claiming completion.

## Task 4: Review, focused regression and report

- [x] Run P1/P2 tests, new HFEM/study/snapshot/runner/wrapper tests, old snapshot/runner focused Python tests and existing 6D/7D focused Fortran tests. No full long-arc FPM suite in this stage.
- [x] Inspect complete diff, request independent code review if available, fix material findings and rerun covering tests. Check `git diff --check`.
- [x] Record smoke measurements, publication provenance and reconstructed snapshot errors in `docs/superpowers/reports/2026-10-14-boxwing-modes-p3-results.md`; update design phase status without claiming P4 or merge readiness.

Execution uses the existing isolated worktree. Preserve P1/P2 data, runtime symlinks, prior build archive and unrelated main-checkout campaign files. Tests and FPM builds are serialized to avoid competing writes to `.mod` files.

## Completion checkpoint — 2026-10-14

21 Fortran focused targets and 40 Python tests passed, including 16 isolated native error paths. Independent read-only review found and re-reviewed three corrected contracts: DA completion-array orientation, default raw metadata hash, and 6/7-active-in-8-storage reconstructed CSV padding. Many-step and incomplete-integration ownership regressions passed. Final audited one-hour publication is `SRP/boxwing_modes_DRO-2_init_p3_smoke_20261014_verified`; the earlier smoke remains preserved with its own historical provenance. Initial/final snapshots independently reconstructed and all published artifact hashes verified. No P4, complete long-arc suite, commit or merge performed.
