# Box-Wing Modes P1 Implementation Plan

> **For agentic workers:** Use superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox syntax for tracking.

**Goal:** Produce audited Earth/Moon SRP error modes evaluated on one Cannonball reference trajectory for DRO-2 over 15 days, and report the 0.3 norm gate.

**Architecture:** An independent Fortran executable propagates only Cannonball and evaluates all three SRP laws on each saved state. A standard-library Python auditor recomputes modes, validates identities and sampling, uses the executable's SPICE epoch inspector, and publishes a complete output directory with a frozen executable by rename. A failed physical gate remains a valid diagnostic dataset but cannot be finalized for C2.

**Tech Stack:** Fortran/FPM, existing RKF78/SPICE/frame/force modules, Python 3 standard library and unittest.

## Global Constraints

- Reference frame J2000, center EARTH, force epochs absolute TDB seconds; km, km/s and km/s2.
- Full-sun, no random errors, no DA state, no changes to shared force/integrator/ADS modules.
- RKF78 relative/absolute tolerances 1e-12; minimum step 1e-6 s; maximum step 3600 s.
- Default duration 360 h and sampling 1 h (361 rows); allow 0.5 h sampling (721 rows).
- Both maximum mode norms must be <=0.3 to finalize. Failed gates preserve staging diagnostics.
- Existing published files must never be overwritten implicitly.

## Task 1: Auditor contract and negative cases

**Files:** Create `test/test_audit_boxwing_mode_table.py` and `tools/audit_boxwing_mode_table.py`.

**Interfaces:** `audit_table(csv_path, meta_path, expected_hours=360.0) -> dict`; `finalize_table(csv_path, meta_path, output_dir, expected_hours=360.0) -> dict`. CLI: `--csv FILE --meta FILE [--finalize --output-dir DIR] [--expected-hours HOURS]`.

- [x] Write synthetic, algebraically known vector fixtures and test recomputed norms/maxima, 0.5 h cadence, corrupt header/time/mode/norm, nonfinite data, input and finalized CSV hash tampering, metadata identity, RTN singularity, gate failure, existing destination rejection, and no partial publication. Add whole-table epoch shift and frozen executable identity regressions.
- [x] Run `python3 -m unittest discover -s test -p test_audit_boxwing_mode_table.py`; confirm the missing auditor produces the expected failure.
- [x] Implement strict schema `boxwing_srp_modes_v1`, identity/sampling/numeric/RTN checks, hashes, physical gate and atomic directory publication.
- [x] Re-run the same tests; all 13 pass. The extractor must be built first for the authoritative SPICE conversion used by these fixtures.

## Task 2: Extractor contract and integration test

**Files:** Create `test/test_boxwing_mode_extract.py` and `app/run_boxwing_mode_extract.f90`.

**Interfaces:** `fpm run run_boxwing_mode_extract -- -opm FILE -cfg FILE -o PREFIX [--hours H --sample-hours H]`; writes `PREFIX.csv` and `PREFIX_meta.json`, with `audit_complete=false`. `--utc-to-tdb UTC_STRING` returns a SPICE-converted epoch for independent identity checks.

- [x] Write CLI and 1 h / 0.5 h smoke tests. Check shared states/epochs, independent historical Cannonball acceleration and reference trajectory (accounting for the recorded pressure rounding), exact mode algebra and RTN invariance via auditor; require no overwrite and reject invalid cadence.
- [x] Run `python3 -m unittest discover -s test -p test_boxwing_mode_extract.py`; confirm target missing is the expected failure.
- [x] Implement one Cannonball trajectory with production Real SRP evaluation, restoring Cannonball before every integration segment. Record effective physical parameters, integrator settings, epoch, fallback diagnostics and sampling metadata.
- [x] Build and re-run all 4 integration tests; run existing Box-wing Real/config/geometry focused regression tests.

## Task 3: Full P1 dataset and physical gate

**Files:** Generate `SRP/boxwing_modes_DRO-2_init/` only after the gate passes; retain failed candidates in a uniquely named staging directory.

- [x] Build the extractor, run 360 h / 1 h against the repository DRO-2 OPM and complete `config/config_box_wing_earth.txt` geometry config, force Cannonball for propagation, and retain the complete log. `config/config.txt` contains sentinel Box-wing geometry and is not usable for extraction.
- [x] Audit all 361 rows, inputs, reference identity, finite values, times, vector algebra and RTN invariance.
- [x] Report both maxima and epochs. Finalize only when both are <=0.3; both passed.
- [x] Run a separate 0.5 h extraction; both maxima remain stable. Preserve and archive earlier audit versions without losing raw data.

## Task 4: Verification and handoff

**Files:** Add a P1 usage/result note in `docs/superpowers/`; update the sampling contract in the existing design document.

- [x] Re-run both new suites and focused existing tests; inspect `git diff --check` and the exact diff. Independent review identified the epoch identity gap; its regression was confirmed failing and then passing after the SPICE fix.
- [x] Independently verify finalized hashes, frozen executable identity, SPICE epoch, source-file hashes, baseline commit and 361/721 data row counts; both physical gates passed.
- [x] Deliver source files, dataset/report paths, maximum norms, gate result and verification scope. Do not start P2.

No commit, merge or cleanup is needed to run P1; leave the isolated checkout available for review.
