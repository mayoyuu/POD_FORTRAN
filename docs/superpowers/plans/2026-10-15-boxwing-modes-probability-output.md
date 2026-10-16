# Box-Wing Modes Probability Output Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox syntax for tracking.

**Goal:** Extend the existing C2 ADS runner with audited hourly 50,000-sample moments, fixed 2,000-point visualization clouds, full keyframes, and hourly reconstructable snapshots.

**Architecture:** Opt-in probability output consumes fixed global 8D coordinates and evaluates the same live ADS history. A dedicated Fortran module owns strict input, full-population statistics and layered CSV output; Python generates the population, audits/streams it to a new HDF5 schema, and publishes using the existing C2 transaction.

**Tech Stack:** Fortran/FPM/DACE, Python 3 stdlib, isolated NumPy/h5py environment.

## Global Constraints

- Approved spec: main repository `docs/superpowers/specs/2026-10-15-boxwing-modes-probability-output-design.md` (user approved 2026-10-15).
- Work only in `/home/songyu/fortran_test/POD_Fortran/.worktrees/boxwing-modes-p1`; preserve P1–P4 changes. On 2026-10-16 the user explicitly authorized a local commit followed by the full local probability run. This supersedes the development-phase no-commit restriction only; no new branch, merge, push, C0/C1 rerun, or historical output mutation is authorized.
- Independent componentwise ±3 truncated standard normals for z1:6; q=z/3; independent uniform s/u; q7:8=2[s,u]-1; xi_E=s, xi_M=(1-s)u. PCG64DXSM seed 20261015, two spawned streams, stable prefixes.
- N=50000 computes statistics independently of V=min(2000,N). Fixed IDs use 1+floor(j*(N-1)/(V-1)+.5). All 36 covariance elements, ddof=1, linear quantiles, existing corrected skew/kurtosis definitions.
- Full keyframes at 0/24/72/168/240/360 h, actual endpoints, and first split/depth/failure with preceding hour. No hourly full-population CSV staging.
- Retain original 512/337 validation and zero gates; separate min(64,N) Real probability probes; hourly full v2 snapshots. Same physics and integrator, no relaxed tolerances.
- New isolated results only. Formal acceptance separates numerical validation and 20k→50k statistical convergence. Small runs cannot certify 50k sufficiency.
- Original CLI remains stdlib-only. New Python dependencies load only for opt-in probability output. Never install into system Python.

## File map and interfaces

- `tools/generate_srp_modes_probability_samples.py`: `generate_samples(count,seed)` → (N,19) canonical numeric rows + manifest; `write_sample_bundle(csv,manifest,...)`; `load_sample_bundle(csv,manifest,expected_count)` strict validated population; `choose_ids(N,V)`.
- `src/lib/uncertainty/io/pod_srp_modes_probability_io_module.f90`: `probability_population_type` and `probability_output_type`; `load_probability_samples`, `probability_open`, `probability_write_hour`, `probability_close`; reuses `compute_group_stats`, never edits integrator.
- `tools/pack_srp_modes_probability_hdf5.py`: strict staged frame reader, full-keyframe independent statistics audits, HDF5 pack/audit, convergence, probe/snapshot audits, disk preflight.
- `app/run_srp_modes_ads.f90`: opt-in parsing, initial probe mapping, probability Real advancement, all-N ADS evaluation and patch lookup, hourly output/snapshot and raw report counts.
- `tools/run_srp_modes_ads.py`: opt-in arguments/preflight, frozen samples/manifest, provenance, probability audit/pack before atomic publication.
- `tools/requirements-srp-modes-probability.txt`: NumPy/h5py only.
- `test/test_srp_modes_probability_samples.py`, `test/test_srp_modes_probability_io.f90`, `test/test_srp_modes_probability_output.py`, `test/test_srp_modes_probability_runner.py`: focused red/green tests including negative paths.

### Task 1: Fixed, validated 8D prior population

- [x] Write failing tests for deterministic prefix stability, no clipping, known truncated/uniform moments, exact mappings, fixed IDs, manifest/hash/schema/duplicate/nonfinite/domain rejection and no overwrite.

```python
small,_=generate_samples(73,20261015)
large,_=generate_samples(50000,20261015)
np.testing.assert_array_equal(small,large[:73])
assert choose_ids(50000,2000)[[0,-1]].tolist()==[1,50000]
```

- [x] Run `.venv-probability/bin/python -m unittest discover -s test -p test_srp_modes_probability_samples.py -v`; expect missing generator assertion.
- [x] Implement generator, strict loader and requirements; repeat command until green. Manifest records complete prior, seed/spawn keys/rejection count/version/hash/ID digest and prefixes.

### Task 2: Native complete moments and layered output

- [x] Add synthetic native test with correlated position/velocity components, constant component flags, strict sample loader negatives and event keyframe selection.

```fortran
call probability_open(output,prefix,population,2,status)
call probability_write_hour(output,population,0.0_DP,epoch,states,nominal,patch_ids, &
    .true.,.true.,.false.,.false.,.false.,1,0,0.0_DP,0.0_DP,0.0_DP,0.0_DP,.false.,status)
```

- [x] Run `source setup_env.sh && fpm test --target test_srp_modes_probability_io`; expect missing new module.
- [x] Implement strict 19-field CSV reader, prefix moments writer, fixed V output, cached prior frame for event keyframes, close/final frame handling; check finite values and all IDs before output.
- [x] Rerun native test; compare exported covariance/mean/quantiles against NumPy in Python output tests. In particular require nonzero `cov_14` and symmetric `cov_41`.

### Task 3: Strict HDF5 packer and statistical sufficiency

- [x] Write failing fixture tests: all hourly V vs keyframe N shapes/weights, all 36 covariances, matching shared IDs, keyframe moment recomputation, missing/duplicate/nonfinite frames, false count/reference/covariance, degenerate moments and convergence failure.

```python
with h5py.File(path) as data:
    assert data['visualization/state'].shape==(3,2,6)
    assert data['keyframes/state'].shape==(2,8,6)
    np.testing.assert_allclose(data['population/weight'][:].sum(),1)
```

- [x] Run `.venv-probability/bin/python -m unittest discover -s test -p test_srp_modes_probability_output.py -v`; expect absent packer assertion.
- [x] Implement two-pass bounded-memory validation/packing, attributes, gzip4/shuffle/chunks, full-keyframe NumPy statistics audit and dimensionless PSD checks; implement block-scaled 20k→50k convergence with explicit degenerate/smoke states.
- [x] Rerun tests and reject each corruption fixture without publishing/deleting input CSVs.

### Task 4: Opt-in C2 native/wrapper integration

- [x] Add failing native/wrapper tests requiring paired sample+manifest arguments, count validation, unchanged default source behavior, preserved original counts, separate probability probes and all-hour snapshots.

```python
report=run_modes(opm,cfg,out,hours=2,probability_samples=samples,
    probability_manifest=manifest,probability_count=50000,
    probability_visualization_count=2000,probability_validation_count=64)
assert report['shape_count']==512 and report['validation_count']==337
assert report['probability_outputs']['statistics_sample_count']==50000
```

- [x] Run `.venv-probability/bin/python -m unittest discover -s test -p test_srp_modes_probability_runner.py -v`; expect unsupported options.
- [x] Extend native opt-in path only, evaluate all N with deterministic `mf_find_patch`, integrate independent probability probes, save all-hour snapshots and preserve original report contract.
- [x] Extend wrapper strict preflight/dependency/disk checks, freeze inputs, bind new source hashes, run/audit probes and 64 live/snapshot points, pack HDF5/convergence, recheck hashes, then publish. Owned staging point CSVs removable only after verified HDF5; keep on error.
- [x] Rerun new tests plus existing `test_run_srp_modes_ads.py`, `test_srp_modes_ads_runner.py`, modes snapshot/study/8D split tests.

### Task 5: Actual 50k/2k short propagation and handoff

- [x] Generate canonical 50k population in a fresh owned input directory. Run two-hour probability case into a new exclusive result directory (no full 15-day run before smoke passes).
- [x] Verify 3 hourly V=2000 frames, 2 endpoint N=50000 frames unless real events add frames, 9 stats rows, full covariance, separate 64 probes/hour, 3 snapshots, zero handle leak and all original gates.
- [x] Reuse the already published identical-configuration 15-day baseline's first 0/1/2 h instead of rerunning that baseline; compare nominal/512 states with required ≤1e-6 km/1e-9 km/s limits. Also run a small-N opt-out case and repeat probability output with another V; verify statistics unchanged.
- [x] Re-run targeted tests and `git diff --check`; record exact commands/results and remaining full-case acceptance work. Do not claim 360 h complete from smoke.

## Plan self-review / execution choice

Spec coverage: Tasks 1–4 implement §§2–11, Task 5 verifies short-case acceptance in §12. Full 360 h run/convergence is a later numerical acceptance after the short gate, not something implied by writing code. File interfaces and schema names follow the approved spec. No placeholder code steps or altered math/physics are introduced. The user explicitly requested implementation in this conversation: execute inline with TDD and checkpoints, without starting new agents or asking for a redundant workflow choice.

## Execution record (2026-10-16)

Tasks 1–5 are implemented and verified. Smoke published to `/home/songyu/fortran_test/POD_Fortran/SRP/boxwing_modes_DRO-2_init_c2_probability_smoke_20261016`: 3 hourly subset frames, 2 full endpoint frames, 9 statistics rows, and 192 live/snapshot probability probe checks. The existing 15-day baseline's first 2 h nominal/512 differences are exactly zero. Read-only review findings were fixed with red/green regressions and re-reviewed. Final targeted regression: 25 Python tests and 7 native targets passed; `git diff --check` passed. The full 360 h probability campaign is not started. Detailed evidence and command: main repo `docs/superpowers/reports/2026-10-16-boxwing-modes-probability-implementation-smoke.md`.

## Authorized next execution (2026-10-16)

The user requested “先提交再启动”: verify and commit the current related P1–P4/probability implementation and documentation, then launch one new local 360 h DRO-2_init case with N=50000, V=2000 and 64 independent probability probes. Reuse existing input samples and mode table; do not push or overwrite old results. The source revision in the new run must be the new commit. The actual launch/run status belongs in runtime logs and the resulting C2 report, not a claim that 360 h acceptance is already complete.
