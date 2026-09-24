# SRP ADS Snapshot Reconstruction Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a standard-library Python CLI that densely reconstructs saved 6D or 7D ADS patch snapshots and verifies them against the runner's saved point cloud.

**Architecture:** A focused module parses and validates the four snapshot CSVs into patch objects, generates deterministic local samples, evaluates sparse polynomials, and streams reconstructed rows. Optional parent-run discovery supplies nominal states and a 512-point cross-check.

**Tech Stack:** Python 3 standard library, `unittest`, existing Fortran snapshot CSV schema.

## Global Constraints

- Do not change treatment/control propagation or the 512-point hourly design.
- Support both active 7D treatment and padded 6D eta-zero snapshots.
- Default to 4096 deterministic samples per patch with 50% of extra points on patch faces.
- Treat reconstructed output as the saved ADS polynomial approximation, not a physical continuous-domain guarantee.
- Fail saved-cloud verification above `1e-8 km` position or `1e-11 km/s` velocity.
- Do not add NumPy or another runtime dependency.

---

### Task 1: Snapshot parser and polynomial evaluator

**Files:**
- Create: `tools/reconstruct_srp_ads_snapshot.py`
- Create: `test/test_srp_ads_reconstruct.py`

**Interfaces:**
- Produces: `load_snapshot(prefix: Path) -> Snapshot`
- Produces: `Snapshot.evaluate_global(point: tuple[float, ...]) -> tuple[int, tuple[float, ...]]`
- Produces: `replay_split_history(entries: list[int], n_variables: int) -> tuple[list[float], list[float]]`

- [ ] **Step 1: Write failing synthetic parser/evaluator tests**

Create a 7D unsplit fixture whose component 1 is
`10 + 2*v1 - 3*v2*v7` and whose other five components are constants. Assert
evaluation at a known global point and assert metadata/map fields.

- [ ] **Step 2: Run the test and confirm RED**

Run:

```sh
python3 test/test_srp_ads_reconstruct.py
```

Expected: import failure because `tools/reconstruct_srp_ads_snapshot.py` does not exist.

- [ ] **Step 3: Implement minimal parser and evaluator**

Parse map, metadata, split rows, and coefficients. Validate active dimension
from zero/nonzero seventh widths and eta sigma. Use first-containing-patch
selection and `local=2*(global-center)/width`.

- [ ] **Step 4: Run the test and confirm GREEN**

Run the same command. Expected: parser/evaluator test passes.

### Task 2: Deterministic per-patch sampling and CSV output

**Files:**
- Modify: `tools/reconstruct_srp_ads_snapshot.py`
- Modify: `test/test_srp_ads_reconstruct.py`

**Interfaces:**
- Produces: `generate_local_samples(n_variables: int, count: int, boundary_fraction: float) -> list[Sample]`
- Produces: `reconstruct_snapshot(...) -> dict`

- [ ] **Step 1: Write failing sampling tests**

Assert exact sample count, all corners/axes/center, deterministic repeated
output, active coordinates within `[-1,1]`, padded seventh coordinate zero
for 6D, and both `face` and `interior` samples for a sufficiently large
count.

- [ ] **Step 2: Run and confirm the expected missing-function failure**

Run the focused unittest and verify RED.

- [ ] **Step 3: Implement sample generation and streaming CSV writer**

Use prime-base Halton coordinates. Map local to global coordinates, compute
physical eta, evaluate six sparse polynomials, and write the documented stable
columns. Track extrema without retaining the full reconstructed cloud.

- [ ] **Step 4: Run and confirm GREEN**

Verify deterministic sampling and CSV schema tests pass.

### Task 3: Structural and saved-cloud verification

**Files:**
- Modify: `tools/reconstruct_srp_ads_snapshot.py`
- Modify: `test/test_srp_ads_reconstruct.py`

**Interfaces:**
- Produces: `discover_parent_files(snapshot_prefix: Path) -> ParentFiles | None`
- Produces: `verify_saved_cloud(snapshot: Snapshot, parent: ParentFiles, ...) -> dict`

- [ ] **Step 1: Write failing split and saved-cloud tests**

Create two 1D-split-in-7D patches with histories `-1` and `1`. Assert
replayed centers `-0.5/+0.5`, widths 1, correct boundary ownership, and
near-zero comparison against synthetic initial/history/shape CSVs.

- [ ] **Step 2: Write malformed-fixture tests**

Assert explicit failures for metadata/split disagreement, unknown coefficient
patch, missing component, nonzero power7 in a 6D snapshot, and verification
above tolerance.

- [ ] **Step 3: Run and confirm RED**

Expected failures identify the missing verification and validation behavior.

- [ ] **Step 4: Implement structural checks and saved-cloud verification**

Replay every split history, compare centers/widths, locate the snapshot history
row, evaluate all saved global points, and compare to nominal plus deviations.

- [ ] **Step 5: Run and confirm GREEN**

All synthetic reconstruction tests pass without warnings.

### Task 4: CLI, real snapshot integration, and documentation

**Files:**
- Modify: `tools/reconstruct_srp_ads_snapshot.py`
- Modify: `SRP_ADS_HISTORY.md`
- Modify: `test/test_srp_ads_reconstruct.py`

**Interfaces:**
- CLI: `--prefix`, `--points-per-patch`, `--boundary-fraction`,
  `--output`, `--report`, `--history`, `--absolute-only`, and
  `--no-verify`.

- [ ] **Step 1: Write a failing CLI test**

Invoke the tool on a synthetic fixture and assert output/report paths, row
count, schema version, active dimension, and verification status.

- [ ] **Step 2: Implement argument validation and JSON report**

Reject invalid fractions/counts and missing history unless
`--absolute-only`. Write source SHA-256 hashes, sampling settings, extrema,
and verification results atomically.

- [ ] **Step 3: Add usage documentation**

Document treatment and eta-zero examples, output columns, accuracy semantics,
and the fact that only saved snapshot hours can be reconstructed.

- [ ] **Step 4: Run a real short-run integration**

Reconstruct a runner-generated start or end snapshot with at least 512 points
per patch. Expected maximum differences are below `1e-8 km` and
`1e-11 km/s`.

- [ ] **Step 5: Run the complete verification gate**

```sh
python3 -m py_compile tools/reconstruct_srp_ads_snapshot.py
python3 test/test_srp_ads_reconstruct.py
python3 test/test_srp_ads_batch.py
python3 test/test_srp_ads_pair_compare.py
git diff --check
```

Expected: all commands exit zero.
