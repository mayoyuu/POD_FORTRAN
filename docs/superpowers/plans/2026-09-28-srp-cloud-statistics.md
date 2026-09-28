# SRP Cloud Statistics Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Extend the portable 15-orbit SRP paired-analysis script with km/km/s outputs, fixed-design cloud moments, first-split comparison, relative maximum error figures, and an evidence-based audit of whether saved ADS data support hourly Gaussian resampling.

**Architecture:** Keep `tools/analyze_srp_non_nrho15.py` as the portable entry point. Add pure numerical helpers for moment and split calculations, integrate their results into per-orbit CSVs, and isolate plotting from calculation. Audit scheme-two inputs with a separate read-only function that inventories saved snapshots and writes a machine-readable and human-readable conclusion.

**Tech Stack:** Python 3.10+, NumPy, pandas, matplotlib, unittest.

## Global Constraints

- Include only the explicit 15-orbit whitelist and exclude `L2NRHO-1/2`.
- Position units are km and velocity units are km/s.
- Fixed 512-point skewness and kurtosis are geometric indicators, not probability normality tests.
- Exclude the first failed/unconfirmed checkpoint and every later checkpoint from statistics and figures.
- Retain raw invalid values only in the audit history with explicit exclusion columns.
- First split is resolved to the saved one-hour interval from cumulative split counters.

---

### Task 1: Cloud moment calculations and unit conversion

**Files:**
- Modify: `tools/analyze_srp_non_nrho15.py`
- Create: `test/test_analyze_srp_non_nrho15.py`

**Interfaces:**
- Produces: `component_shape_moments(values: np.ndarray) -> tuple[np.ndarray, np.ndarray]`
- Produces: `cloud_metrics(group: pd.DataFrame, suffix: str) -> dict[str, object]`
- Produces columns named `delta_r_*_km`, `delta_v_*_km_s`, `*_mean_origin_*`, `*_longest_rms_axis_*`, `*_max_abs_skewness`, `*_max_abs_excess_kurtosis`, and `*_relative_max_*`.

- [ ] **Step 1: Write failing moment tests**

Use symmetric arrays to assert zero skewness, Gaussian-like seeded arrays to assert near-zero excess kurtosis, asymmetric arrays to assert the correct maximum component, and zero-variance arrays to assert missing values rather than division by zero.

- [ ] **Step 2: Run the focused test and verify failure**

Run `python -m unittest test.test_analyze_srp_non_nrho15 -v`; expect missing helper failures.

- [ ] **Step 3: Implement moment helpers and km/km/s fields**

Use population central moments (`ddof=0`). Compute position and velocity centroid norms relative to the origin, separate longest covariance RMS axes, all six signed skewness/excess-kurtosis values, maximum absolute values plus component labels, and maximum norm divided by the corresponding longest RMS axis.

- [ ] **Step 4: Run focused tests and verify pass**

Expect all numerical tests to pass with finite tolerances.

### Task 2: First split comparison and plotting

**Files:**
- Modify: `tools/analyze_srp_non_nrho15.py`
- Modify: `test/test_analyze_srp_non_nrho15.py`

**Interfaces:**
- Produces: `first_split_record(history: pd.DataFrame, prefix: str, dimensions: Iterable[int]) -> dict[str, object]`
- Produces: `first_split_comparison.csv`
- Produces: `figures/first_split_comparison.png`
- Produces per orbit `cloud_moments` and `relative_max_error` PNGs.

- [ ] **Step 1: Write failing split tests**

Create cumulative counters with first nonzero values at known hours, multiple first-hour directions, and no-split data. Assert interval start/end, direction labels, null handling, and control dimension-seven zero validation.

- [ ] **Step 2: Run the focused test and verify failure**

Expect missing split-record helper failure.

- [ ] **Step 3: Implement comparison CSV and plots**

Create one row per orbit with treatment/control first saved hour, interval, direction set, and hour difference. Plot valid treatment/control cloud metrics using shared axes; stop lines at the first invalid checkpoint. Mark Gaussian reference zero lines on skewness and excess-kurtosis panels.

- [ ] **Step 4: Run focused tests and verify pass**

Expect moment and split tests to pass.

### Task 3: Campaign integration and real-data verification

**Files:**
- Modify: `tools/analyze_srp_non_nrho15.py`
- Update generated artifact: `C:/Users/13113/Desktop/科研/误差传播/SRP影响/analyze_srp_non_nrho15.py`
- Regenerate: `C:/Users/13113/Desktop/科研/误差传播/SRP影响/analysis_non_nrho15/`

**Interfaces:**
- Produces `cloud_statistics/<orbit>_cloud_statistics.csv`, `cloud_orbit_summary.csv`, updated population history, per-orbit figures, and updated manifest.

- [ ] **Step 1: Integrate helpers into per-orbit processing**

Write all 361 raw rows with validity flags and a valid-only cloud-statistics view. Rename all m/mm/s derived columns and plot labels to km/km/s.

- [ ] **Step 2: Run the 15-orbit campaign to a new build directory**

Expect 15/15 completed, 361 checkpoints per raw history, 512 points per checkpoint, and no L2NRHO paths.

- [ ] **Step 3: Verify scientific and file invariants**

Assert `L1PL-2` first invalid hour 259, all three DRO tracks valid through 360 h, figures use only valid rows, no included statistic exceeds its source valid maximum, and all requested PNG/CSV files exist and parse.

- [ ] **Step 4: Atomically replace the current analysis directory**

Preserve the prior directory as a timestamped backup; verify both paths resolve under the intended `SRP影响` directory before moving.

### Task 4: Scheme-two hourly Gaussian sample feasibility audit

**Files:**
- Modify: `tools/analyze_srp_non_nrho15.py`
- Create generated outputs: `analysis_non_nrho15/gaussian_sampling_feasibility.json` and `.md`

**Interfaces:**
- Produces: `audit_gaussian_sampling_feasibility(input_root: Path) -> dict`

- [ ] **Step 1: Inventory real saved artifacts**

Inspect initial point weights/design, snapshot hours, patch-domain maps, sparse DA coefficients, reconstruction metadata, and validity coverage separately for treatment and control.

- [ ] **Step 2: Classify feasibility from evidence**

Return exactly one of `hourly_direct`, `snapshot_hours_only`, or `requires_new_propagation_or_hourly_snapshots`, with per-orbit evidence and missing fields.

- [ ] **Step 3: Write JSON and Markdown audit outputs**

State whether arbitrary Gaussian samples can be assigned to patches and evaluated at each hour, what probability assumptions are missing, and what minimal new output would make it possible.

- [ ] **Step 4: Verify audit claims against files**

Check that every claimed available snapshot exists and every claimed missing hourly series was searched for; do not infer availability from filenames alone.

### Task 5: Final verification and commit

**Files:**
- Modify: `tools/analyze_srp_non_nrho15.py`
- Create: `test/test_analyze_srp_non_nrho15.py`
- Modify: `tools/requirements-srp-analysis.txt` only if dependencies change.

- [ ] **Step 1: Run syntax and unit tests**

Run `python -m py_compile tools/analyze_srp_non_nrho15.py` and the unittest module; expect success.

- [ ] **Step 2: Inspect representative images and CSVs**

Visually inspect DRO-1, L1PL-2, first split, and population figures; verify labels, units, truncation, and legends.

- [ ] **Step 3: Commit implementation**

Add only the analysis script, requirements, and tests; exclude generated data and `__pycache__`. Commit with a focused message.
