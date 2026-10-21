# DRO-2_init P4 Complete Case Implementation Plan

> **For agentic workers:** Use superpowers:subagent-driven-development task-by-task, TDD and verification-before-completion. The existing design is approved and the user explicitly asks to continue until one complete orbit case is verified. Preserve uncommitted P1/P2/P3 and completed results; do not commit, merge, push, or rerun C0/C1 domains.

**Goal:** Finish a reproducible 360-hour DRO-2_init case with a validated C2 baseline, three convergence variants, independent Box-Wing proxy-error assessment, reused control/subdomain evidence, full relevant regression and a classified final report.

**Architecture:** Extend only the existing C2 launcher tolerance interface; keep dynamics/integrators/ADS algorithms unchanged. Reuse the immutable completed baseline and existing deterministic truth when its Real trajectory and identity are provable. A separate analysis tool checks each published run, aligns hourly states, quantifies convergence of sampled envelopes and physical SRP effects, and separates ADS numerical error from proxy-model error. Run minimal independent Real endpoint trajectories only when existing truth is not valid for reuse; never silently replace or relabel controls.

**Tech Stack:** Existing Fortran/FPM, DACE/SPICE, Python standard library; WSL Linux processes and frozen executables.

**Verified status 2026-10-15:** The complete 360-hour orbit case, all three convergence variants, independent Real truth, historical-control restriction, numerical gates and immutable snapshot re-evaluation are complete. Available Python modules: 67 tests passed. Full native suite: 101 attempted, 91 passed, 7 timed out at 300 s, 3 missing fixtures; 6 Python modules lack numpy/pandas/h5py. Independent child review unavailable due quota/capacity; parent review only. Therefore the complete orbit case is achieved, but full-repository green regression / overall design §17 remains open. Final main-repository report: `docs/superpowers/reports/2026-10-15-boxwing-modes-p4-complete-case-results.md`. The detailed durable checklist is `.superpowers/sdd/progress.md`; unchecked full-regression and independent-review requirements below must not be read as passed.

## Global Constraints

- Worktree `/home/songyu/fortran_test/POD_Fortran/.worktrees/boxwing-modes-p1`, branch `codex/boxwing-modes-p1`; preserve all prior dirty changes and runtime links.
- Baseline `/home/songyu/fortran_test/POD_Fortran/SRP/boxwing_modes_DRO-2_init_c2_15d_20261014` is already formal, 360 h, 361 hourly rows, DA order 4/depth 8, domain_sigma 3, RKF78 rel/abs tolerance 1e-12. Do not overwrite or rerun this baseline.
- C2 ADS gates remain position <=0.1 km, velocity <=1e-6 km/s; zero-mode gates remain 1e-6 km and 1e-9 km/s. No relaxation to obtain a pass.
- Mode sampling variant uses existing immutable 0.5 h/721-node P1 data; tighter integration uses rel_tol=abs_tol=1e-13; higher DA uses order 5. Each variant runs the full 360 h on the same 512 shape/337 Real validation grid, with one-hour output.
- Convergence key metrics are the hourly sampled uncertainty-envelope position/velocity maxima and nominal-state physical Earth/Moon SRP offsets where available, not percentage changes of near-zero numerical residuals. Compare max absolute metric difference across hours to max absolute baseline metric; <=5% means confirmed, otherwise explicitly unconfirmed. Every variant must independently meet numerical gates.
- `ads_numerical_error`, `surrogate_model_error`, and `conditional_srp_offset` are distinct. Physical full Box-Wing truth must be independent Real nominal propagation, not MC mean or DA constant without an independent consistency proof. A large proxy error is a model limitation to report, not a reason to alter the physical formula or numerical gates.
- Reuse C0 and C1 result artifacts. Historical eta±0.9 is C1-H; an honest ±0.3 restriction of saved polynomials may support final-time C1-30, but cannot invent unavailable hourly polynomials or claim a rerun of the smaller domain.
- No new 6D/7D control-domain propagation. Small Real Deta/Earth/Moon/proxy endpoint trajectories are in scope if no reusable independent truth exists.
- Each run uses new exclusive output/staging and frozen runner plus verified inputs/source hashes. No source edits while a run whose source identity includes those files is active. Serialize FPM builds/tests; do not alter unrelated active campaigns.
- Final deliverable reports actual 360-hour coverage, gates, handle cleanup, precision/depth flags, convergence, proxy error, snapshots, input/executable/source identities and control-reuse limitations. Do not call the whole P4 complete if a required experiment is missing.

## Task 1: Expose and bind C2 integrator tolerances

**Files:** Modify only `app/run_srp_modes_ads.f90`, `tools/run_srp_modes_ads.py`, `test/test_srp_modes_ads_runner.py`, `test/test_run_srp_modes_ads.py`.

**Interfaces:** Native optional `--rel-tol` and `--abs-tol`, default 1e-12, feed existing `options%rel_tol/abs_tol`. Python `run_modes(..., rel_tol=1e-12, abs_tol=1e-12)` and corresponding argparse flags. The existing raw report already has `rel_tol/abs_tol`; formal wrapper must verify them against requested values with no fixed absolute tolerance that would hide a 1e-13 identity mismatch. Existing callers remain unchanged.

- [ ] Add failing native help/invalid finite-positive tolerance tests and short run proving explicit 1e-13 values reach both the native report and frozen formal result.
- [ ] Add wrapper malformed tolerance and tiny-value report-identity regression; use real validation logic rather than assertions about mocks.
- [ ] Observe RED before implementation; implement only parsing, finite/positive validation, argv propagation and strict settings binding.

Native modifications:

```fortran
! Append '--rel-tol','--abs-tol' to the existing numeric option cases.
case('--rel-tol'); read(arg,*,iostat=status) options%rel_tol
case('--abs-tol'); read(arg,*,iostat=status) options%abs_tol
! Include both options in the existing ieee_is_finite and positive-limit guards.
```

Python modifications:

```python
# Add positive finite rel_tol/abs_tol to run_modes validation and native argv.
# Require relative-only agreement for requested tolerance identity.
require(math.isclose(float(report['rel_tol']), rel_tol, rel_tol=1e-12, abs_tol=0.0), 'run setting rel_tol mismatch')
require(math.isclose(float(report['abs_tol']), abs_tol, rel_tol=1e-12, abs_tol=0.0), 'run setting abs_tol mismatch')
```

- [ ] Run both focused runner modules; exact commands `source setup_env.sh && python3 -m unittest discover -s test -p test_srp_modes_ads_runner.py -v` and `python3 -m unittest discover -s test -p test_run_srp_modes_ads.py -v`; default and explicit tolerances pass.
- [ ] Task-scoped read-only review of diff, tests and provenance; fix material findings and re-run covering tests. Record source/executable hash, no commit.

## Task 2: Identify reused truth and prepare convergence inputs

**Files:** Read existing reports/CSV/HDF5/writer semantics only; create `docs/superpowers/reports/2026-10-14-p4-reuse-audit.md` and new `config/config_srp_modes_DRO-2_half_hour.txt` using the existing complete C2 configuration with only CSV/metadata paths changed to the half-hour immutable dataset.

- [ ] Audit C0 identity, 361 hours, numerical status and nominal states; inspect C1-H final sparse polynomial for an eta±0.3 subdomain without rerunning.
- [ ] Inspect completed `20261005SRP/results_box_wing_initial_only/{earth,moon}/DRO-2` and other existing artifacts. Establish whether every stored nominal state is an independent Real trajectory with matching OPM/config/epoch/common physics/tolerances. If unavailable, name the missing proof and the minimal independent Real truth interface needed.
- [ ] Create the half-hour config by apply_patch; assert all non-path keys exactly match baseline config and run the existing wrapper preflight. Do not re-finalize mode metadata or modify original data.

## Task 3: Complete independent Real endpoints and analysis tool

**Files:** Create `tools/analyze_srp_modes_case.py`, `test/test_analyze_srp_modes_case.py`; only if needed add a focused Real-only executable/test and wrapper with immutable input identity. Do not alter core SRP formulas or integrators.

**Interfaces:** Analysis consumes baseline/half-hour/tighter/higher-DA formal directories and independent Real endpoint truth identified by Task 2. It validates report flags, input/artifact hashes, variable roles/grid/common initial identity, unique contiguous 361-hour records and same physical domain; aligns by hour/point identity, emits classified JSON plus inspectable per-hour CSV in a new exclusive analysis directory.

- [ ] TDD synthetic complete tiny runs: same-state cases yield zero comparisons; known offsets yield known norms; changed mode sampling is allowed only for that variant; OPM/domain/grid/non-SRP mismatches, missing/duplicate hours, nonfinite values, tampered hashes and failed numerical gates are rejected.
- [ ] Compute convergence across identical 512 study points relative to a common C0 nominal trajectory, so numerical differences of reference subtraction do not fake a change in domain size. Emit the <=5% metric checks and actual maximum state differences, while keeping numerical residuals separately classified.
- [ ] Quantify independent Real proxy Earth/Moon versus physical Box-Wing Earth/Moon every hour and full-arc maxima. Compare ADS endpoint evaluations from saved final snapshots against matching proxy Real final endpoints. Record reused truth provenance; do not imply proxy-model error is a software gate.
- [ ] If no valid truth can be reused, TDD and run one Real-only multi-scenario endpoint check with Deta eta=+0.3, proxy Earth (s=1,u=.5), proxy Moon (s=0,u=1), full Earth/Moon. Reuse P1/C0 for D0 instead of a new control domain. Run 360 h and validate final-time completion, finite states, identities and no DA ownership.
- [ ] Restrict the saved final C1-H polynomial to eta±0.3 on the matched state grid only when map identity and numerical validity support it; label this as final-time subdomain reconstruction, not a newly propagated or hourly C1-30 campaign.
- [ ] Run focused synthetic/native tests and task review before claiming analysis results.

## Task 4: Execute three full-arc convergence runs and regressions

- [ ] After launch code is frozen, launch half-hour mode table, rel/abs 1e-13, and DA order 5 C2 runs with distinct output directories under main `SRP/`. Reuse baseline; max depth/output/grid unchanged. Independent jobs may run concurrently if machine capacity allows; no duplicate runs.
- [ ] Observe each actual process/log through completion; hash checks and `formal_result=true` are necessary, not just PID exit. Preserve failed staging and diagnose rather than relax thresholds.
- [ ] Run full Fortran test suite and relevant Python tests with serialized FPM access. Investigate failures with systematic-debugging and regress only scoped corrections; distinguish missing old fixtures/environment from propagation defects.
- [ ] Reconstruct/audit final snapshots and check all 361 hourly records, shape/probe counts and handle baseline. Run Task 3 analysis after all inputs have finished publication.

## Task 5: Final complete-case gate and durable report

- [ ] Independent broad review of current case tooling and numerical methodology; resolve material findings with covering tests.
- [ ] Independently recheck immutable outputs/input hashes, actual last epoch/hour, gates, precision flags and convergence classification. Report any proxy limitation or missing reusable control evidence explicitly.
- [x] Save final complete-case report at main `docs/superpowers/reports/2026-10-15-boxwing-modes-p4-complete-case-results.md` and update design/plan/run ledger. Include exact commands, metrics and source identities, and links to results.
- [ ] Only hand off as a completed full orbit case after this evidence exists; no merge/commit/PR or unrelated campaign modifications.
