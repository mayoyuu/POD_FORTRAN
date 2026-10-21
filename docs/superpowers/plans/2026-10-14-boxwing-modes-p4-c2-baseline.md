# Box-Wing Modes P4 C2 Baseline Execution Plan

> **For agentic workers:** Use superpowers:executing-plans inline, with using-git-worktrees and verification-before-completion. This is an experiment execution plan, not a new feature implementation. The user explicitly requests direct C2 execution and reuse of existing controls; do not launch C0/C1 or change core code, commit, merge, or push.

**Goal:** Run the approved 360-hour DRO-2_init C2 baseline through the audited P3 pipeline, while preserving and inventorying existing control results without rerunning them.

**Architecture:** Reuse the isolated `codex/boxwing-modes-p1` worktree and P3's frozen executable. The existing Python wrapper validates immutable inputs, freezes the runner into exclusive staging, executes C2, independently audits outputs/snapshots, and atomically publishes only a valid result. A detached Linux process and exclusive outer log preserve the long run after the launching shell exits.

**Tech Stack:** Existing Fortran/DACE/SPICE executable, Python standard-library wrapper, WSL Linux process management. No new package installation or physics implementation.

## Global Constraints

- Worktree: `/home/songyu/fortran_test/POD_Fortran/.worktrees/boxwing-modes-p1`; preserve all uncommitted P1/P2/P3 changes and runtime symlinks.
- OPM: `OPM/DRO-2/DRO-2_init.opm.json`; config: `config/config_srp_modes_DRO-2.txt`.
- Use immutable finalized 1-hour mode data, 361 nodes over 360 hours; do not modify P1 tables or metadata.
- C2: 8 active variables, s/u in [0,1], DA order 4, max depth 8, domain sigma 3, RKF78 rel/abs tol 1e-12, save every hour, 512 study points and 337 independent Real probes.
- ADS numerical gates: position <=0.1 km, velocity <=1e-6 km/s. Zero-mode gates: position <=1e-6 km, velocity <=1e-9 km/s.
- Output: `/home/songyu/fortran_test/POD_Fortran/SRP/boxwing_modes_DRO-2_init_c2_15d_20261014`; never overwrite an existing directory or launch a duplicate.
- Frozen executable: `/home/songyu/fortran_test/POD_Fortran/SRP/boxwing_modes_DRO-2_init_p3_smoke_20261014_verified/runner`, SHA-256 `89c03449d2d5d72da311cf2a989a721488cf5181299ba50b7ec6770c89a0ef62`.
- Outer log: `/home/songyu/fortran_test/POD_Fortran/SRP/logs/boxwing_modes_DRO-2_init_c2_15d_20261014.session.launch.log`; native hourly progress is in the wrapper-owned staging `run.log`.
- Do not rerun C0/C1-30, alter historical outputs, relabel a 6D result as 7D C1-30, or treat historical ±0.9 eta as ±0.3. Missing comparison identity is a reporting limitation, not permission to rerun.
- Do not stop or reconfigure unrelated existing campaigns. No full-suite builds or core edits while C2 runs.
- This run alone does not complete all P4 proxy-truth comparisons or convergence experiments.

## Task 1: Verify reused baseline and inputs

- [x] Confirm linked worktree/branch and preserve the existing dirty state.
- [x] Compare every current wrapper source hash against the successful P3 report; confirm frozen runner hash and formal P3 status.
- [x] Run wrapper preflight: require finalized P1 gate, matching OPM/common physics/non-SRP config, 361 nodes and 360-hour coverage.
- [x] Confirm output does not exist, sufficient disk space, and no existing C2 run. Existing unrelated Box-Wing jobs remain untouched.

Verification command (from the worktree):

```bash
python3 -c 'import sys,json; from pathlib import Path; sys.path.insert(0,"tools"); import run_srp_modes_ads as m; root=Path("/home/songyu/fortran_test/POD_Fortran/SRP/boxwing_modes_DRO-2_init_p3_smoke_20261014_verified"); p=json.loads((root/"c2_report.json").read_text()); assert p["formal_result"]; assert m.source_hashes()==p["provenance"]["source_files_sha256"]; assert m.file_sha256(root/"runner")==p["provenance"]["executable"]["sha256"]; m.preflight(Path("OPM/DRO-2/DRO-2_init.opm.json"),Path("config/config_srp_modes_DRO-2.txt"))'
```

Expected: exit 0. Preflight completed 2026-10-14; source/executable identity unchanged. Free disk space exceeds 900 GiB.

## Task 2: Start only the C2 baseline

- [x] Launch the command below detached, with unbuffered Fortran output and an exclusively created outer log; record the actual wrapper/child PID and staging path in a separate run note.
- [x] Verify the native command line has `--hours 360`, 8D MODES initialization and actual progress. Do not infer startup success from a PID alone.

Underlying run command:

```bash
env GFORTRAN_UNBUFFERED_ALL=y python3 tools/run_srp_modes_ads.py \
  --opm OPM/DRO-2/DRO-2_init.opm.json \
  --config config/config_srp_modes_DRO-2.txt \
  --output-dir /home/songyu/fortran_test/POD_Fortran/SRP/boxwing_modes_DRO-2_init_c2_15d_20261014 \
  --executable /home/songyu/fortran_test/POD_Fortran/SRP/boxwing_modes_DRO-2_init_p3_smoke_20261014_verified/runner \
  --hours 360 --save-hours 1 --da-order 4 --max-depth 8 \
  --pos-tol-km 0.1 --vel-tol-kms 1e-6
```

Expected: one new staging directory, copied runner, input points and logs. The final output directory remains absent until successful audit/publication.

## Task 3: Observe and verify the result

- [x] Inspect only this process tree, outer log, staging `run.log` and checkpoint CSVs. Record hour/patch/error progress; preserve diagnostics on failure and do not automatically rerun or relax tolerances. Initial cross-call verification reached 23 h with one patch and passing checkpoint errors; final verification remains pending.
- [ ] Once the wrapper ends, require final directory and `formal_result=true`; independently recompute all artifact/input hashes, verify 361 hourly rows, final 360 h, DA handles back to baseline, zero/ADS gates and precision flags, and reconstructed snapshot checks.
- [ ] If the result cannot publish, report the specific numerical/structural failure with its preserved staging path; do not claim a 15-day pass.

## Task 4: Read-only reuse inventory and reporting

- [ ] Inspect existing DRO-2_init reports/manifests for C0 and C1 candidates. Compare actual variable count, eta domain, OPM epoch/state/covariance, force settings, tolerances, duration and identity evidence before claiming matched comparison.
- [ ] Record reuse eligibility or missing evidence without changing historical files or running replacement controls.
- [ ] Update run note and design status to distinguish baseline running/completed/failed from overall P4 completion. Full D1/D2 proxy comparisons and three convergence checks are separate remaining P4 tasks.

No commit, merge, push, PR, control rerun, or unrelated campaign mutation is part of this plan.

## Running checkpoint — 2026-10-14

Wrapper PID 68442, native PID 68466, independent SID 68442; staging suffix `sa3ewcec`. The initial nohup-only process did not persist and produced no staging/output. The existing `nohup setsid --wait` pattern was tested and then launched with verified cross-call survival and actual 8D progress. Empty first-attempt log remains preserved. See `/home/songyu/fortran_test/POD_Fortran/docs/superpowers/reports/2026-10-14-boxwing-modes-p4-c2-run.md`.

Read-only reuse inventory found a matching-OPM 360-hour 6D C0 with matching preserved config and non-SRP settings, plus historical scalar eta±0.9 (C1-H). C1-30 identity remains unconfirmed; no replacement control is authorized or launched. C2 final publication, comparison eligibility and the rest of P4 remain pending.
