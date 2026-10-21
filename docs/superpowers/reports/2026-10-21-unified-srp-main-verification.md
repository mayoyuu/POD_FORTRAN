# Main-branch integration verification (2026-10-21)

Repository: `/home/songyu/fortran_test/POD_Fortran`.

## Preserved Git state

- Before integration: main `5a3b67e`, three local modifications. A safety branch
  `backup/main-before-srp-unified-20261021` points to that main commit.
- Local modifications preserved separately in `8192085` (zero-SRP eta=-1
  boundary, its regression test, and the analysis-directory ignore rule).
- `codex/boxwing-modes-p1` merged into main with merge commit `5aeb190`.
- No repository/worktree/library deletion, reset, history rewrite or remote push.
- Original fpm.toml retained; automatic executable discovery includes the new
  driver. Subsequent builds target only the new executable with `fpm run ...
  -- --help`, not all historical apps. Build logs are under ignored `build/`.

## Focused verification

1. Python: `python -m unittest discover -s test -p 'test_srp_unified_*.py' -v`:
   15 tests, passed. Includes joint moment recovery, common populations, HZ,
   build-only parsing, batch dry-run matrix, GNU-only wrapper and no-overwrite guard.
2. Bash syntax checks: `bash -n` individually on app/run.sh, .sh/halo_15d.sh
   and .sh/halo_15d_parallel.sh, passed. Dry-run
   produces 12 commands, each 360 hours, N=50000, save interval 1 hour and a
   single child job, without creating output or starting propagation.
3. Main release compilation completed successfully; `app/run.sh --build-only`
   returned an executable. Single-target rebuild/help also succeeded.
4. Main normal-tolerance smoke: C1/T2 x UT/ADS, 0..0.05 hours, N=50000,
   two workers, normal 0.1 km / 1e-6 km/s tolerances. All four reports completed
   with `all_hours_numerically_valid=True`. ADS 2000 saved members at both
   epochs; UT native weighted sigma points. This smoke took seconds, not a
   multi-hour or 15-day run.

Evidence is outside Git in the Windows workspace:
`SRP_unified_results_20261021/main_integration_smoke/campaign_summary.csv`.
Earlier 2-hour UT/ADS six-case evidence is in sibling `smoke_2h/`.

GNU parallel is not installed on the current WSL host. Its batch command was
syntax/dry-run checked and its options checked against the official tutorial;
actual GNU scheduler execution is not claimed. The built-in worker backend is
the same one used for the four real main-branch smoke jobs.

The extra first-order/strict-tolerance pressure test is discontinued and is not
acceptance evidence. No new stress test, 24-hour campaign or 15-day propagation
was launched. No independent review sign-off is claimed. Fifteen-day numerical
accuracy and adaptive DAMC escalation are not certified by this short check.
