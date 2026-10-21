# Main-branch integration, 2026-10-21

User-directed update: discontinue the additional pressure test and do not run the
15-day campaign on the user's behalf. Preserve the existing Git repository,
worktrees, old entry points, and all three pre-existing local modifications.

- [x] Record the current main commit with a safety branch. Commit the three
  existing modifications separately so their provenance remains inspectable.
- [x] Merge `codex/boxwing-modes-p1` into main without resetting or deleting files.
- [x] Add only the unified driver, its Python helpers/tests, portable kernel-path
  change, mode inputs and documentation. Preserve the main fpm manifest and let
  its existing executable discovery find the new driver.
- [x] Add `--build-only` and `.sh/halo_15d.sh`: build once, then dispatch the
  UT/ADS x six-case matrix using GNU parallel or the existing Python job pool.
  Use 360 hours, hourly output, 50000 population, 2000 saved non-UT points,
  order 4, 100 m / 1 mm/s checks. No one-order/strict-tolerance pressure test.
- [x] Run focused unit tests, Bash syntax/dry-run checks, a release build and
  one short ADS/T2 UT check in main. Do not launch the 15-day calculation.
- [x] Commit the integration and verify main contains the modes branch, the
  working tree is clean for tracked files, and the original worktrees remain.

No remote push, worktree removal, library deletion, or history rewrite is part
of this task. Historical outputs from the discontinued pressure test are left
untouched but are not acceptance evidence.
