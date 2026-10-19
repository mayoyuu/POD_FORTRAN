# Box-Wing modes probability results: portable export and analysis

## Approved scope

Copy, do not move, the three existing probability campaign directories from
`POD_Fortran/SRP` into the new desktop directory
`C:\Users\13113\Desktop\科研\误差传播\SRP影响\20261005SRP\boxwing_modes_probability_20261019`.
Preserve original directory names, all saved hours, snapshots, logs, failed runs,
and interrupted staging directories. Do not overwrite the older desktop experiment.
Do not restart propagation, change tolerances, or promote diagnostic results to
formally accepted results. NRHO remains on hold; only the original 15 orbits are
analyzed. Copy verification compares relative file names, sizes, and SHA-256.

The user approved the analysis scope in chat on 2026-10-19: default window 0–300 h,
50,000-sample hourly statistics, saved 2,000-point clouds, continuous numerically
passing prefixes for distribution plots, and separate diagnostics for failed data.

## Portable interface and files

Deliver `analyze_boxwing_probability.py`, dependency requirements, and a Chinese
README in the export root. The script resolves inputs relative to itself, accepts
`--root`, `--output`, and `--max-hour` (default 300), and runs on Windows Python
without WSL, Fortran, DACE, or SPICE. Outputs go to `analysis/` by default.
Keep the maintained source and tests in the existing linked source worktree.

Discover the 15 known orbits only in the three named campaign directories. Accept
published DRO directories and explicitly identified failed/interrupted staging
directories; reject duplicate orbit candidates instead of selecting silently.
Never use the older desktop paired-experiment data as a substitute. A manifest
records source paths, source revision, copy checks, and input/output conventions.

## Interpretation and analysis

Read `c2_probability_hourly_stats.csv`: select the 50,000-sample rows for primary
means, complete 6×6 covariance, standard deviations and quantiles. Keep 10,000 and
20,000 prefixes as supporting convergence evidence, not extra hourly observations.
Join `c2_history.csv` by hour for original-probe and zero-mode diagnostics. Use
`ads_accuracy_valid`, `precision_confirmed`, finite values, and continuity to
identify the recorded passing prefix. Missing hours or malformed records must be
reported; a trailing incomplete record must never be silently considered complete.
The resulting prefix is a recorded numerical-check prefix, not a new formal
acceptance or a proof of accuracy throughout the complete uncertainty domain.

Record original-probe and probability-probe position/velocity errors (limits
0.1 km and 1e-6 km/s), zero-mode errors (1e-6 km and 1e-9 km/s), depth-limited
patches, first invalid hour, and the directly observed cause(s). If a cumulative
failure flag cannot be explained by available rows, report an unknown recorded
failure rather than inventing a cause. Keep completed accepted, completed failed,
and interrupted/no-final-report states distinct.

Generate per-orbit summary and hourly diagnostic CSV files, the full hourly
means/covariance in a NumPy archive, and static figures for numerical errors with
tolerance lines, mean offsets and standard deviations, and probability quantiles.
Distribution plots use only the continuous passing prefix and state its last hour.
Plot representative saved 2,000-point position and velocity projections at that
last passing hour, with covariance ellipses from the 50,000-sample statistics.
An ellipse is a covariance contour, not a guaranteed probability contour for a
non-Gaussian distribution. Read accepted DRO point files from HDF5 and failed or
interrupted runs from their raw CSVs. Stream large CSV clouds rather than loading
all hourly points into memory. Do not recompute statistics from the 2,000 points.

Clearly distinguish physical dispersion from ADS approximation errors. These
errors compare ADS against Real integration of the same SRP surrogate; they do
not measure error relative to the full live Box-Wing model. The 50,000 samples
are evaluations of ADS maps, not 50,000 independent Real integrations. Preserve
covariance units: position km², velocity (km/s)², cross blocks km²/s.

## Failure handling and verification

Invalid segments remain in exported raw data and hourly diagnostics, but do not
appear as trustworthy distribution results. A missing final report on an
interrupted case is expected and explicitly labeled. Missing required files,
duplicate hours/orbits, malformed rows, nonfinite statistics, and unsupported point
schemas produce a clear diagnostic; no silent zero filling or extrapolation.

Use unit tests with small CSV/HDF5 fixtures for discovery, 50,000-row selection,
continuity, individual failure causes including zero-mode-only failure, missing
reports, malformed/truncated data, complete covariance order, and point schemas.
Run the script against the real copied 15 cases, inspect generated figures, and
verify desktop paths and output counts before claiming delivery. Preserve raw
inputs and distinguish copy verification from numerical validation.
