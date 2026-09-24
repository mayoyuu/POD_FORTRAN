# SRP–ADS uncertainty-domain history and paired control

## One-orbit runner

Build once, then run from the project root:

```sh
source setup_env.sh
fpm build
fpm run --target run_srp_ads_history -- \
  -opm OPM/L1Halo-1/L1Halo-1_init.opm.json \
  -o SRP/example/L1Halo-1_init \
  --eta-sigma 0.3 \
  --days 15 --save-hours 1 \
  --da-order 4 --max-depth 8 \
  --pos-tol-km 0.1 --vel-tol-kms 1e-6 \
  -cfg config/config.txt
```

The input is the repository's JSON `.opm` or `.opm.json` format with a
complete Cartesian 6×6 covariance. State units are km and km/s, using the
project's GCRF/J2000 convention and configured cannonball SRP model.

With `--eta-sigma 0.3`, the ADS domain is seven dimensional: the whitened
initial-state ±3σ box and zero-mean relative SRP multiplier η in
`[-0.9,+0.9]`. With `--eta-sigma 0`, the same nominal cannonball SRP remains
active, physical η is fixed at zero, and ADS has six variables. The saved
`u7` remains only a stable pairing label, so the two groups retain the same
512 point IDs.

Defaults are 15 days, hourly checkpoints, DA order 4, maximum split depth 8,
0.1 km position validation tolerance, and `1e-6` km/s velocity validation
tolerance. Optional duration flags are `--days` and `--hours`.

## Per-orbit outputs

For output prefix `PREFIX`:

- `PREFIX_history.csv`: 361 rows for the default run, including the nominal
  state, sampled component ranges, sampled position/velocity maxima and
  principal axes, patch/depth counts, seven normalized split counts, queue
  peak, validation discrepancies, and step runtime.
- `PREFIX_shape.csv`: the full 512-point error cloud at every checkpoint,
  keyed by `(hour,point_id)`.
- `PREFIX_initial_points.csv`: stable `point_id`, `u1...u7`, physical
  initial state, physical η, and validation membership. In an η=0 run,
  `u7` is a pairing tag while the `eta` column is zero.
- `PREFIX_split_events.csv`: one row for every actual split. It stores the
  saved propagation interval, parent depth, selected state component and
  uncertainty direction, six truncation estimates/tolerances/positive
  excesses, seven direction scores, and seven-dimensional padded parent
  center/width. A 6D control always has zero events in direction 7.
- `PREFIX_report.json`: format version, ADS dimension, nominal/uncertain SRP
  flags, first arbitrary split, first η-direction split, first maximum-depth
  warning, first validation failure and previous passing hour, plus the full
  decision diagnostics for the first events. Batch execution adds OPM,
  configuration, executable hashes and the source revision.
- `PREFIX_patch_sXXXXXXX_{meta,coeff,splits,map}.csv`: sparse patch
  polynomials and split histories at the start, end, first split, first η
  split, first validation failure, or first unconfirmed checkpoint and its
  predecessor. Seven columns are retained in both groups; the seventh
  coordinate and powers are zero for a 6D control.

A split is taken when at least one DACE truncation estimate exceeds its state
component tolerance. The component is selected by the largest positive
estimate-minus-tolerance and the split direction by the largest top-order
direction score. Its time is reported as the saved interval
`[previous_hour,current_hour]`.

The 207 independent real-integrator validation points are 128 7D corners, 14
axis points, the center, and 64 fixed Halton points. A failure means an ADS
difference above the configured validation threshold at a saved time and
point. It does not measure the physical SRP effect. Reaching maximum depth is
reported separately and does not by itself mean validation failed.

All ranges, principal axes, RMS values, and 95% quantiles described here are
fixed-design sampling metrics. They are not probability confidence intervals
or continuous-domain guarantees.

## Paired 15-day batch

The batch scanner accepts only
`OPM/<one-orbit-directory>/*_init.opm[.json]`. It therefore ignores nested
supplementary batch directories. It also rejects OPMs whose covariance RSS is
not 10 km and 0.03 m/s.

Audit current treatment results for the new schema, hashes, and split log:

```sh
bash tools/run_srp_ads_batch.sh \
  --output-root SRP/SRP_UCERTAINTY_PROP_WITHINIT \
  --eta-sigma 0.3 --audit-only
```

Refresh only missing, legacy, or incompatible treatment results. Each orbit is
written to a unique temporary prefix and replaces the old result only after
the new report passes the complete compatibility audit:

```sh
bash tools/run_srp_ads_batch.sh \
  --output-root SRP/SRP_UCERTAINTY_PROP_WITHINIT \
  --eta-sigma 0.3 --jobs 4 --refresh-incompatible
```

Run the η=0 control campaign:

```sh
bash tools/run_srp_ads_batch.sh \
  --output-root SRP/SRP_UCERTAINTY_PROP_WITHINIT_ETA0 \
  --eta-sigma 0 --jobs 4 --refresh-incompatible
```

The batch stages an immutable copy of the executable so a later build cannot
change binaries midway through a running campaign. Its root files are
`batch_preflight.csv`, `batch_compatibility_audit.csv`,
`batch_results.csv`, `batch_summary.json`, and `batch_build.log`.

To run a campaign in the background:

```sh
mkdir -p SRP/logs
nohup bash tools/run_srp_ads_batch.sh \
  --output-root SRP/SRP_UCERTAINTY_PROP_WITHINIT_ETA0 \
  --eta-sigma 0 --jobs 4 --refresh-incompatible \
  > SRP/logs/eta0_batch.log 2>&1 &
echo $! > SRP/logs/eta0_batch.pid
```

Monitor it with:

```sh
tail -f SRP/logs/eta0_batch.log
ps -fp "$(cat SRP/logs/eta0_batch.pid)"
```

## Paired effect comparison

After both roots are complete:

```sh
python3 tools/compare_srp_ads_pairs.py \
  --treatment-root SRP/SRP_UCERTAINTY_PROP_WITHINIT \
  --control-root SRP/SRP_UCERTAINTY_PROP_WITHINIT_ETA0 \
  --output-root SRP/SRP_UCERTAINTY_EFFECT_WITHINIT
```

The comparator rejects a pair if dimensions, epochs, point IDs, propagation
parameters, OPM/config hashes, or checkpoint count disagree. When both
executable hashes or source revisions are present, those must also agree. It
streams the two shape files by `(hour,point_id)` and does not save a third
pointwise difference cloud.

Each orbit produces:

- `PREFIX_paired_history.csv`: treatment/control bounds and widths, absolute
  increments and amplification factors, axes and volume proxies, both
  centroids and their six-component difference, six-component SRP increment
  extrema, maximum/RMS/design-95% Δr and Δv with worst IDs, 64 u7 corner-pair
  odd/even terms, and both groups' ADS diagnostics and runtime.
- `PREFIX_paired_report.json`: integrity results, sources and hashes, group
  event times and time differences, complete first η split criterion, nominal
  center agreement, and endpoint/peak/peak-hour/time-average/trapezoid
  summaries for the physical-effect metrics.

The comparison root also contains `paired_orbit_summary.csv`,
`paired_population_history.csv`, and `paired_summary.json`.


## Offline reconstruction of a saved ADS domain

A complete patch snapshot consists of four files with a shared prefix:

- `PREFIX_meta.csv`: patch centers, widths, and split depth.
- `PREFIX_coeff.csv`: sparse DA coefficients for six absolute state
  components.
- `PREFIX_splits.csv`: signed split history for every patch.
- `PREFIX_map.csv`: DA order, whitening map, domain sigma, and eta sigma.

Reconstruct a dense deterministic cloud from one saved snapshot by passing the
prefix before those suffixes:

```sh
python3 tools/reconstruct_srp_ads_snapshot.py \
  --prefix SRP/SRP_UCERTAINTY_PROP_WITHINIT/L1Halo-1/run_patch_s1296000
```

The default writes `PREFIX_reconstructed.csv` and
`PREFIX_reconstruction.json`, with 4096 samples per patch. Every patch
includes all local hypercube corners, both signed axes, its center, deterministic
Halton face points, and deterministic Halton interior points. Increase the
density without rerunning orbit propagation:

```sh
python3 tools/reconstruct_srp_ads_snapshot.py \
  --prefix /path/to/run_patch_s1296000 \
  --points-per-patch 16384 \
  --boundary-fraction 0.75 \
  --output /path/to/hour360_dense.csv \
  --report /path/to/hour360_dense.json
```

The tool automatically handles 7D treatment snapshots and padded 6D eta=0
control snapshots. It replays every signed split and rejects metadata that
cannot be reproduced by the split history. When the snapshot follows the
normal `RUN_patch_sNNNNNNN` naming convention, it reads
`RUN_history.csv`, evaluates the original saved design points from
`RUN_initial_points.csv`, and requires agreement with `RUN_shape.csv`
before reporting success. Use `--history` only when the history file was
moved. For a standalone snapshot, `--absolute-only --no-verify` writes
absolute states without the parent run files.

The reconstructed CSV contains patch-local coordinates, original global
whitened coordinates, physical eta, six absolute state components, and, when
history is available, six errors relative to the nominal center trajectory.
The JSON records source hashes, dimensions, sampling parameters, extrema,
structural checks, and saved-cloud verification errors.

This reconstructs the saved ADS polynomial domain. It can make the point cloud
arbitrarily denser, while its physical accuracy remains limited by the DA
order, ADS splitting, maximum depth, and the independent real-integrator
validation recorded by the propagation run. A time can be reconstructed only
when its complete four-file patch snapshot was exported.
