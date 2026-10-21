# Unified SRP propagation (main branch)

## Location and Git integration

The active project is `/home/songyu/fortran_test/POD_Fortran`, branch `main`.
The Windows `SRP_UNIFIED_20261021` folder was a staging/source snapshot, not a
replacement Git repository. Its minimal fpm.toml is **not** installed over the
main manifest. Existing executables, Git history, libraries and worktrees remain.
The modes branch is merged, then the new files are committed to main. There is
no remote push in this workflow.

New common entry: `app/run.sh` -> `tools/run_srp_unified.py` ->
`app/run_srp_unified.f90`. Older SRP runners remain as research/compatibility
entry points; new campaigns need not copy their per-case loops.

## Run 15 days yourself

In WSL/Linux, from the main project:

```bash
cd /home/songyu/fortran_test/POD_Fortran
JOBS=3 bash .sh/halo_15d.sh "/mnt/c/Users/13113/Desktop/科研/误差传播/SRP影响/SRP_unified_results_20261021/halo_15d_$(date +%Y%m%d_%H%M%S)"
```

This builds once and runs UT/ADS x C0,C1,C2,C3,T1,T2: 12 independent processes,
at most JOBS simultaneously, each single-threaded. 15 days = 360 hours; outputs
are at 0,1,...,360 h. Defaults: population 50000, seed 20261005, order 4,
64 fixed Real probes, tolerances 0.1 km / 1e-6 km/s. These are normal settings,
not the discontinued first-order/nanometre pressure test. Jobs with invalid
precision are recorded as `numerically_invalid` and return nonzero; a completed
process alone is not a 15-day accuracy certificate.

The default backend uses GNU parallel if present, otherwise the existing
Python process pool. A dedicated GNU-only file is provided for a server that
already has parallel; **nothing is installed by these scripts**:

```bash
JOBS=3 bash .sh/halo_15d_parallel.sh /absolute/new/output
```

GNU backend uses `parallel --jobs "$JOBS" --joblog ... --halt soon,fail=1`:
failed jobs stop new dispatch, already running jobs finish. Each task calls
the common app/run.sh with `--jobs 1 --executable ...`; no concurrent builds.
See the [GNU parallel tutorial](https://www.gnu.org/software/parallel/parallel_tutorial.html).

Without installing parallel, force the built-in alternative:

```bash
JOBS=3 BACKEND=python bash .sh/halo_15d.sh /absolute/new/output
```

Or directly run the common entry with its built-in process concurrency:

```bash
bash app/run.sh --methods UT ADS --cases C0 C1 C2 C3 T1 T2 \
  --hours 360 --save-hours 1 --count 50000 --seed 20261005 --jobs 3 \
  --output /absolute/new/output
```

For a command-only inspection (no build, no propagation, no output creation):

```bash
DRY_RUN=1 bash .sh/halo_15d.sh /absolute/new/output
```

Use a **new** output directory each time. Existing results are never overwritten.
JOBS can be increased within available CPU/RAM, but this script does not guess
server capacity. Do not additionally enable nested OpenMP/BLAS parallelism.
No 15-day run is started automatically by the integration or unit tests.

## Results and numerical meaning

- GNU backend: `output/C0_ADS/L1Halo-1/C0/ADS/`, etc.; top-level
  `parallel_joblog.tsv`, each case-method has `campaign_summary.csv` and
  `campaign.json`. Python backend: `output/L1Halo-1/C0/ADS/`, etc., and a
  combined top-level campaign summary.
- `hourly_moments.csv`: all-population mean and full 6x6 covariance; Real-probe
  position/velocity errors, map precision flags, order, patch count and timers.
- `points.csv`, `results.h5`: fixed 2000 sample IDs at every hour for ADS/MC/DAMC
  (or all if N<2000). UT saves its native 13/15/17 weighted sigma points instead.
- `splits.csv`: ADS split interval, parent depth, split direction/component,
  six errors/tolerances/excesses, direction scores, parent centre/width. A
  header-only file means no split events, not a missing diagnostic.
- `hz.csv`: HZ test on the fixed saved random sample, normally N=2000, **not**
  the full population. UT is explicitly not applicable. HZ at h=0 is retained;
  initial truncation can already cause nonnormality. Singular covariance is
  reported, not assigned a fabricated p-value.
- `timing.csv`, `report.json`: propagation, map evaluation, full-population
  moments, validation, I/O; Python HZ and postprocessing timings are separate.
  `compute_excluding_io_seconds` is the measured Fortran computation sum, not
  total process elapsed minus a guessed output cost. Bootstrap/internal kernel
  I/O and parsing overhead are not assigned to pure propagation time. Final
  timing-row write latency cannot be included in that same row.
- Shared initial design and hashes, resolved config, command replay and logs
  are saved. UT matches the shared population's full joint first/second moments
  and propagates the same sigma trajectories without hourly re-Gaussianization.

Case definitions: C0 no SRP; C1 cannonball; C2 moon-pointing box-wing; C3
earth-pointing box-wing; T1 cannonball eta error (sigma=0.3); T2 existing bounded
time-varying attitude mode **surrogate**, not a live physical box-wing attitude
error model. L1Halo-1 mode data cover exactly 360 hours and retain original audit
hashes. UT/MC outputs have no independent truncation-error certificate; ADS uses
finite probes and a depth cap, not a global mathematical error guarantee.

## Direct interface and migration

```bash
bash app/run.sh --help
bash app/run.sh --build-only
bash app/run.sh --methods UT ADS --cases C0 C1 C2 C3 T1 T2 \
  --hours 360 --save-hours 1 --count 50000 --jobs 3 --output /absolute/new/output

# Single custom physical model; explicit case label required for provenance:
bash app/run.sh --methods UT ADS --cases C3 --srp-model box_wing \
  --attitude earth --srp-error none --hours 360 --output /absolute/new/custom

# Optional additional methods, not launched by the default 15-day script:
METHODS='MC DAMC-4 DAMC-auto' JOBS=3 bash .sh/halo_15d.sh /absolute/new/other_methods
```

`--state-scale` scales the OPM uncertainty standard deviation; full 6x6
correlations are retained. `--samples input.csv --count N` accepts standardized
`z1..z6` (optional z7,s,u,sample_id), or physical `x,y,z,vx,vy,vz` (optional
eta,s,u,sample_id), in km / km/s. The row count must match N; input must lie in
the declared uncertainty domain. Missing SRP parameter columns use the shared
generated parameter design. Explicit physical state samples are not rescaled.
Use `--opm` / `--opm-root`, `--mode-root`, `--kernel-dir` as needed. The shipped
mode table is specific to L1Halo-1; do not relabel it for another orbit. Its CSV,
audit metadata and the required L1Halo-1 OPM are explicitly versioned in main.

On a server, install fpm/gfortran and provide compatible SPICE, DACE C++ wrapper
and OpenBLAS libraries. Preserve/rebuild `lib/` for the target ABI; binary
libraries and kernels are not necessarily part of a Git clone. Install Python
dependencies into a local virtual environment:

```bash
python3 -m venv .venv-srp
.venv-srp/bin/python -m pip install -r tools/requirements-srp-unified.txt
export DACE_ROOT=/your/dace/prefix
export OPENBLAS_ROOT=/your/openblas/prefix
export SPICE_LIB_DIR=/your/spice/lib
export POD_KERNEL_DIR=/your/kernels
JOBS=3 BACKEND=python bash .sh/halo_15d.sh /your/new/output
```

Explicit `PYTHON=/path/to/python` overrides discovery. app/run.sh otherwise
tries local .venv-srp/.venv/.venv-probability and the preserved modes worktree's
existing environment, then python3. Existing ~/local_dace and ~/local_openblas
are detected only if their lib directories exist; explicit environment variables
take precedence. Kernel path selection also covers gravity model files.

Focused verification only: `python -m unittest discover -s test -p
'test_srp_unified_*.py' -v`, release build, shell syntax/dry-run and short
normal-tolerance main-branch smoke. Earlier 2-hour UT/ADS six-case evidence is
under the Windows SRP_unified_results_20261021/smoke_2h directory. No extra
pressure tests or 24-hour/15-day certification is claimed.
