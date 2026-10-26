# Unified SRP propagation (main branch)

## Required compiler: gfortran >= 13.3.0

**GNU Fortran 13.3.0 or newer is required. Older versions are unsupported;
do not use them for long-running or batch DA/ADS/DAMC propagation because
missing finalization can cause severe, sustained native-memory leaks.**

The scalar operator-result/defined-assignment diagnostic reproduced missing
temporary finalizers on gfortran 11.4.0 (100 expressions, 100 unreleased mock
handles). The same source released all temporary handles with gfortran 13.3.0
at both `-O0` and `-O3`. Other versions below 13.3.0 have not all been tested;
13.3.0 is the project's supported minimum, not a claim that every older release
has the same defect. Numerically correct output is not proof of memory safety.

Choose the actual compiler before building or launching batch scripts:

```bash
gfortran-13 --version                   # Must report >= 13.3.0
export FPM_FC="$(command -v gfortran-13)"   # Or use its full installed path
"$FPM_FC" --version
```

Check the version itself: a binary named `gfortran-13` can still be older than
13.3.0. Keep existing dependency/link settings, back up the old build directory,
and rebuild all project Fortran sources, including the DA module, with this
compiler. Do not mix old `.mod`, `.o`, project archives or executables into the
new build. Rebuilding only the C++ `dace_wrapper` does not restore missing
Fortran finalizer calls. Already-built binaries and running processes do not
change when a compiler is installed or an environment variable is exported.
First verify real DA handle cleanup with a short run before resuming a long
campaign; the version requirement alone does not certify all paths leak-free.
See the compiler and memory-leak warning in [the main README](README.md).

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
earth-pointing box-wing; T1 cannonball scale error `eta`; T2 existing bounded
time-varying attitude mode **surrogate**, not a live physical box-wing attitude
error model. L1Halo-1 mode data cover exactly 360 hours and retain original audit
hashes. UT/MC outputs have no independent truncation-error certificate; ADS uses
finite probes and a depth cap, not a global mathematical error guarantee.

T1 `eta` prior (`--eta-prior`, default `normal`): `--eta-prior normal` keeps the
legacy truncated-normal design, `z7 ~ N(0,1)`, `|z7| <= --domain-sigma`, so
`eta = --eta-sigma * z7` and the SRP scale factor is `1+eta`. `--eta-prior
uniform` makes `z7` uniform on `[-domain-sigma, +domain-sigma]`, i.e. `eta`
uniform on `+-domain-sigma*--eta-sigma`; the support must stay inside the
physical bound `eta >= -1`, enforced as `domain-sigma*eta-sigma <= 1`. With the
common `--domain-sigma 3 --eta-sigma 1/3` the support is exactly `[-1,+1]`, so
SRP=0 (5% of the mass within 10% of zero SRP) is represented instead of sitting
on the 3-sigma edge. The prior only changes `z7`; `z1..z6,s,u` and therefore the
C0-C3 populations are bit-identical to the normal design, and `eta_prior` plus
`eta_support` are recorded in `campaign.json` and every `population.json`.
Uniform-prior T1 results are not directly comparable with legacy sigma=0.3 T1
results.

## Verifying a campaign

Two read-only checkers audit an existing campaign directory without propagating
anything; both exit nonzero on failure and can write a JSON report.

Before propagation (works on a `--prepare-only` tree):

```bash
python3 tools/check_t1_eta_prior.py --root /absolute/output --eta-sigma 0.3333333333333333
```

It checks the recorded `eta_prior`/`eta_support`, the `[-1,1]` physical bound and
the ADS domain, a Kolmogorov-Smirnov test of eta against `U(-a,+a)` (default
alpha 1e-3, so 17 orbits still have a small false-alarm rate), the uniform
standard deviation, and that the SRP~0 tail at `eta=-1` is populated at the rate
`(a-0.99)/(2a)`. It also proves C0-C3 keep `eta=0` and share one bit-identical
6D design hash with T1, so the prior switch cannot hide inside the states.

After propagation, with `hourly_moments.csv` present:

```bash
python3 tools/check_t1_c0_coverage.py --root /absolute/output --methods UT ADS
```

Per hour it reports the 6D Mahalanobis shift of the C0 mean (m6), the smallest
eigenvalue of the C0-whitened `Sigma_T1 - Sigma_C0` (0 = equal variance, +1 =
T1 doubles the C0 variance, negative = T1 tighter than C0) and the covariance
volume ratio. Defaults: m6 <= 6, whitened lambda >= -0.05 for Monte-Carlo noise,
volume ratio >= 0.9. A T1 prior that truly reaches SRP~0 must contain C0, so a
negative `lambda` beyond noise or `m6` above the limit is a real regression.

Both tools are covered by `test/test_srp_unified_checks.py`.

## Direct interface and migration

```bash
bash app/run.sh --help
bash app/run.sh --build-only
bash app/run.sh --methods UT ADS --cases C0 C1 C2 C3 T1 T2 \
  --hours 360 --save-hours 1 --count 50000 --jobs 3 --output /absolute/new/output

# T1 with a bounded uniform eta prior on [-1,1] instead of the legacy normal one:
bash app/run.sh --methods UT ADS --cases C0 C1 C2 C3 T1 --eta-prior uniform \
  --eta-sigma 0.3333333333333333 --domain-sigma 3 \
  --hours 360 --save-hours 1 --count 50000 --jobs 5 --output /absolute/new/uniform

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

On a server, install fpm and gfortran >= 13.3.0, and provide compatible SPICE, DACE C++ wrapper
and OpenBLAS libraries. Preserve/rebuild `lib/` for the target ABI; binary
libraries and kernels are not necessarily part of a Git clone. Two data files are
mandatory and are fatal when absent: `kernels/gravity_models/GGM05C.GEO` (Earth,
WGS84/GGM05C, fully normalized, first record `GM R` in SI units) and
`kernels/gravity_models/gggrx_0660pm_sha.tab` (Moon). `kernels/spk/de440.bsp`
alone is not enough for the default `earth_degree = 10`, `moon_degree = 10`
configuration. A `libspicelib.a` built by an older compiler can link but crash
inside `furnsh`; rebuild it from the matching SPICE toolkit sources with the
same gfortran as the application. Install Python dependencies into a local
virtual environment:

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
