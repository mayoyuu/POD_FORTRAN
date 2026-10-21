#!/usr/bin/env bash
# Usage: JOBS=3 bash .sh/halo_15d_fpm_parallel.sh /absolute/new/output
# Requires an existing GNU parallel installation. Does not install anything.
set -euo pipefail
ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
cd -- "$ROOT"
export OUT="${1:-$ROOT/SRP/halo_15d_fpm_$(date +%Y%m%d_%H%M%S)}"
JOBS="${JOBS:-3}"
[[ "$JOBS" =~ ^[1-9][0-9]*$ ]] || { printf 'Invalid JOBS\n' >&2; exit 2; }
if [[ "${DRY_RUN:-0}" == 1 ]]; then
  sed -n '/^fpm run .*--cases /p' "${BASH_SOURCE[0]}"
  exit 0
fi
command -v parallel >/dev/null || { printf 'GNU parallel is unavailable; nothing started.\n' >&2; exit 2; }
[[ ! -e "$OUT" ]] || { printf 'Output exists: %s\n' "$OUT" >&2; exit 2; }
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export FPM_FFLAGS=-ffree-line-length-none
export DACE_ROOT="${DACE_ROOT:-$HOME/local_dace}"
export OPENBLAS_ROOT="${OPENBLAS_ROOT:-$HOME/local_openblas}"
export LIBRARY_PATH="$ROOT/lib:$DACE_ROOT/lib:$OPENBLAS_ROOT/lib:${SPICE_LIB_DIR:-$ROOT/lib}:${LIBRARY_PATH:-}"
export LD_LIBRARY_PATH="$ROOT/lib:$DACE_ROOT/lib:$OPENBLAS_ROOT/lib:${SPICE_LIB_DIR:-$ROOT/lib}:${LD_LIBRARY_PATH:-}"
mkdir -p -- "$(dirname -- "$OUT")"
mkdir -- "$OUT"
# Let fpm prepare this target once before concurrent jobs; do not edit sources during a campaign.
command fpm run --target run_srp_unified --profile release -- --help >"$OUT/fpm_prepare.log" 2>&1
# --runner retains the common sample/UT/HZ interface and passes fpm's executable to it.
nohup parallel --jobs "$JOBS" --joblog "$OUT/joblog.tsv" >"$OUT/nohup.log" 2>&1 <<'COMMANDS' &
fpm run --target run_srp_unified --profile release --runner 'bash app/run.sh --executable' -- --methods UT --cases C0 --hours 360 --save-hours 1 --count 50000 --jobs 1 --output "$OUT/C0_UT"
fpm run --target run_srp_unified --profile release --runner 'bash app/run.sh --executable' -- --methods ADS --cases C0 --hours 360 --save-hours 1 --count 50000 --jobs 1 --output "$OUT/C0_ADS"
fpm run --target run_srp_unified --profile release --runner 'bash app/run.sh --executable' -- --methods UT --cases C1 --hours 360 --save-hours 1 --count 50000 --jobs 1 --output "$OUT/C1_UT"
fpm run --target run_srp_unified --profile release --runner 'bash app/run.sh --executable' -- --methods ADS --cases C1 --hours 360 --save-hours 1 --count 50000 --jobs 1 --output "$OUT/C1_ADS"
fpm run --target run_srp_unified --profile release --runner 'bash app/run.sh --executable' -- --methods UT --cases C2 --hours 360 --save-hours 1 --count 50000 --jobs 1 --output "$OUT/C2_UT"
fpm run --target run_srp_unified --profile release --runner 'bash app/run.sh --executable' -- --methods ADS --cases C2 --hours 360 --save-hours 1 --count 50000 --jobs 1 --output "$OUT/C2_ADS"
fpm run --target run_srp_unified --profile release --runner 'bash app/run.sh --executable' -- --methods UT --cases C3 --hours 360 --save-hours 1 --count 50000 --jobs 1 --output "$OUT/C3_UT"
fpm run --target run_srp_unified --profile release --runner 'bash app/run.sh --executable' -- --methods ADS --cases C3 --hours 360 --save-hours 1 --count 50000 --jobs 1 --output "$OUT/C3_ADS"
fpm run --target run_srp_unified --profile release --runner 'bash app/run.sh --executable' -- --methods UT --cases T1 --hours 360 --save-hours 1 --count 50000 --jobs 1 --output "$OUT/T1_UT"
fpm run --target run_srp_unified --profile release --runner 'bash app/run.sh --executable' -- --methods ADS --cases T1 --hours 360 --save-hours 1 --count 50000 --jobs 1 --output "$OUT/T1_ADS"
fpm run --target run_srp_unified --profile release --runner 'bash app/run.sh --executable' -- --methods UT --cases T2 --hours 360 --save-hours 1 --count 50000 --jobs 1 --output "$OUT/T2_UT"
fpm run --target run_srp_unified --profile release --runner 'bash app/run.sh --executable' -- --methods ADS --cases T2 --hours 360 --save-hours 1 --count 50000 --jobs 1 --output "$OUT/T2_ADS"
COMMANDS
pid=$!
printf '%s\n' "$pid" >"$OUT/parallel.pid"
printf 'Background PID=%s\nOutput=%s\nLog=%s/nohup.log\n' "$pid" "$OUT" "$OUT"
