#!/usr/bin/env bash
# UT/ADS x C0,C1,C2,C3,T1,T2, 15 days. Does not overwrite an existing run.
# Usage: JOBS=3 BACKEND=auto bash .sh/halo_15d.sh /absolute/new/output
set -euo pipefail
SRP_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
OUT="${1:-$SRP_ROOT/SRP/unified_15d_$(date +%Y%m%d_%H%M%S)}"
JOBS="${JOBS:-3}"
BACKEND="${BACKEND:-auto}"
DRY_RUN="${DRY_RUN:-0}"
[[ "$JOBS" =~ ^[1-9][0-9]*$ ]] || { printf 'JOBS must be a positive integer\n' >&2; exit 2; }
[[ ! -e "$OUT" ]] || { printf 'Output exists; choose a new directory: %s\n' "$OUT" >&2; exit 2; }
[[ "$BACKEND" =~ ^(auto|parallel|python)$ ]] || { printf 'BACKEND must be auto, parallel or python\n' >&2; exit 2; }
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
read -r -a cases <<< "${CASES:-C0 C1 C2 C3 T1 T2}"
read -r -a methods <<< "${METHODS:-UT ADS}"
for case_name in "${cases[@]}"; do
  [[ "$case_name" =~ ^(C0|C1|C2|C3|T1|T2)$ ]] || { printf 'Invalid case: %s\n' "$case_name" >&2; exit 2; }
done
for method in "${methods[@]}"; do
  [[ "$method" =~ ^(UT|ADS|MC|DAMC|DAMC-auto|DAMC-[1-9][0-9]*)$ ]] || { printf 'Invalid method: %s\n' "$method" >&2; exit 2; }
done
common=(--orbits L1Halo-1 --hours 360 --save-hours 1 --count "${COUNT:-50000}"
        --seed "${SEED:-20261005}" --order 4 --max-order 6 --probes 64
        --pos-tol-km 0.1 --vel-tol-kms 1e-6
        --kernel-dir "${POD_KERNEL_DIR:-$SRP_ROOT/kernels}")
if [[ "$DRY_RUN" == 1 ]]; then
  printf '# Build once: bash %q --build-only\n' "$SRP_ROOT/app/run.sh"
  for case_name in "${cases[@]}"; do
    for method in "${methods[@]}"; do
      printf '%q ' bash "$SRP_ROOT/app/run.sh" --executable '<built-executable>' "${common[@]}" \
        --jobs 1 --output "$OUT/${case_name}_${method}" --cases "$case_name" --methods "$method"
      printf '\n'
    done
  done
  exit 0
fi
if [[ "$BACKEND" == auto ]]; then
  if command -v parallel >/dev/null 2>&1; then BACKEND=parallel; else BACKEND=python; fi
fi
if [[ "$BACKEND" == parallel ]] && ! command -v parallel >/dev/null 2>&1; then
  printf 'GNU parallel is not installed. Install it or set BACKEND=python.\n' >&2
  exit 2
fi
# Build in the parent only: child jobs never race on fpm/build files.
exe="$(bash "$SRP_ROOT/app/run.sh" --build-only)"
[[ -x "$exe" ]] || { printf 'Build returned no executable: %s\n' "$exe" >&2; exit 2; }
printf 'Backend=%s jobs=%s output=%s\n' "$BACKEND" "$JOBS" "$OUT"
if [[ "$BACKEND" == python ]]; then
  exec bash "$SRP_ROOT/app/run.sh" --executable "$exe" "${common[@]}" \
    --output "$OUT" --jobs "$JOBS" --cases "${cases[@]}" --methods "${methods[@]}"
fi
mkdir -p -- "$OUT"
# Quote fixed paths; only validated case/method tokens are substituted.
printf -v job_command '%q ' bash "$SRP_ROOT/app/run.sh" --executable "$exe" "${common[@]}"
printf -v quoted_out '%q' "$OUT"
job_command+="--jobs 1 --output ${quoted_out}/{1}_{2} --cases {1} --methods {2}"
parallel --jobs "$JOBS" --joblog "$OUT/parallel_joblog.tsv" --halt soon,fail=1 \
  "$job_command" ::: "${cases[@]}" ::: "${methods[@]}"
printf 'Finished. Reports: %s/*/campaign_summary.csv\n' "$OUT"
