#!/usr/bin/env bash
set -Eeuo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

tmp_root="$(mktemp -d)"
trap 'rm -rf "$tmp_root"' EXIT

scripts=(
    scripts/run_halo_srp_sweep_parallel.sh
    scripts/run_halo_nominal_srp_sweep_parallel.sh
)

for script in "${scripts[@]}"; do
    bash -n "$script"
    grep -q -- '--zero-init-cov' "$script"
done

ZERO_INIT_COV=1 DRY_RUN=1 DT_LIST=60 \
OUT_ROOT="$tmp_root/srp" bash scripts/run_halo_srp_sweep_parallel.sh

grep -q '^zero_init_cov=1$' "$tmp_root/srp/batch.log"
grep -q '_zeroInitCov' "$tmp_root/srp/queue.tsv"
grep -q 'zero_init_cov' "$tmp_root/srp/summary.tsv"

ZERO_INIT_COV=TRUE DRY_RUN=1 DT_LIST=60 \
OUT_ROOT="$tmp_root/nominal" bash scripts/run_halo_nominal_srp_sweep_parallel.sh

grep -q '^zero_init_cov=1$' "$tmp_root/nominal/batch.log"
grep -q '_zeroInitCov' "$tmp_root/nominal/queue.tsv"
grep -q 'zero_init_cov' "$tmp_root/nominal/summary.tsv"

DRY_RUN=1 DT_LIST=60 \
OUT_ROOT="$tmp_root/default" bash scripts/run_halo_nominal_srp_sweep_parallel.sh

if grep -q '_zeroInitCov' "$tmp_root/default/queue.tsv"; then
    echo 'FAIL: default mode unexpectedly used zero-covariance output suffix'
    exit 1
fi

echo 'test_zero_init_cov_scripts passed'
