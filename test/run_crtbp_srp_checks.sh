#!/usr/bin/env bash
# From any directory: bash /path/to/POD_Fortran/test/run_crtbp_srp_checks.sh
set -eo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
source ./setup_env.sh
set -u
compiler="${FPM_FC:-gfortran-13}"
out="$PWD/build/crtbp_srp_validation"
mkdir -p "$out"
args=(--compiler "$compiler" --profile debug --flag -ffree-line-length-none)
fpm test "${args[@]}" --target test_crtbp_srp | tee "$out/force.log"
fpm test "${args[@]}" --target test_crtbp_srp_propagation -- "$out/15day.csv" | tee "$out/propagation.log"
for scenario in invalid_units sun_collision invalid_parameters da_sun_collision; do
    if fpm test "${args[@]}" --target test_crtbp_srp -- "$scenario" > "$out/$scenario.log" 2>&1; then
        echo "FAIL: $scenario unexpectedly succeeded" >&2
        exit 1
    fi
    case "$scenario" in
        invalid_units) expected='DU, TU, AU, radius and period must be finite and positive' ;;
        sun_collision|da_sun_collision) expected='Sun-spacecraft distance is singular or non-finite' ;;
        invalid_parameters) expected='initialize DACE with exactly the number of active SRP parameters' ;;
    esac
    grep -F "$expected" "$out/$scenario.log" >/dev/null
    echo "PASS: $scenario rejected with expected diagnostic"
done
echo "CRTBP SRP validation outputs: $out"
