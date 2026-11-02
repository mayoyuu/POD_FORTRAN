#!/usr/bin/env bash
# Run from any directory. Logs live in build/, not in the orbit input catalogue.
set -eo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
source ./setup_env.sh
set -u
compiler="${FPM_FC:-gfortran-13}"
out="$PWD/build/crtbp_flatplate_validation"
mkdir -p "$out"
args=(--compiler "$compiler" --profile debug --flag -ffree-line-length-none)
fpm test "${args[@]}" --target test_crtbp_flatplate | tee "$out/force.log"
fpm test "${args[@]}" --target test_crtbp_flatplate_propagation | tee "$out/propagation.log"
for scenario in bad_domain bad_tilt polar_basis wrong_dimension; do
    if fpm test "${args[@]}" --target test_crtbp_flatplate -- "$scenario" > "$out/$scenario.log" 2>&1; then
        echo "FAIL: $scenario unexpectedly succeeded" >&2
        exit 1
    fi
    case "$scenario" in
        bad_domain) expected='chi bounds must lie in [0,1]' ;;
        bad_tilt) expected='maximum tilt must be in [0,90) degrees' ;;
        polar_basis) expected='tangent axis parallel to Sun direction' ;;
        wrong_dimension) expected='DACE variable count must match active parameters' ;;
    esac
    grep -F "$expected" "$out/$scenario.log" >/dev/null
    echo "PASS: $scenario rejected with expected diagnostic"
done
echo "CRTBP flat-plate validation outputs: $out"
