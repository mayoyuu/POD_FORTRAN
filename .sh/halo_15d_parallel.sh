#!/usr/bin/env bash
# Dedicated GNU parallel entry; shares every scientific option with halo_15d.sh.
set -euo pipefail
export BACKEND=parallel
exec bash "$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)/halo_15d.sh" "$@"
