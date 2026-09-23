#!/usr/bin/env bash
set -eo pipefail
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$PROJECT_ROOT"
source ./setup_env.sh
set -u
exec python3 tools/run_srp_ads_batch.py "$@"
