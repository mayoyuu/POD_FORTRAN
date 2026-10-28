#!/usr/bin/env bash
set -euo pipefail
FALCON_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$FALCON_ROOT"
if [[ -z "${PYTHON:-}" ]]; then
  PYTHON=python3
  for candidate in "$FALCON_ROOT/.venv-srp/bin/python" "$FALCON_ROOT/.venv/bin/python" \
    "$FALCON_ROOT/.venv-probability/bin/python" \
    "$FALCON_ROOT/.worktrees/boxwing-modes-p1/.venv-probability/bin/python"; do
    if [[ -x "$candidate" ]]; then PYTHON="$candidate"; break; fi
  done
fi
if [[ -z "${DACE_ROOT:-}" && -d "$HOME/local_dace/lib" ]]; then export DACE_ROOT="$HOME/local_dace"; fi
if [[ -z "${OPENBLAS_ROOT:-}" && -d "$HOME/local_openblas/lib" ]]; then export OPENBLAS_ROOT="$HOME/local_openblas"; fi
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
exec "$PYTHON" "$FALCON_ROOT/.sh/falcon9_impact.py" "$@"
