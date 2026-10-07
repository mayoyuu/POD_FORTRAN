#!/usr/bin/env bash
set -eo pipefail
project_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
python_bin="${SRP_MC_PYTHON:-/home/songyu/venvs/astro/bin/python}"
experiment_root="${SRP_MC_EXPERIMENT_ROOT:-${project_root}/SRP/20261005SRP}"
mkdir -p "${experiment_root}/logs"
stamp="$(date +%Y%m%d_%H%M%S)"
log="${experiment_root}/logs/campaign_${stamp}.log"
pid_file="${experiment_root}/logs/campaign.pid"
if [[ -f "${pid_file}" ]]; then
    old_pid="$(cat "${pid_file}")"
    if [[ "${old_pid}" =~ ^[0-9]+$ ]] && kill -0 "${old_pid}" 2>/dev/null; then
        echo "Campaign already running with PID ${old_pid}" >&2
        exit 2
    fi
fi
source "${project_root}/setup_env.sh" >/dev/null
set -u
nohup setsid --wait "${python_bin}" -u "${project_root}/tools/run_srp_mc_campaign.py" \
    --project-root "${project_root}" \
    --experiment-root "${experiment_root}" "$@" \
    >"${log}" 2>&1 </dev/null &
pid=$!
printf '%s\n' "${pid}" >"${pid_file}"
printf 'PID=%s\nLOG=%s\nEXPERIMENT_ROOT=%s\n' "${pid}" "${log}" "${experiment_root}"
