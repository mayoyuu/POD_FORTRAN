#!/usr/bin/env bash
set -euo pipefail
project_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
experiment_root="${SRP_MC_EXPERIMENT_ROOT:-${project_root}/SRP/20261005SRP}"
pid_file="${experiment_root}/logs/campaign.pid"
if [[ ! -f "${pid_file}" ]]; then
    echo "No campaign PID file: ${pid_file}"
    exit 1
fi
pid="$(cat "${pid_file}")"
if [[ ! "${pid}" =~ ^[0-9]+$ ]]; then
    echo "Invalid PID file: ${pid_file}" >&2
    exit 2
fi
if ! kill -0 "${pid}" 2>/dev/null; then
    echo "Campaign PID ${pid} is not running"
    exit 1
fi
command_line="$(tr '\0' ' ' <"/proc/${pid}/cmdline")"
if [[ "${command_line}" != *run_srp_mc_campaign.py* ]]; then
    echo "PID ${pid} does not supervise run_srp_mc_campaign.py" >&2
    exit 2
fi
target_pid="${pid}"
children_file="/proc/${pid}/task/${pid}/children"
if [[ "${command_line}" == *"setsid --wait"* && -r "${children_file}" ]]; then
    child_pids="$(cat "${children_file}")"
    if [[ -n "${child_pids}" ]]; then
        target_pid="${child_pids%% *}"
    fi
fi
kill -TERM -- "-${target_pid}"
echo "Sent SIGTERM to campaign process group ${target_pid} (supervisor ${pid})"
