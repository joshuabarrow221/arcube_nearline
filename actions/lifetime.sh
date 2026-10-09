#!/usr/bin/env bash
# One watcher action for completed FLOW files and immutable slow-control snapshots.
# init.inc.sh supplies paths plus errexit/pipefail. Scientific fit failures
# are status rows in history; process/IO failures still fail the worker task.
# The configured profile owns all paths. Carry its absolute config in the
# queued command so worker environments need not inherit watcher shell vars.
python_cmd=${ARCUBE_NEARLINE_PYTHON:-python3}
central="$(dirname "${BASH_SOURCE[0]}")/lifetime/lifetime.py"
if [[ "${1:-}" == --monitor-config ]]; then
    exec "$python_cmd" "$central" "$@"
elif [[ "${2:-}" == --monitor-config ]]; then
    input=$1
    shift
    exec "$python_cmd" "$central" "$@" --input-file "$input"
elif [[ -n "${ARCUBE_NEARLINE_LIFETIME_CONFIG:-}" ]]; then
    configured=(--monitor-config "$ARCUBE_NEARLINE_LIFETIME_CONFIG")
    if [[ $# -gt 0 ]]; then configured+=(--input-file "$1"); shift; fi
    exec "$python_cmd" "$central" "${configured[@]}" "$@"
fi
# No config preserves the historical command-line interface for old workflows.
stage=lifetime
source "$(dirname "${BASH_SOURCE[0]}")/../lib/init.inc.sh"
script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/lifetime" && pwd)
python_cmd=${ARCUBE_NEARLINE_PYTHON:-python3}
inbase="$data_root/$nearline_name/flowed_charge"
jsonpath="$json_outbase/json_elifetime.json"
globalplotpath="$plot_outbase/elifetime_time_series.png"
mkdir -p "$json_outbase" "$plot_outbase" "$log_outbase"
# The central command publishes both the mean and _raw plots under this name.
args=(--output_file_json "$jsonpath" --output-timeseries "$globalplotpath")
logpath="$log_outbase/refresh.log"
if [[ $# -gt 0 ]]; then
    inpath=$(realpath "$1")
    shift
    if [[ "$inpath" == *.json ]]; then
        args+=(--slow-controls "$inpath")
    else
        # Read only segment datasets from FLOW; total file size is not a fit limit.
        reldir=$(dirname "${inpath#"$inbase/"}")
        reldir=${reldir#/}
        plotpath="$plot_outbase/$reldir/$(basename "$inpath").lifetime.png"
        logpath="$log_outbase/$reldir/$(basename "$inpath").lifetime.log"
        mkdir -p "$(dirname "$plotpath")" "$(dirname "$logpath")"
        args+=(--input_file "$inpath" --output_file_plot "$plotpath")
        if [[ "${ARCUBE_NEARLINE_SELECT_MUONS:-0}" == 1 ]]; then
            args+=(--select-muons)
        fi
        if [[ "${ARCUBE_NEARLINE_WRITE_LIFETIME_METADATA:-0}" == 1 ]]; then
            args+=(--write-hdf5-metadata)
        fi
        # Packet pairing follows the existing mirrored directory tree. A late
        # packet requires a deliberate retry; watcher.py deduplicates by path.
        if [[ "$inpath" == "$inbase/"* ]]; then
            packet_relative=${inpath#"$inbase/"}
            packet_relative=${packet_relative%.FLOW.hdf5}.h5
            packetpath="$data_root/$nearline_name/packet/$packet_relative"
            if [[ -f "$packetpath" ]]; then
                args+=(--packet-file "$packetpath")
                # Preserve the entire acquisition hierarchy: Step1/trial1
                # and Step2/trial1 must not silently share a packet pool.
                # An operator's verified cohort remains an explicit override.
                if [[ -z "${ARCUBE_NEARLINE_PACKET_COHORT:-}" ]]; then
                    args+=(--packet-cohort "$reldir")
                fi
            else
                echo "Matching native packet file unavailable: $packetpath"
            fi
        fi
    fi
fi
snapshot_dir=${ARCUBE_NEARLINE_PURITY_SNAPSHOT_DIR:-$data_root/$nearline_name/slow_controls}
# Aggregate across snapshot boundaries: a three-hour plateau or a daily mean
# can cross export-file boundaries. Keep snapshots in this flat directory.
shopt -s nullglob
for snapshot in "$snapshot_dir"/*.json; do
    args+=(--slow-controls "$snapshot")
done
if [[ -n "${ARCUBE_NEARLINE_PACKET_EXT_IO:-}" ]]; then
    args+=(--ext-io "$ARCUBE_NEARLINE_PACKET_EXT_IO")
fi
if [[ -n "${ARCUBE_NEARLINE_PACKET_PEDESTAL:-}" ]]; then
    args+=(--pedestal "$ARCUBE_NEARLINE_PACKET_PEDESTAL")
fi
if [[ -n "${ARCUBE_NEARLINE_PACKET_PEDESTAL_SOURCE:-}" ]]; then
    args+=(--packet-pedestal-source "$ARCUBE_NEARLINE_PACKET_PEDESTAL_SOURCE")
fi
if [[ -n "${ARCUBE_NEARLINE_PACKET_COHORT:-}" ]]; then
    args+=(--packet-cohort "$ARCUBE_NEARLINE_PACKET_COHORT")
fi
if [[ -n "${ARCUBE_NEARLINE_PACKET_SAMPLES:-}" ]]; then
    args+=(--packet-samples "$ARCUBE_NEARLINE_PACKET_SAMPLES")
fi
args+=(--packet-window-hours "${ARCUBE_NEARLINE_PACKET_WINDOW_HOURS:-24}")
if [[ -n "${ARCUBE_NEARLINE_GAS_CONVERSION_CONFIG:-}" ]]; then
    args+=(--gas-conversion-config "$ARCUBE_NEARLINE_GAS_CONVERSION_CONFIG")
fi
if [[ -n "${ARCUBE_NEARLINE_GAS_QUALITY_CONFIG:-}" ]]; then
    args+=(--gas-quality-config "$ARCUBE_NEARLINE_GAS_QUALITY_CONFIG")
fi
if [[ -n "${ARCUBE_NEARLINE_TRACK_REVIEW:-}" ]]; then
    args+=(--track-review "$ARCUBE_NEARLINE_TRACK_REVIEW")
fi
"$python_cmd" "$script_dir/lifetime.py" "${args[@]}" "$@" 2>&1 | tee "$logpath"
