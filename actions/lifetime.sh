#!/usr/bin/env bash
# One watcher action for completed FLOW files and immutable slow-control snapshots.
stage=lifetime
source "$(dirname "${BASH_SOURCE[0]}")/../lib/init.inc.sh"
script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/lifetime" && pwd)
python_cmd=${ARCUBE_NEARLINE_PYTHON:-python3}
inbase="$data_root/$nearline_name/flowed_charge"
jsonpath="$json_outbase/json_elifetime.json"
globalplotpath="$plot_outbase/elifetime_time_series.png"
mkdir -p "$json_outbase" "$plot_outbase" "$log_outbase"
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
        if [[ "$inpath" == "$inbase/"* ]]; then
            packet_relative=${inpath#"$inbase/"}
            packet_relative=${packet_relative%.FLOW.hdf5}.h5
            packetpath="$data_root/$nearline_name/packet/$packet_relative"
            if [[ -f "$packetpath" ]]; then
                args+=(--packet-file "$packetpath")
            else
                echo "Matching native packet file unavailable: $packetpath"
            fi
        fi
    fi
fi
snapshot_dir=${ARCUBE_NEARLINE_PURITY_SNAPSHOT_DIR:-$data_root/$nearline_name/slow_controls}
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
if [[ -n "${ARCUBE_NEARLINE_GAS_QUALITY_CONFIG:-}" ]]; then
    args+=(--gas-quality-config "$ARCUBE_NEARLINE_GAS_QUALITY_CONFIG")
fi
"$python_cmd" "$script_dir/lifetime.py" "${args[@]}" "$@" 2>&1 | tee "$logpath"
