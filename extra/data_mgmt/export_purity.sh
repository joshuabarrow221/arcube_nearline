#!/usr/bin/env bash
# Run explicitly on a database-reachable host (for example DAQ05). Database
# credentials are environment variables referenced by the source config; this
# wrapper never writes them into a snapshot, command argument, or log.
set -euo pipefail
if [[ $# != 4 ]]; then
    echo 'Usage: export_purity.sh SOURCE_CONFIG START_ISO END_ISO EXPORT_DIRECTORY' >&2
    exit 2
fi
config=$1; start=$2; end=$3; destination=$4
mkdir -p "$destination"
# Unique filenames permit overlapping exports to be deduplicated at import.
# The producer writes under a hidden work directory, then moves the completed
# snapshot into the flat archive. A failed query can never replace good data.
attempt=$(mktemp -d "$destination/.export-XXXXXXXX")
trap 'rm -rf -- "$attempt"' EXIT
python_cmd=${ARCUBE_NEARLINE_PYTHON:-python3}
"$python_cmd" "$(dirname "${BASH_SOURCE[0]}")/../../actions/lifetime/lifetime.py" \
    --source-config "$config" --start "$start" --end "$end" \
    --export-slow-controls "$attempt/snapshot.json"
name="slowcontrols-$(date -u +%Y%m%dT%H%M%SZ)-${attempt##*.export-}.json"
mv -- "$attempt/snapshot.json" "$destination/$name"
echo "Exported $destination/$name"
