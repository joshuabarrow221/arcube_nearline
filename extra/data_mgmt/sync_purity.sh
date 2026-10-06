#!/usr/bin/env bash
# Run on the receiving host. Export immutable timestamped snapshots on DAQ05 first.
set -euo pipefail
: "${ARCUBE_NEARLINE_PURITY_RSYNC_SOURCE:?Set to user@host:/path/to/purity_exports/}"
: "${ARCUBE_NEARLINE_PURITY_SNAPSHOT_DIR:?Set the local snapshot directory}"
mkdir -p "$ARCUBE_NEARLINE_PURITY_SNAPSHOT_DIR"
# Normal rsync temp-file + rename behavior keeps partial transfers out of *.json.
rsync -av --include='*/' --include='*.json' --exclude='*' \
    "$ARCUBE_NEARLINE_PURITY_RSYNC_SOURCE" "$ARCUBE_NEARLINE_PURITY_SNAPSHOT_DIR/"
bash "$(dirname "${BASH_SOURCE[0]}")/../../actions/lifetime.sh"
