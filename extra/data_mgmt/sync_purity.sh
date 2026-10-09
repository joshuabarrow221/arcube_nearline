#!/usr/bin/env bash
# Run on the receiving host. Immutable snapshots are exported on DAQ05 first.
# No SSH connection details or credentials are inferred by this wrapper.
set -euo pipefail
: "${ARCUBE_NEARLINE_PURITY_RSYNC_SOURCE:?Set to user@host:/path/to/purity_exports/}"
python_cmd=${ARCUBE_NEARLINE_PYTHON:-python3}
script_dir=$(cd -- "$(dirname "${BASH_SOURCE[0]}")" && pwd)
if [[ -n "${ARCUBE_NEARLINE_LIFETIME_CONFIG:-}" ]]; then
    # Use the same destination as watcher/worker, preventing a successful
    # transfer into an archive that the configured monitor never reads.
    destination=$("$python_cmd" - "$script_dir/../../actions/lifetime" "$ARCUBE_NEARLINE_LIFETIME_CONFIG" <<'PY'
import sys
sys.path.insert(0, sys.argv[1])
from monitor_config import load_config
print(load_config(sys.argv[2])['snapshot_root'])
PY
)
else
    : "${ARCUBE_NEARLINE_PURITY_SNAPSHOT_DIR:?Set the local snapshot directory}"
    destination=$ARCUBE_NEARLINE_PURITY_SNAPSHOT_DIR
fi
mkdir -p "$destination"
# Normal rsync temp-file + rename behavior keeps partial transfers out of the
# flat *.json archive. Never use --inplace or --delete for these snapshots.
# RSYNC_RSH can specify a site-approved ProxyJump when one is required.
rsync -av --exclude='.*' --include='slowcontrols-*.json' --exclude='*' \
    "$ARCUBE_NEARLINE_PURITY_RSYNC_SOURCE" "$destination/"
# The configured path refreshes all active six-hour pools and daily controls,
# even with no FLOW input. The legacy environment path remains supported.
bash "$script_dir/../../actions/lifetime.sh"
