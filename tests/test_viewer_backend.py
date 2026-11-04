"""Exercise CLI subprocesses with notebook-only plotting settings inherited."""
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

import h5py
import pytest

from test_track_viewer import make_receipt, POOL_KEY


SCRIPTS = Path(__file__).resolve().parents[1] / 'actions/lifetime'
INLINE_BACKEND = 'module://matplotlib_inline.backend_inline'


def export_viewer(tmp_path, environment):
    summary, receipt = make_receipt(tmp_path)
    output = tmp_path / 'viewer.html'
    payload = tmp_path / 'viewer.json'
    result = subprocess.run([
        sys.executable, str(SCRIPTS / 'track_viewer.py'), str(summary),
        '--output', str(output), '--json-output', str(payload),
    ], env=environment, capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stderr
    assert 'Verified 1/1:' in result.stderr
    assert output.is_file()
    saved = json.loads(payload.read_text())
    # Do not solve the backend failure by silently disabling the direct fit.
    assert saved['charge_scale']['direct']['n_segments'] == receipt['n_valid_segments']
    assert saved['charge_scale']['direct']['status'] == 'failed'  # Tiny fixture cannot fit.
    assert len(saved['tracks']) == 2


@pytest.mark.parametrize('backend', [INLINE_BACKEND, 'missing-notebook-backend'])
def test_viewer_cli_overrides_unusable_inherited_backend(tmp_path, backend):
    environment = dict(os.environ, MPLBACKEND=backend,
                       MPLCONFIGDIR=str(tmp_path / 'mpl-cache'))
    export_viewer(tmp_path, environment)


def test_notebook_child_environment_workaround_preserves_parent(tmp_path, monkeypatch):
    monkeypatch.setenv('MPLBACKEND', INLINE_BACKEND)
    environment = dict(os.environ, MPLBACKEND='Agg',
                       MPLCONFIGDIR=str(tmp_path / 'mpl-cache'))
    export_viewer(tmp_path, environment)
    assert os.environ['MPLBACKEND'] == INLINE_BACKEND


def test_charge_monitor_cli_direct_fit_with_notebook_backend(tmp_path):
    summary, receipt = make_receipt(tmp_path)
    with h5py.File(receipt['pooled_segments'], 'r') as pool:
        digest = hashlib.sha256(pool[POOL_KEY][:].tobytes()).hexdigest()
    receipt['fit'] = {'flow_datasets': [{'path': POOL_KEY, 'sha256': digest}]}
    receipt['reconstruction_signature'] = {'reconstruction': {'run_info': {'crs_ticks': .1}}}
    summary.write_text(json.dumps(receipt))
    output, payload = tmp_path / 'monitor.html', tmp_path / 'monitor.json'
    result = subprocess.run([
        sys.executable, str(SCRIPTS / 'charge_scale.py'), str(summary), '--direct',
        '--output', str(output), '--json-output', str(payload),
    ], env=dict(os.environ, MPLBACKEND=INLINE_BACKEND,
                MPLCONFIGDIR=str(tmp_path / 'mpl-cache')),
        capture_output=True, text=True, timeout=60)
    # The small sample has no accepted intercept; a complete diagnostic report
    # is still written. A backend import failure would leave no report at all.
    assert result.returncode == 2, result.stderr
    assert output.is_file(), result.stderr
    record = json.loads(payload.read_text())['records'][0]
    assert record['direct']['status'] == 'failed'
    assert record['direct']['n_segments'] == receipt['n_valid_segments']
    assert not record['direct'].get('input_error')
