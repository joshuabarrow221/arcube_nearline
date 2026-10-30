"""Offline SSH/export contracts, not evidence of NERSC network connectivity.

The fake SSH executable runs the transmitted bootstrap in another Python process
with the actual shipped exporter modules. Only SQLAlchemy's network/SQL execution
boundary is replaced. Thus SQL construction, source clocks, export filtering,
remote credential loading and local atomic publication remain production code.
"""
from argparse import Namespace
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys

import pytest


REPO = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location('purity_fetch_tested', REPO/'extra/data_mgmt/fetch_purity.py')
fetcher = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(fetcher)


FAKE_SSH = r'''
import contextlib
from datetime import datetime, timedelta, timezone
import io
import json
import os
from pathlib import Path
import sys
import time
import types

program = sys.stdin.read()
Path(os.environ['FAKE_SSH_ARGS']).write_text(json.dumps(sys.argv[1:]))
Path(os.environ['FAKE_SSH_INPUT']).write_text(program)
mode = os.environ.get('FAKE_SSH_MODE', 'success')
if mode == 'partial_failure':
    print('{"schema_version":', flush=True)
    sys.exit(255)
if mode == 'malformed':
    print('not JSON')
    sys.exit(0)
if mode == 'timeout':
    print('{', flush=True)
    time.sleep(10)
    sys.exit(0)

audit = []
class Connection:
    def __enter__(self): return self
    def __exit__(self, *args): pass
    def begin(self): return self
    def execute(self, query, params=None):
        audit.append({'query': query, 'params': params})
        if not query.startswith('SELECT'):
            return []
        if mode == 'database_failure':
            raise RuntimeError('SENTINEL_PRIVATE_PASSWORD must not leak')
        if mode == 'empty':
            return []
        if 'query_cryo_float' in query:
            # Ignition epochs are already UTC; the PRM correction must not
            # be applied to this second source.
            first = params['start'].replace(tzinfo=timezone.utc)
            return [(int(first.timestamp()*1000), 11.5),
                    (int(params['end'].replace(tzinfo=timezone.utc).timestamp()*1000), 999.)]
        first = params['start']
        return [(first, .00025), (first+timedelta(minutes=30), -.00001),
                (first+timedelta(minutes=45), float('nan')), (params['end'], .0008)]

class Engine:
    def connect(self): return Connection()
    def dispose(self): audit.append({'disposed': True})

def create_engine(url, **kwargs):
    # This assertion proves the remote bootstrap read its credential INI.
    assert url['password'] == 'SENTINEL_PRIVATE_PASSWORD'
    assert kwargs == {'connect_args': {'connect_timeout': 15}}
    audit.append({'connected': True})
    return Engine()

sqlalchemy = types.ModuleType('sqlalchemy')
sqlalchemy.__version__ = '2.0.0'
sqlalchemy.URL = types.SimpleNamespace(create=lambda driver, **values: values)
sqlalchemy.create_engine = create_engine
sqlalchemy.text = lambda value: value
sys.modules['sqlalchemy'] = sqlalchemy

# No repository on the simulated remote cwd or PYTHONPATH. The bootstrap must
# materialize/import the two modules transmitted over stdin itself.
os.chdir(os.environ['FAKE_SSH_REMOTE_ROOT'])
capture = io.StringIO()
try:
    with contextlib.redirect_stdout(capture):
        exec(compile(program, '<ssh-stdin>', 'exec'), {'__name__': '__main__'})
finally:
    Path(os.environ['FAKE_SSH_AUDIT']).write_text(json.dumps(audit, default=str))
payload = json.loads(capture.getvalue())
if mode == 'wrong_unit': payload['measurements'][0]['unit'] = 'ms'
if mode == 'wrong_window': payload['end'] = '2027-01-01T00:00:00+00:00'
if mode == 'wrong_hash': payload['export_code']['module_sha256']['purity_sources.py'] = 'wrong'
if mode == 'missing_source': payload['sources'].pop()
if mode == 'nonfinite': payload['measurements'][0]['value'] = float('nan')
if mode == 'numeric_timestamp': payload['measurements'][0]['timestamp'] = 1234567890
if mode == 'malformed_container': payload = []
print(json.dumps(payload))
'''


@pytest.fixture
def transport(tmp_path, monkeypatch):
    """Replace SSH through PATH only; real subprocess/streaming still executes."""
    binaries = tmp_path/'bin'
    binaries.mkdir()
    ssh = binaries/'ssh'
    ssh.write_text('#!' + sys.executable + '\n' + FAKE_SSH)
    ssh.chmod(0o755)
    remote = tmp_path/'remote without checkout'
    remote.mkdir()
    credentials = remote/'database.ini'
    credentials.write_text('[secrets]\nhost=private-db\nport=5432\ndatabase=purity\n'
                           'user=readonly\npassword=SENTINEL_PRIVATE_PASSWORD\n')
    source = tmp_path/'sources.json'
    source.write_text(json.dumps({'sources': [
        dict(name='PSQL PRM', kind='postgres', url_env='UNSET_PRIVATE_URL',
             table='prm_table', column='prm_lifetime', naive_timezone='America/Chicago',
             quantity='prm_lifetime', unit='s'),
        dict(name='O2', kind='ignition_function', url_env='UNSET_PRIVATE_URL',
             function='query_cryo_float', tag_id=1890, naive_timezone='UTC',
             quantity='o2', unit='ppb')]}))
    monkeypatch.setenv('PATH', str(binaries) + os.pathsep + os.environ['PATH'])
    monkeypatch.setenv('ARCUBE_NEARLINE_LOCAL_OUTPUT', '1')
    monkeypatch.delenv('PYTHONPATH', raising=False)
    for key, path in [('ARGS',tmp_path/'ssh-args.json'), ('INPUT',tmp_path/'ssh-input.py'),
                      ('AUDIT',tmp_path/'sql-audit.json'), ('REMOTE_ROOT',remote)]:
        monkeypatch.setenv('FAKE_SSH_' + key, str(path))
    args = Namespace(host='reader@daq.invalid', jump_host='reader@gateway.invalid',
                     remote_python='python3', source_config=source,
                     credential_ini=str(credentials), credential_section='secrets',
                     start='2026-10-30T13:00:00Z', end='2026-10-30T14:00:00Z',
                     snapshot_directory=tmp_path/'received', timeout_seconds=30)
    return args, tmp_path


@pytest.mark.parametrize('day, hour', [('2026-10-30', 13), ('2026-12-01', 14)])
def test_real_exporter_over_fake_ssh_corrects_prm_clock_and_retains_audit_values(transport, day, hour):
    args, root = transport
    args.start = f'{day}T{hour:02}:00:00Z'
    args.end = f'{day}T{hour+1:02}:00:00Z'
    output = fetcher.fetch(args)
    payload = json.loads(output.read_text())
    assert output.name.startswith('slowcontrols-') and output.suffix == '.json'
    assert list(args.snapshot_directory.iterdir()) == [output]
    prm = [r for r in payload['measurements'] if r['source'] == 'PSQL PRM']
    assert [r['value'] for r in prm] == [.00025, -.00001]
    assert prm[0]['timestamp'] == f'{day}T{hour:02}:00:00+00:00'
    assert prm[0]['source_naive_timezone'] == 'America/Chicago'
    assert prm[0]['unit'] == 's'
    oxygen = next(r for r in payload['measurements'] if r['source'] == 'O2')
    assert oxygen['timestamp'] == prm[0]['timestamp']
    assert oxygen['value'] == 11.5 and oxygen['unit'] == 'ppb'
    # NaN was omitted from the serialized rows, not from the query receipt.
    assert payload['sources'] == [{'source':'PSQL PRM','count':3}, {'source':'O2','count':1}]
    audit = json.loads((root/'sql-audit.json').read_text())
    queries = [r for r in audit if 'query' in r]
    for sql in ('SET TRANSACTION READ ONLY', "SET LOCAL timezone = 'UTC'",
                "SET LOCAL statement_timeout = '60s'"):
        assert sum(r['query'] == sql for r in queries) == 2
    prm_query = next(r for r in queries if 'SELECT timestamp' in r['query'])
    assert prm_query['params'] == {'start': f'{day} 08:00:00', 'end': f'{day} 09:00:00'}
    assert ':start' in prm_query['query'] and 'timestamp<:end' in prm_query['query']
    gas_query = next(r for r in queries if 'query_cryo_float' in r['query'])
    assert gas_query['params']['tag'] == 1890
    assert gas_query['params']['start'] == f'{day} {hour:02}:00:00'
    assert payload['export_code']['module_sha256'] == {
        name: hashlib.sha256((fetcher.SOURCE/name).read_text().encode()).hexdigest()
        for name in fetcher.MODULES}
    for content in (output.read_text(), (root/'ssh-input.py').read_text(), (root/'ssh-args.json').read_text()):
        assert 'SENTINEL_PRIVATE_PASSWORD' not in content
    assert list((root/'remote without checkout').iterdir()) == [Path(args.credential_ini)]
    command = json.loads((root/'ssh-args.json').read_text())
    assert 'StrictHostKeyChecking=yes' in command and 'BatchMode=yes' in command
    assert command[command.index('-J')+1] == args.jump_host
    assert '-L' not in command and '-N' not in command


def test_empty_successful_export_is_published_without_fabricated_measurements(transport, monkeypatch):
    args, _ = transport
    monkeypatch.setenv('FAKE_SSH_MODE', 'empty')
    payload = json.loads(fetcher.fetch(args).read_text())
    assert payload['measurements'] == []
    assert all(source['count'] == 0 for source in payload['sources'])


@pytest.mark.parametrize('mode', ['partial_failure', 'malformed', 'timeout',
                                 'wrong_unit', 'wrong_window', 'wrong_hash',
                                 'missing_source', 'nonfinite', 'database_failure',
                                 'numeric_timestamp', 'malformed_container'])
def test_failed_fetch_preserves_previous_snapshots_and_cleans_temporary_files(transport, monkeypatch, mode, capfd):
    args, _ = transport
    args.snapshot_directory.mkdir()
    previous = args.snapshot_directory/'slowcontrols-previous.json'
    previous.write_text('{"previous": "unchanged"}\n')
    old = previous.read_bytes()
    monkeypatch.setenv('FAKE_SSH_MODE', mode)
    if mode == 'timeout':
        args.timeout_seconds = .2
    with pytest.raises((ValueError, KeyError, subprocess.SubprocessError)):
        fetcher.fetch(args)
    assert previous.read_bytes() == old
    assert list(args.snapshot_directory.iterdir()) == [previous]
    assert 'SENTINEL_PRIVATE_PASSWORD' not in capfd.readouterr().err


def test_remote_interpreter_is_one_shell_quoted_executable(transport):
    args, _ = transport
    args.remote_python = '/remote/venv with spaces/py;$(touch unwanted)'
    command = fetcher.ssh_command(args)
    assert command[-3:-1] == ['--', args.host]
    assert shlex.split(command[-1]) == [args.remote_python, '-']


@pytest.mark.parametrize('start,end', [('2026-10-30T13:00:00', '2026-10-30T14:00:00Z'),
                                     ('2026-10-30T14:00:00Z', '2026-10-30T13:00:00Z')])
def test_invalid_query_interval_fails_before_any_ssh_or_publication(transport, start, end):
    args, root = transport
    args.start, args.end = start, end
    with pytest.raises(ValueError):
        fetcher.fetch(args)
    assert not (root/'ssh-args.json').exists()
    assert not args.snapshot_directory.exists()


@pytest.mark.parametrize('invalid', ['duplicate', 'nonobject'])
def test_invalid_source_configuration_fails_before_query(transport, invalid):
    args, root = transport
    config = json.loads(args.source_config.read_text())
    if invalid == 'duplicate':
        config['sources'].append(dict(config['sources'][0]))
    else:
        config['sources'] = [1]
    args.source_config.write_text(json.dumps(config))
    with pytest.raises(ValueError):
        fetcher.fetch(args)
    assert not (root/'ssh-args.json').exists()
    assert not args.snapshot_directory.exists()


def test_transport_client_help_needs_no_site_packages():
    # -S excludes venv site-packages: a transport-only NERSC client must not
    # import numpy/pandas/SQLAlchemy before SSH starts the remote exporter.
    result = subprocess.run([sys.executable, '-S', str(SPEC.origin), '--help'],
                            capture_output=True, text=True, check=True)
    assert '--snapshot-directory' in result.stdout
    assert '--credential-ini' in result.stdout
