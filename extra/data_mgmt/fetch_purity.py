#!/usr/bin/env python3
"""Fetch completed purity snapshots through SSH from NERSC or another receiver.

The client needs only Python's standard library and OpenSSH. It sends the two
repository exporter modules to the database host, where the existing read-only
query implementation runs. Credentials are resolved there, never on the client.
No remote repository checkout, persistent tunnel, rsync, or detector fitting is
needed. Only an approved remote Python environment and database access are used.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import shlex
import subprocess
import sys
import tempfile
import uuid


SOURCE = Path(__file__).resolve().parents[2] / 'actions' / 'lifetime'
# Keep this explicit: never send a repository tree, environment, or credential
# file to the remote host. Both modules are shared with the ordinary exporter.
MODULES = ('purity_sources.py', 'lifetime_io.py')

REMOTE_PROGRAM = r'''
import contextlib
import hashlib
import json
from pathlib import Path
import sys
import tempfile

try:
    # The temporary directory lives only for this export. No services, checkout,
    # system configuration, or acquisition files are changed on the DAQ host.
    with tempfile.TemporaryDirectory(prefix='nearline-purity-') as directory:
        root = Path(directory)
        hashes = {}
        for name, source in request['modules'].items():
            (root / name).write_text(source, encoding='utf-8')
            hashes[name] = hashlib.sha256(source.encode('utf-8')).hexdigest()
        sys.path.insert(0, directory)
        config = root / 'sources.json'
        config.write_text(json.dumps(request['config']), encoding='utf-8')
        # Import the same implementation used by lifetime.py. Do not copy SQL
        # or clock-conversion logic into this transport helper.
        with contextlib.redirect_stdout(sys.stderr):
            from purity_sources import export_snapshot
            payload = export_snapshot(config, request['start'], request['end'],
                                      root / 'snapshot.json')
        payload['export_code'] = {'module_sha256': hashes}
        print(json.dumps(payload, allow_nan=False), flush=True)
except Exception as error:
    # Connection exceptions can contain credentials supplied by a driver.
    # Expose the exception type, not its message or a public traceback. Diagnose
    # further in a protected terminal on the export host when necessary.
    print('Purity export failed (' + type(error).__name__ + '). Check remote '
          'Python dependencies, credential configuration, database access and '
          'source timestamps on the export host.', file=sys.stderr)
    sys.exit(1)
'''


def timestamp(value):
    """Require an explicit offset; never infer the receiver's local clock."""
    if not isinstance(value, str):
        raise ValueError('timestamps must be ISO strings')
    stamp = datetime.fromisoformat(value.replace('Z', '+00:00'))
    if stamp.tzinfo is None or stamp.utcoffset() is None:
        raise ValueError('timestamps must include Z or an explicit UTC offset')
    return stamp.astimezone(timezone.utc)


def make_request(args):
    """Collect non-secret configuration and exactly the reviewed exporter code."""
    start, end = timestamp(args.start), timestamp(args.end)
    if end <= start:
        raise ValueError('end must be after start')
    config = json.loads(args.source_config.read_text())
    if not isinstance(config, dict):
        raise ValueError('source configuration must be a JSON object')
    sources = config['sources']
    if not isinstance(sources, list) or not all(isinstance(s, dict) for s in sources):
        raise ValueError('sources must be a list of source objects')
    if not any(s.get('enabled', True) for s in sources):
        raise ValueError('source configuration must enable at least one source')
    names = [s['name'] for s in sources if s.get('enabled', True)]
    if len(names) != len(set(names)):
        raise ValueError('enabled source names must be unique')
    if args.credential_ini:
        if not args.credential_ini.startswith('/'):
            raise ValueError('--credential-ini must name an absolute path on the remote host')
        # The known DAQ05 PRM/cryogenics queries share a database configuration.
        # For different databases, omit this override and set credential_ini or
        # url_env separately in each source. Influx configuration is unaffected.
        for source in sources:
            if source['kind'] != 'influx':
                source['credential_ini'] = args.credential_ini
                source['credential_section'] = args.credential_section
                source.pop('url_env', None)
    modules = {name: (SOURCE / name).read_text(encoding='utf-8') for name in MODULES}
    return dict(start=start.isoformat(), end=end.isoformat(), config=config, modules=modules)


def ssh_command(args):
    """Quote the remote shell command separately from the local argument list."""
    if not args.host or args.host.startswith('-') or any(c.isspace() for c in args.host):
        raise ValueError('--host must be a host or user@host, not SSH options')
    if not args.remote_python or args.remote_python.startswith('-'):
        raise ValueError('--remote-python must name an interpreter or executable wrapper')
    # These options apply to the destination. ProxyJump gateways also need
    # BatchMode/host-key policy in SSH config; see lifetime-remote-export.md.
    command = ['ssh', '-T', '-o', 'BatchMode=yes', '-o', 'ConnectTimeout=15',
               '-o', 'StrictHostKeyChecking=yes', '-o', 'ServerAliveInterval=15',
               '-o', 'ServerAliveCountMax=3']
    if args.jump_host:
        command += ['-J', args.jump_host]
    # SSH passes its command through the remote shell even with shell=False on
    # this machine. shlex.join protects interpreter paths containing spaces.
    return command + ['--', args.host, shlex.join([args.remote_python, '-'])]


def validate_snapshot(payload, request):
    """Reject truncated, wrong-window, or malformed responses before publication.

    Empty successful queries and finite negative instrument values are retained;
    scientific validity cuts belong to the established importer, not transport.
    Query counts can exceed serialized counts because nonfinite values are omitted.
    """
    start, end = timestamp(request['start']), timestamp(request['end'])
    if not isinstance(payload, dict):
        raise ValueError('snapshot must be a JSON object')
    if payload.get('schema_version') != 1:
        raise ValueError('unsupported snapshot schema')
    if timestamp(payload['start']) != start or timestamp(payload['end']) != end:
        raise ValueError('snapshot query interval does not match request')
    timestamp(payload['exported_at'])
    expected = {s['name']: s for s in request['config']['sources'] if s.get('enabled', True)}
    statuses = payload['sources']
    if not isinstance(statuses, list) or not all(isinstance(s, dict) for s in statuses):
        raise ValueError('snapshot sources must be a list of status objects')
    if len(statuses) != len(expected) or {s['source'] for s in statuses} != set(expected):
        raise ValueError('snapshot does not report every enabled source')
    if any(type(s['count']) is not int or s['count'] < 0 for s in statuses):
        raise ValueError('invalid source query count')
    hashes = {name: hashlib.sha256(source.encode('utf-8')).hexdigest()
              for name, source in request['modules'].items()}
    receipt = payload.get('export_code')
    if not isinstance(receipt, dict) or receipt.get('module_sha256') != hashes:
        raise ValueError('exporter code receipt does not match request')
    if not isinstance(payload['measurements'], list):
        raise ValueError('snapshot measurements must be a list')
    for row in payload['measurements']:
        if not isinstance(row, dict):
            raise ValueError('snapshot measurements must be objects')
        source = expected[row['source']]
        if row['quantity'] != source['quantity'] or row['unit'] != source['unit']:
            raise ValueError('snapshot quantity/unit does not match source configuration')
        if row.get('use_for_lifetime', True) != source.get('use_for_lifetime', True):
            raise ValueError('snapshot use_for_lifetime flag does not match configuration')
        if type(row['value']) not in (int, float) or not math.isfinite(row['value']):
            raise ValueError('snapshot contains a nonfinite or nonnumeric measurement')
        if not start <= timestamp(row['timestamp']) < end:
            raise ValueError('measurement falls outside requested interval')


def fetch(args):
    """Publish one immutable snapshot only after SSH and validation both succeed."""
    request = make_request(args)
    command = ssh_command(args)
    # repr of the JSON string makes configuration data a string literal in the
    # transmitted program, never executable Python or shell interpolation.
    program = 'import json\nrequest = json.loads(' + repr(json.dumps(request)) + ')\n' + REMOTE_PROGRAM
    args.snapshot_directory.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix='.fetch-', suffix='.tmp', dir=args.snapshot_directory)
    try:
        with os.fdopen(fd, 'w') as stream:
            subprocess.run(command, input=program, text=True, stdout=stream,
                           check=True, timeout=args.timeout_seconds)
        with open(temporary) as stream:
            # Python's JSON parser normally accepts nonstandard NaN/Infinity.
            # Reject those anywhere in the response, including metadata.
            def invalid_constant(value):
                raise ValueError('nonfinite JSON constant: ' + value)
            payload = json.load(stream, parse_constant=invalid_constant)
        validate_snapshot(payload, request)
        sys.path.insert(0, str(SOURCE))
        from lifetime_io import public_permissions
        public_permissions(temporary)
        name = 'slowcontrols-' + datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
        output = args.snapshot_directory / (name + '-' + uuid.uuid4().hex + '.json')
        os.replace(temporary, output)
        return output
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--host', required=True, help='Database-reachable SSH destination, user@host or SSH alias')
    parser.add_argument('--jump-host', help='Optional SSH -J gateway; SSH config is also honored')
    parser.add_argument('--remote-python', default='python3', help='Remote Python 3.10+ interpreter or absolute venv executable')
    parser.add_argument('--source-config', type=Path, required=True, help='Local non-secret source JSON; credentials resolve remotely')
    parser.add_argument('--credential-ini', help='Existing remote INI path, used for every PostgreSQL source')
    parser.add_argument('--credential-section', default='secrets')
    parser.add_argument('--start', required=True, help='Inclusive ISO timestamp with UTC offset')
    parser.add_argument('--end', required=True, help='Exclusive ISO timestamp with UTC offset')
    parser.add_argument('--snapshot-directory', type=Path, required=True, help='Receiver directory watched as snapshot_root')
    parser.add_argument('--timeout-seconds', type=int, default=600, help='Whole SSH/export timeout; default 600 seconds')
    args = parser.parse_args()
    if args.timeout_seconds <= 0:
        parser.error('--timeout-seconds must be positive')
    try:
        output = fetch(args)
    except (ValueError, KeyError, TypeError, OSError, subprocess.SubprocessError) as error:
        # Do not turn a failed export into a fresh-looking file or dump source
        # configuration into an operational log.
        print('No snapshot published (' + type(error).__name__ + '). Check SSH '
              'access, the source configuration and remote exporter diagnostics.', file=sys.stderr)
        return 1
    print('Exported ' + str(output))
    return 0


if __name__ == '__main__':
    sys.exit(main())
