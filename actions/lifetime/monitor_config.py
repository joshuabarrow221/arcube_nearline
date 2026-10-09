"""Portable deployment configuration and exact mirrored input pairing.

This module uses only the standard library so watcher discovery can run without
loading the numerical fitting environment. No credentials belong in this file.
A config is explicit about input, private state, and publication directories;
there is no automatic write to the production website when testing locally.
"""
import hashlib
import json
import os
from pathlib import Path
from zoneinfo import ZoneInfo


def load_config(path):
    """Resolve relative paths against the config, then reject recursive layouts."""
    path = Path(path).resolve()
    config = json.loads(path.read_text())
    required = ('flow_root', 'packet_root', 'snapshot_root', 'state_root', 'publish_root')
    for key in required:
        value = os.path.expandvars(config[key])
        if '$' in value:
            raise ValueError(f'unresolved environment variable in {key}')
        p = Path(value).expanduser()
        config[key] = str((path.parent/p).resolve() if not p.is_absolute() else p.resolve())
    for key in ('gas_conversion_config', 'gas_quality_config', 'calibration_root'):
        if config.get(key):
            value = os.path.expandvars(config[key])
            if '$' in value:
                raise ValueError(f'unresolved environment variable in {key}')
            p = Path(value).expanduser()
            config[key] = str((path.parent/p).resolve() if not p.is_absolute() else p.resolve())
    watched = [Path(config[k]) for k in ('flow_root', 'packet_root', 'snapshot_root')]
    for key in ('state_root', 'publish_root'):
        p = Path(config[key])
        if any(p == root or root in p.parents or p in root.parents for root in watched):
            raise ValueError(f'{key} must be outside the watched input trees')
    state, public = Path(config['state_root']), Path(config['publish_root'])
    if state == public or state in public.parents or public in state.parents:
        raise ValueError('private state and public plots must not contain each other')
    config.setdefault('timezone', 'America/Chicago')
    ZoneInfo(config['timezone'])  # Fail at startup for an invalid timezone.
    # Window boundaries and display clock are separate decisions. Shifters
    # can compare UTC across detector systems without moving a PRM day or pool.
    config.setdefault('display_timezone', config['timezone'])
    ZoneInfo(config['display_timezone'])
    config.setdefault('window_hours', 6)
    if config['window_hours'] != 6:
        raise ValueError('this deployment profile uses fixed six-hour detector windows')
    config.setdefault('bootstrap', 100)
    if config['bootstrap'] < 10:
        raise ValueError('packet bootstrap count must be at least 10')
    config.setdefault('packet_enabled', True)
    config.setdefault('packet_ext_io', 6)
    config.setdefault('snapshot_pattern', 'slowcontrols-*.json')
    config.setdefault('display_days', 35)
    config.setdefault('maximum_display_us', 10000)
    if config['display_days'] <= 0 or config['maximum_display_us'] <= 0:
        raise ValueError('display period and lifetime range must be positive')
    config.setdefault('diagnostic_patterns', ['*HotPixelsHunt*', '*BadNoise*', '*induced_noise*'])
    config.setdefault('cohort_rules', [])
    config['_path'] = str(path)
    # File-content versions allow queued tasks to notice a configuration edit.
    # Included auxiliary configs must also affect a discovery revision.
    blobs = [path.read_bytes()]
    for key in ('gas_conversion_config', 'gas_quality_config'):
        if config.get(key): blobs.append(Path(config[key]).read_bytes())
    # A deliberate software update also makes discovery revisit existing
    # inputs. Extraction receipts use the narrower processing revision below,
    # so a plot-only edit does not repeat packet decoding.
    for name in ('monitor.py', 'track_shards.py', 'packet_pooling.py', 'packet_lifetime.py',
                 'packet_pedestal.py', 'calibration_pedestal.py', 'lifetime.py',
                 'lifetime_funcs.py', 'purity_sources.py'):
        blobs.append(Path(__file__).with_name(name).read_bytes())
    config['_revision'] = hashlib.sha256(b'\n'.join(blobs)).hexdigest()
    processing = {k: config.get(k) for k in ('flow_root', 'packet_root', 'calibration_root',
        'packet_enabled', 'packet_ext_io', 'diagnostic_patterns', 'cohort_rules', 'reprocess_token')}
    config['_processing_revision'] = hashlib.sha256(json.dumps(processing, sort_keys=True).encode()).hexdigest()
    return config


def file_version(path):
    """Identity of an immutable completed file, without hashing multi-GB inputs.

    Size and nanosecond mtime detect normal replacements. In-place edits that
    preserve both are unsupported: operators must publish immutable inputs or
    deliberately change the config/reprocess token. Scientific shard bytes are
    independently hash-checked before fitting.
    """
    p = Path(path)
    try:
        stat = p.stat()
    except FileNotFoundError:
        return None
    return dict(path=str(p.resolve()), size=stat.st_size, mtime_ns=stat.st_mtime_ns)


def paired_input(path, config):
    """Map either side of one acquisition to a canonical FLOW path.

    Preserve the full relative hierarchy. Equal basenames from different runs
    are distinct inputs, and partial/segment-cache files never qualify.
    """
    path = Path(path).resolve()
    flow, packet = Path(config['flow_root']), Path(config['packet_root'])
    if flow in path.parents and path.name.endswith('.FLOW.hdf5'):
        rel = path.relative_to(flow)
        native = packet/rel.parent/(rel.name.removesuffix('.FLOW.hdf5')+'.h5')
        return path, native, str(rel)
    if packet in path.parents and path.name.endswith('.h5') and path.name.startswith('packet-'):
        rel = path.relative_to(packet)
        primary = flow/rel.parent/(rel.name.removesuffix('.h5')+'.FLOW.hdf5')
        return primary, path, str(primary.relative_to(flow))
    return None


def discovery(path, config):
    """Return canonical action input, dependency revision, and stability inputs.

    A packet arriving after FLOW changes the same acquisition revision and
    schedules one retry. Packet-first arrival waits for FLOW. Slow-control
    snapshots independently trigger publication even with no detector files.
    """
    pair = paired_input(path, config)
    if pair:
        flow, packet, relative = pair
        fv, pv = file_version(flow), file_version(packet)
        if fv is None: return None
        dependencies = [fv] + ([pv] if pv else [])
        key = 'flow:'+relative
        action_input = flow
    else:
        path = Path(path).resolve()
        root = Path(config['snapshot_root'])
        if path.parent != root or not path.match(config['snapshot_pattern']): return None
        version = file_version(path)
        if version is None: return None
        key = 'snapshot:'+path.name; action_input = path; dependencies = [version]
    revision = hashlib.sha256(json.dumps([config['_revision'], dependencies], sort_keys=True).encode()).hexdigest()
    namespace = hashlib.sha256(config['_path'].encode()).hexdigest()[:16]
    return dict(key=namespace+':'+key, revision=revision, input=str(action_input), dependencies=dependencies)
