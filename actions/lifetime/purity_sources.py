"""Slow-controls snapshots and explicit, reproducible lifetime aggregation.

Snapshots contain raw observations, not pre-averaged values. Credentials are
resolved only at export time from environment variables and never serialized.
"""
import json
import math
import os
import re
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


def atomic_json(path, payload):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, name = tempfile.mkstemp(prefix=path.name + '.', suffix='.tmp', dir=path.parent)
    try:
        with os.fdopen(fd, 'w') as stream:
            json.dump(payload, stream, indent=2, allow_nan=False)
            stream.write('\n')
        # Snapshots/history contain measurements only and are published by nearline.
        os.chmod(name, 0o644)
        os.replace(name, path)
    finally:
        if os.path.exists(name):
            os.unlink(name)


def aware_time(value, naive_timezone=None):
    stamp = pd.Timestamp(value)
    if pd.isna(stamp):
        raise ValueError('missing timestamp')
    if stamp.tzinfo is None:
        if not naive_timezone:
            raise ValueError(f'timestamp needs a UTC offset: {value!r}')
        stamp = stamp.tz_localize(naive_timezone, ambiguous='raise', nonexistent='raise')
    return stamp.tz_convert('UTC')


def load_measurements(paths):
    """Normalize and deduplicate repeated/overlapping raw snapshot exports."""
    rows = {}
    for path in paths:
        with open(path) as stream:
            payload = json.load(stream)
        for row in payload['measurements']:
            if row.get('valid', True) is not True or not row.get('use_for_lifetime', True):
                continue
            quantity = row['quantity']
            if quantity not in ('prm_lifetime', 'o2', 'h2o'):
                raise ValueError(f'unknown quantity {quantity!r} in {path}')
            value = float(row['value'])
            if not math.isfinite(value) or value < 0 or (quantity == 'prm_lifetime' and value == 0):
                continue
            unit = row['unit']
            scales = {'s': 1e6, 'ms': 1e3, 'us': 1} if quantity == 'prm_lifetime' else {'ppb': 1, 'ppm': 1000}
            if unit not in scales:
                raise ValueError(f'invalid unit {unit!r} for {quantity}')
            stamp = aware_time(row['timestamp'])
            source = row['source']
            key = (source, quantity, stamp.isoformat())
            normalized = dict(timestamp=stamp, quantity=quantity, source=source,
                              value=value * scales[unit])
            if key in rows and rows[key]['value'] != normalized['value']:
                raise ValueError(f'conflicting duplicate measurement: {key}')
            rows[key] = normalized
    return pd.DataFrame(rows.values(), columns=['timestamp', 'quantity', 'source', 'value'])


def aggregate_measurements(paths, timezone_name='America/Chicago', quality_config=None):
    frame = load_measurements(paths)
    results = []
    if frame.empty:
        return results
    frame['timestamp'] = pd.to_datetime(frame.timestamp, utc=True).dt.tz_convert(timezone_name)
    for source, group in frame[frame.quantity == 'prm_lifetime'].groupby('source'):
        for start, values in group.set_index('timestamp').value.resample('1D'):
            if values.empty:
                continue
            end = start + pd.DateOffset(days=1)
            results.append(dict(timestamp=(start + (end-start)/2).isoformat(),
                period_start=start.isoformat(), period_end=end.isoformat(),
                sample='prm', method='prm', source=source, fit_status='ok',
                lifetime_us=float(values.mean()),
                error_us=float(values.sem()) if len(values) > 1 else None,
                uncertainty='SEM of observed daily lifetimes; excludes calibration systematics',
                n_measurements=len(values), first_observed_at=values.index.min().isoformat(),
                last_observed_at=values.index.max().isoformat(), averaging='arithmetic mean, local calendar day'))
    gases = frame[frame.quantity != 'prm_lifetime'].copy()
    if gases.empty:
        return results
    quality = dict(minimum_span_hours=3, plateau_hours=3, max_gap_minutes=30,
                   flat_tolerance_ppb=0, minimum_ppb=0, maximum_ppb=None)
    quality.update(quality_config or {})
    per_source = quality.pop('sources', {})
    gases['period'] = gases.timestamp.map(lambda t: (t.tz_localize(None).normalize() + pd.Timedelta(hours=(t.hour//6)*6)).tz_localize(timezone_name))
    bins = {}
    for (species, source), readings in gases.groupby(['quantity', 'source']):
        config = dict(quality, **per_source.get(source, {}))
        if config['minimum_span_hours'] < 0 or config['plateau_hours'] <= 0 or config['flat_tolerance_ppb'] < 0:
            raise ValueError('gas quality durations/tolerance must be nonnegative; plateau_hours must be positive')
        if config['max_gap_minutes'] is not None and config['max_gap_minutes'] <= 0:
            raise ValueError('max_gap_minutes must be positive or null')
        readings = readings.sort_values('timestamp').copy()
        readings['plateau'] = plateau_flags(readings, config)
        for period, group in readings.groupby('period'):
            reasons = []
            span = (group.timestamp.max()-group.timestamp.min()).total_seconds()/3600
            if span < config['minimum_span_hours']:
                reasons.append('insufficient multi-hour coverage')
            if config['max_gap_minutes'] is not None:
                gaps = group.timestamp.diff().dt.total_seconds().dropna()/60
                if len(gaps) and gaps.max() > config['max_gap_minutes']:
                    reasons.append('sampling gap exceeds quality limit')
            if group.plateau.any():
                reasons.append('sustained analyzer plateau')
            if config['minimum_ppb'] is not None and (group.value < config['minimum_ppb']).any():
                reasons.append('below configured calibration range')
            if config['maximum_ppb'] is not None and (group.value > config['maximum_ppb']).any():
                reasons.append('above configured calibration range')
            bins[(species, source, period)] = dict(mean=float(group.value.mean()), count=len(group),
                reasons=reasons, span_hours=span, quality_config=config,
                first_observed_at=group.timestamp.min().isoformat(),
                last_observed_at=group.timestamp.max().isoformat())
    oxygen_sources = sorted(gases.loc[gases.quantity == 'o2', 'source'].unique())
    water_sources = sorted(gases.loc[gases.quantity == 'h2o', 'source'].unique())
    if len(water_sources) > 1:
        raise ValueError('select one H2O source for each gas-equivalent comparison')
    for period in sorted(gases.period.unique()):
        end = (period.tz_localize(None) + pd.Timedelta(hours=6)).tz_localize(timezone_name)
        water = bins.get(('h2o', water_sources[0], period)) if water_sources else None
        for oxygen_source in oxygen_sources:
            oxygen = bins.get(('o2', oxygen_source, period))
            base = dict(timestamp=(period+(end-period)/2).isoformat(),
                period_start=period.isoformat(), period_end=end.isoformat(),
                source=oxygen_source, last_observed_at=oxygen['last_observed_at'] if oxygen else None, fit_status='unavailable', lifetime_us=None, error_us=None,
                uncertainty='not estimated; conversion and instrument systematics not included',
                fit_quality='provisional_conversion',
                averaging='six-hour concentration means, then conversion',
                concentration_ppb={}, n_measurements={}, analyzer_quality={})
            for species, item in [('o2', oxygen), ('h2o', water)]:
                if item:
                    base['concentration_ppb'][species] = item['mean']
                    base['n_measurements'][species] = item['count']
                    base['analyzer_quality'][species] = item
            o2_only = dict(base, sample='gas_o2', method='gas_o2',
                conversion='Eva O2 term only: tau_us = 299 / mean_O2_ppb')
            if oxygen and not oxygen['reasons'] and oxygen['mean'] > 0:
                o2_only.update(lifetime_us=299/oxygen['mean'], fit_status='ok')
            else:
                o2_only['fit_message'] = '; '.join(oxygen['reasons']) if oxygen and oxygen['reasons'] else 'O2 unavailable or nonpositive'
            results.append(o2_only)
            combined = dict(base, sample='gas', method='gas',
                conversion='Eva: tau_us = 1000 / (mean_O2_ppb/0.299 + mean_H2O_ppb/17)',
                water_source=water_sources[0] if water_sources else None)
            reasons = []
            for species, item in [('O2', oxygen), ('H2O', water)]:
                if item is None:
                    reasons.append(species+' missing')
                elif item['reasons']:
                    reasons.extend(species+': '+reason for reason in item['reasons'])
            if not reasons:
                combined['last_observed_at'] = min(oxygen['last_observed_at'], water['last_observed_at'], key=aware_time)
                rate = oxygen['mean']/.299 + water['mean']/17
                if rate > 0:
                    combined.update(lifetime_us=1000/rate, fit_status='ok')
                else:
                    reasons.append('zero concentrations: finite lifetime not determined')
            if reasons:
                combined['fit_message'] = '; '.join(reasons)
            results.append(combined)
    return results


def plateau_flags(readings, config):
    """Flag complete sustained flat runs, including their initial samples.

    Tolerance is in ppb after unit conversion; zero catches exactly repeated
    digital values. A gap resets the run. No forward filling is performed.
    """
    flags = np.zeros(len(readings), dtype=bool)
    if not len(readings):
        return flags
    times = readings.timestamp.to_list()
    values = readings.value.to_numpy()
    start = 0
    low = high = values[0]
    for i in range(1, len(values)+1):
        at_end = i == len(values)
        gap = not at_end and config['max_gap_minutes'] is not None and (times[i]-times[i-1]).total_seconds()/60 > config['max_gap_minutes']
        changed = not at_end and max(high, values[i])-min(low, values[i]) > config['flat_tolerance_ppb']
        if at_end or gap or changed:
            if i-start >= 2 and (times[i-1]-times[start]).total_seconds()/3600 >= config['plateau_hours']:
                flags[start:i] = True
            if not at_end:
                start = i
                low = high = values[i]
        else:
            low, high = min(low, values[i]), max(high, values[i])
    return flags


def identifier(value):
    if not re.fullmatch(r'[A-Za-z_][A-Za-z_0-9]*', value):
        raise ValueError(f'invalid database identifier {value!r}')
    return value


def query_postgres(config, start, end):
    import sqlalchemy as sa
    if 'credential_ini' in config:
        import configparser
        credentials = configparser.ConfigParser()
        if not credentials.read(config['credential_ini']):
            raise FileNotFoundError('credential INI is not readable on this host')
        c = credentials[config.get('credential_section', 'secrets')]
        url = sa.URL.create('postgresql+psycopg2', username=c['user'], password=c['password'],
                            host=c['host'], port=int(c['port']), database=c['database'])
    else:
        url = os.environ[config['url_env']]
    engine = sa.create_engine(url, connect_args={'connect_timeout': 15})
    names = []
    if config['kind'] == 'ignition':
        month = start.tz_convert('UTC').tz_localize(None).to_period('M')
        last = (end-pd.Timedelta(nanoseconds=1)).tz_convert('UTC').tz_localize(None).to_period('M')
        while month <= last:
            names.append(f"{identifier(config['table_prefix'])}_{month.year}_{month.month:02d}")
            month += 1
    else:
        names = [identifier(config.get('table', 'prm_table'))]
    ranges = []
    cursor = start
    step = pd.Timedelta(hours=config.get('query_hours', 6))
    if step <= pd.Timedelta(0):
        raise ValueError('query_hours must be positive')
    while cursor < end:
        stop = min(cursor+step, end) if config['kind'] == 'ignition_function' else end
        ranges.append((cursor, stop))
        cursor = stop
    observations = []
    from itertools import product
    try:
        with engine.connect() as connection:
            with connection.begin():
                connection.execute(sa.text('SET TRANSACTION READ ONLY'))
                connection.execute(sa.text("SET LOCAL timezone = 'UTC'"))
                connection.execute(sa.text("SET LOCAL statement_timeout = '60s'"))
                for name, (query_start, query_end) in product(names, ranges):
                    if config['kind'] == 'ignition_function':
                        function = identifier(config.get('function', 'query_cryo_float'))
                        query = sa.text(f'SELECT t_stamp, floatvalue FROM {function}'
                                        '(CAST(:start AS timestamp(3)), CAST(:end AS timestamp(3)), :tag) ORDER BY t_stamp')
                        params = dict(tag=config['tag_id'], start=query_start.tz_convert('UTC').tz_localize(None).to_pydatetime(),
                                      end=query_end.tz_convert('UTC').tz_localize(None).to_pydatetime())
                    elif config['kind'] == 'ignition':
                        query = sa.text(f'SELECT t_stamp, floatvalue FROM "{name}" WHERE tagid=:tag '
                                        'AND t_stamp>=:start AND t_stamp<:end ORDER BY t_stamp')
                        params = dict(tag=config['tag_id'], start=int(start.timestamp()*1000), end=int(end.timestamp()*1000))
                    else:
                        column = identifier(config.get('column', 'prm_lifetime'))
                        query = sa.text(f'SELECT timestamp, "{column}" FROM "{name}" '
                                        'WHERE timestamp>=:start AND timestamp<:end ORDER BY timestamp')
                        params = dict(start=start.to_pydatetime(), end=end.to_pydatetime())
                    for stamp, value in connection.execute(query, params):
                        if value is None:
                            continue
                        if config['kind'] == 'ignition':
                            stamp = pd.to_datetime(stamp, unit='ms', utc=True)
                        elif config['kind'] == 'ignition_function' and isinstance(stamp, (int, float)):
                            stamp = pd.to_datetime(stamp, unit='ms', utc=True)
                        else:
                            stamp = aware_time(stamp, config.get('naive_timezone'))
                        if query_start <= stamp < query_end:
                            observations.append((stamp, value))
    finally:
        engine.dispose()
    return observations


def query_influx(config, start, end):
    from influxdb import InfluxDBClient
    kwargs = dict(config.get('connection', {}))
    for key, env in config.get('connection_env', {}).items():
        kwargs[key] = os.environ[env]
    client = InfluxDBClient(**kwargs, timeout=30)
    measurement = identifier(config['measurement']) if '-' not in config['measurement'] else config['measurement']
    if not re.fullmatch(r'[A-Za-z_0-9-]+', measurement):
        raise ValueError('invalid influx measurement')
    field = identifier(config['field'])
    try:
        query = (f'SELECT "{field}" FROM "{measurement}" WHERE '
                 f"time >= '{start.isoformat()}' AND time < '{end.isoformat()}'")
        result = client.query(query, database=config['database'])
        return [(aware_time(p['time']), p[field]) for p in result.get_points() if p.get(field) is not None]
    finally:
        client.close()


def export_snapshot(config_path, start, end, output):
    start, end = aware_time(start), aware_time(end)
    if end <= start:
        raise ValueError('end must be after start')
    config = json.loads(Path(config_path).read_text())
    rows, status = [], []
    for source in config['sources']:
        if not source.get('enabled', True):
            continue
        # Fail the export rather than atomically replacing a good snapshot with a partial one.
        values = (query_influx(source, start, end) if source['kind'] == 'influx'
                  else query_postgres(source, start, end))
        for stamp, value in values:
            value = float(value)
            if math.isfinite(value):
                rows.append(dict(timestamp=aware_time(stamp).isoformat(), quantity=source['quantity'],
                                 value=value, unit=source['unit'], source=source['name'],
                                 use_for_lifetime=source.get('use_for_lifetime', True)))
        status.append(dict(source=source['name'], count=len(values)))
    payload = dict(schema_version=1, exported_at=datetime.now(timezone.utc).isoformat(),
                   start=start.isoformat(), end=end.isoformat(), sources=status, measurements=rows)
    atomic_json(output, payload)
    return payload
