"""Tests use synthetic fixtures only; they are never published as detector data."""
import importlib.util
import json
from pathlib import Path
import subprocess
import sys

import h5py
import numpy as np
import pandas as pd
from PIL import Image
import pytest

SCRIPTS = Path(__file__).resolve().parents[1] / 'actions/lifetime'
sys.path.insert(0, str(SCRIPTS))
import lifetime
from purity_sources import aggregate_measurements as aggregate_with_quality, load_measurements, aware_time, atomic_json
from packet_pedestal import flow_ped_map, key


def aggregate_measurements(paths):
    # Unit/conversion fixtures intentionally contain only a few observations.
    return aggregate_with_quality(paths, quality_config={'minimum_span_hours':0, 'max_gap_minutes':None})


def snapshot(tmp_path, rows):
    path = tmp_path/'snapshot.json'
    atomic_json(path, {'measurements': rows})
    return path


def row(quantity, value, timestamp, unit=None):
    return dict(quantity=quantity, value=value, timestamp=timestamp,
                unit=unit or ('s' if quantity == 'prm_lifetime' else 'ppb'), source=quantity)


def test_daily_prm_and_gas_concentrations_before_conversion(tmp_path):
    path = snapshot(tmp_path, [row('prm_lifetime', .001, '2026-09-29T01:00:00-05:00'),
        row('prm_lifetime', .003, '2026-09-29T12:00:00-05:00'),
        row('o2', .001, '2026-09-29T00:01:00-05:00', 'ppm'),
        row('o2', .003, '2026-09-29T05:59:00-05:00', 'ppm'),
        row('h2o', 17, '2026-09-29T03:00:00-05:00')])
    values = {r['sample']: r for r in aggregate_measurements([path, path])}
    assert values['prm']['lifetime_us'] == 2000
    assert values['prm']['error_us'] == 1000
    assert values['prm']['n_measurements'] == 2
    assert values['gas']['lifetime_us'] == pytest.approx(1000/(2/.299+1))
    assert values['gas_o2']['lifetime_us'] == pytest.approx(299/2)
    assert values['gas']['n_measurements'] == {'h2o': 1, 'o2': 2}


def test_missing_water_is_not_zero(tmp_path):
    values = aggregate_measurements([snapshot(tmp_path, [row('o2', 2, '2026-09-29T00:01:00Z')])])
    assert next(v for v in values if v['sample'] == 'gas')['lifetime_us'] is None
    assert next(v for v in values if v['sample'] == 'gas_o2')['lifetime_us'] == 149.5


def test_zero_and_invalid_values_do_not_become_finite_lifetimes(tmp_path):
    path = snapshot(tmp_path, [row('prm_lifetime', -1, '2026-09-29T00:00:00Z'),
        row('o2', 0, '2026-09-29T00:00:00Z'), row('h2o', 0, '2026-09-29T00:00:00Z')])
    result = aggregate_measurements([path])
    assert len(result) == 2 and all(r['fit_status'] == 'unavailable' for r in result)


def test_conflicting_duplicates_fail(tmp_path):
    with pytest.raises(ValueError, match='conflicting'):
        load_measurements([snapshot(tmp_path, [row('o2', 1, '2026-09-29T00:00:00Z'),
                                               row('o2', 2, '2026-09-29T00:00:00Z')])])


def test_naive_time_rejected():
    with pytest.raises(ValueError, match='UTC offset'):
        aware_time('2026-09-29')


@pytest.mark.parametrize('date,offset,hours', [('2026-03-08','-06:00',23), ('2026-11-01','-05:00',25)])
def test_local_day_and_six_hour_boundaries_across_dst(tmp_path, date, offset, hours):
    rows = [row('prm_lifetime', .001, date+'T00:30:00'+offset),
            row('o2', 1, date+'T00:30:00'+offset), row('h2o', 2, date+'T00:30:00'+offset)]
    values = {r['sample']: r for r in aggregate_measurements([snapshot(tmp_path, rows)])}
    prm = values['prm']; gas = values['gas']
    assert (aware_time(prm['period_end'])-aware_time(prm['period_start'])).total_seconds() == hours*3600
    assert pd.Timestamp(gas['period_end']).hour == 6


def test_retry_replaces_results_and_legacy_history_survives(tmp_path):
    path = tmp_path/'history.json'
    atomic_json(path, {'lifetimes':[dict(timestamp='2024-01-01T00:00:00Z', lifetime_us=900, error_us=20)]})
    entry = dict(timestamp='2026-09-29T00:00:00Z', input_file='x.h5', sample='beam',
                 method='track', fit_status='ok', lifetime_us=1000, error_us=100)
    lifetime.update_json(path, [entry])
    entry['lifetime_us'] = 1100
    lifetime.update_json(path, [entry])
    records = json.loads(path.read_text())['lifetimes']
    assert len(records) == 2 and records[-1]['lifetime_us'] == 1100


def test_pool_disjoint_track_samples(tmp_path):
    with h5py.File(tmp_path/'tracks.h5', 'w') as f:
        for _, name in lifetime.SPLIT_SAMPLES:
            f.create_dataset(name, data=np.zeros(1))
        samples = lifetime.find_samples(f)
        assert samples[-1] == ('all_mip', [p for _,p in lifetime.SPLIT_SAMPLES])


def test_static_pedestal_recovery_uses_references_and_channel_specific_slopes(tmp_path):
    path = tmp_path/'flow.h5'
    fields = [('io_group','u1'),('io_channel','u1'),('chip_id','u1'),('channel_id','u1')]
    hits = np.zeros(8, dtype=fields+[('Q_raw','f8')])
    packets = np.zeros(8, dtype=fields+[('packet_type','u1'),('dataword','u1')])
    for i in range(8):
        ch = i//4
        for obj in (hits, packets):
            obj[i]['io_group'], obj[i]['io_channel'], obj[i]['chip_id'], obj[i]['channel_id'] = 1, 2, 3, ch
        packets[i]['dataword'] = 50+i%4
        hits[i]['Q_raw'] = (0.5+ch)*(packets[i]['dataword']-(30+ch*5))
    # Intentionally permute packet storage: equal-index assumptions must fail this fixture.
    with h5py.File(path, 'w') as f:
        f.create_dataset('charge/calib_prompt_hits/data',data=hits)
        f.create_dataset('charge/packets/data',data=packets[::-1])
        f.create_dataset('charge/calib_prompt_hits/ref/charge/packets/ref',data=np.c_[np.arange(8),7-np.arange(8)])
    result = flow_ped_map(path, chunk=3)
    assert result[int(key(1,2,3,0))] == pytest.approx(30)
    assert result[int(key(1,2,3,1))] == pytest.approx(35)


def test_prm_only_cli_without_muon_selection(tmp_path):
    flow = tmp_path/'packet-1-2026_09_29_00_00_00_CDT.FLOW.hdf5'
    with h5py.File(flow, 'w'):
        pass
    raw = snapshot(tmp_path, [row('prm_lifetime', .001, '2026-09-29T00:00:00-05:00')])
    out = tmp_path/'overlay.png'
    result = subprocess.run([sys.executable, str(SCRIPTS/'lifetime.py'), '--input_file',str(flow),
        '--output_file_plot',str(tmp_path/'fit.png'), '--output_file_json',str(tmp_path/'history.json'),
        '--output-timeseries',str(out), '--slow-controls',str(raw)], capture_output=True, text=True)
    assert result.returncode == 0, result.stdout+result.stderr
    assert Image.open(out).size == (3000,2000)
    rows = json.loads((tmp_path/'history.json').read_text())['lifetimes']
    assert any(r['method']=='prm' and r['fit_status']=='ok' for r in rows)
    assert any(r['method']=='track' and r['fit_status']=='unavailable' for r in rows)
    assert Path(str(out)+'.html').exists()
    for published in (out, Path(str(out)+'.html'), tmp_path/'history.json'):
        assert published.stat().st_mode & 0o004, 'web server must be able to read published artifacts'


def test_common_packet_fit_recovers_attenuation():
    from packet_lifetime import common_fit_unweighted
    times = np.arange(25, 180, 10)
    data = pd.DataFrame({'io_group':0, 'time_us':times, 'mpv':20*np.exp(-times/1500)})
    alpha, _, _, _, _ = common_fit_unweighted(data,20,180)
    assert 1/alpha == pytest.approx(1.5,rel=1e-5)
    with pytest.raises(ValueError, match='Insufficient'):
        common_fit_unweighted(data.iloc[:2],20,180)


def test_auxiliary_channels_are_archived_but_not_combined(tmp_path):
    oxygen = row('o2', 1, '2026-09-29T00:00:00Z')
    auxiliary = dict(oxygen, value=100, source='second oxygen', use_for_lifetime=False)
    nitrogen = dict(oxygen, quantity='n2', value=20, unit='ppm', use_for_lifetime=False)
    values = aggregate_measurements([snapshot(tmp_path, [oxygen, auxiliary, nitrogen])])
    assert next(r for r in values if r['sample']=='gas_o2')['lifetime_us'] == 299


def test_export_failure_preserves_previous_snapshot(tmp_path, monkeypatch):
    import purity_sources as sources
    output = tmp_path/'existing.json'
    output.write_text('{"previous":"good"}')
    config = tmp_path/'config.json'
    config.write_text(json.dumps({'sources':[dict(kind='postgres',name='test')]}))
    def fail(*args):
        raise RuntimeError('connection failed')
    monkeypatch.setattr(sources, 'query_postgres', fail)
    with pytest.raises(RuntimeError, match='connection'):
        sources.export_snapshot(config, '2026-09-29T00:00:00Z', '2026-09-30T00:00:00Z', output)
    assert json.loads(output.read_text()) == {'previous':'good'}


def test_ignition_function_uses_bound_tag_and_half_open_window(monkeypatch):
    import sqlalchemy as sa
    from purity_sources import query_postgres
    executed = []
    class Connection:
        def __enter__(self): return self
        def __exit__(self, *args): pass
        def begin(self): return self
        def execute(self, sql, params=None):
            executed.append((str(sql), params))
            if str(sql).startswith('SELECT'):
                return [(pd.Timestamp('2026-09-29T00:00:00Z'),1),
                        (pd.Timestamp('2026-09-30T00:00:00Z'),2)]
            return []
    class Engine:
        def connect(self): return Connection()
        def dispose(self): pass
    monkeypatch.setenv('TEST_PSQL_URL', 'postgresql://localhost/test')
    monkeypatch.setattr(sa,'create_engine', lambda *args,**kwargs:Engine())
    values = query_postgres(dict(kind='ignition_function',tag_id=1893,url_env='TEST_PSQL_URL'),
                           aware_time('2026-09-29T00:00:00Z'),aware_time('2026-09-30T00:00:00Z'))
    assert len(values) == 1
    query, params = next((q,p) for q,p in executed if q.startswith('SELECT'))
    assert 'query_cryo_float' in query and ':tag' in query and params['tag'] == 1893
    assert any('READ ONLY' in q for q,_ in executed)


def test_packet_summary_units_and_missing_uncertainty(tmp_path):
    from argparse import Namespace
    path = tmp_path/'summary.json'
    atomic_json(path, dict(aggregation='detector-wide pooled charge', weighted_tau_ms=1.5,
                          packet_file='packet.h5', n_fit_objects=100))
    args = Namespace(timestamp='2026-09-29T00:00:00Z',input_file=None,packet_file=None,packet_summary=str(path))
    result = lifetime.packet_result(args)
    assert result['lifetime_us'] == 1500 and result['error_us'] is None
    assert result['fit_status'] == 'ok'


def test_plateau_withholds_lifetime_and_retains_reason(tmp_path):
    times = pd.date_range('2026-09-29T00:00:00-05:00', periods=36, freq='10min')
    rows = [row('o2', 1.0, t.isoformat()) for t in times]
    rows += [row('h2o', 20+i*.01, t.isoformat()) for i,t in enumerate(times)]
    result = aggregate_with_quality([snapshot(tmp_path, rows)])
    assert all(r['fit_status']=='unavailable' for r in result)
    assert all('plateau' in r['fit_message'] for r in result)
    assert all(r['lifetime_us'] is None for r in result)


def test_quality_tolerance_and_range_are_configurable(tmp_path):
    times = pd.date_range('2026-09-29T00:00:00-05:00',periods=36,freq='10min')
    rows = [row('o2',1+(i%2)*.001,t.isoformat()) for i,t in enumerate(times)]
    path=snapshot(tmp_path,rows)
    default=aggregate_with_quality([path])
    assert next(r for r in default if r['sample']=='gas_o2')['fit_status']=='ok'
    flat=aggregate_with_quality([path],quality_config={'flat_tolerance_ppb':.01})
    assert 'plateau' in flat[0]['fit_message']
    out_of_range=aggregate_with_quality([path],quality_config={'maximum_ppb':.9})
    assert 'calibration range' in out_of_range[0]['fit_message']


def test_short_history_and_gaps_withhold_gas(tmp_path):
    path=snapshot(tmp_path,[row('o2',1,'2026-09-29T00:00:00Z'),row('o2',2,'2026-09-29T01:00:00Z')])
    result=aggregate_with_quality([path])
    assert all(r['fit_status']=='unavailable' for r in result)
    assert 'coverage' in result[0]['fit_message']


def test_two_oxygen_channels_produce_separate_estimates(tmp_path):
    t='2026-09-29T00:00:00Z'
    path=snapshot(tmp_path,[row('o2',1,t),dict(row('o2',.008,t,'ppm'),source='O2 ppm')])
    result=aggregate_measurements([path])
    values={r['source']:r['lifetime_us'] for r in result if r['sample']=='gas_o2'}
    assert values=={'o2':299,'O2 ppm':299/8}


def test_concurrent_history_writers_preserve_all_rows(tmp_path):
    path=tmp_path/'concurrent.json'
    code='''import sys
sys.path.insert(0,sys.argv[1])
from lifetime import update_json
for i in range(4):
    update_json(sys.argv[2],[dict(timestamp='2026-09-29T00:00:00Z',method='track',sample=sys.argv[3],input_file=str(i),lifetime_us=1000,error_us=1)])
'''
    children=[subprocess.Popen([sys.executable,'-c',code,str(SCRIPTS),str(path),str(i)]) for i in range(3)]
    assert all(child.wait(timeout=30)==0 for child in children)
    assert len(json.loads(path.read_text())['lifetimes'])==12


def test_latest_observation_is_not_averaging_midpoint(tmp_path):
    values=aggregate_measurements([snapshot(tmp_path,[row('prm_lifetime',.001,'2026-09-29T01:00:00-05:00')])])
    assert values[0]['timestamp'] != values[0]['last_observed_at']
    assert aware_time(values[0]['last_observed_at']) == aware_time('2026-09-29T01:00:00-05:00')
