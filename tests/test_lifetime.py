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


def test_retry_replaces_results_and_legacy_history_is_reference_only(tmp_path):
    path = tmp_path/'history.json'
    atomic_json(path, {'lifetimes':[dict(timestamp='2024-01-01T00:00:00Z', lifetime_us=900, error_us=20)]})
    entry = dict(timestamp='2026-09-29T00:00:00Z', input_file='x.h5', sample='beam',
                 method='track', calculation_source='flow_segments', fit_status='ok', lifetime_us=1000, error_us=100)
    lifetime.update_json(path, [entry])
    entry['lifetime_us'] = 1100
    lifetime.update_json(path, [entry])
    records = json.loads(path.read_text())['lifetimes']
    assert len(records) == 1 and records[0]['lifetime_us'] == 1100
    data=json.loads(path.read_text())
    assert len(data['reference_lifetimes']) == 1
    assert all(r['lifetime_us'] != 900 for r in data['plot_lifetimes'])


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
    update_json(sys.argv[2],[dict(timestamp='2026-09-29T00:00:00Z',method='track',calculation_source='flow_segments',sample=sys.argv[3],input_file=str(i),lifetime_us=1000,error_us=1)])
'''
    children=[subprocess.Popen([sys.executable,'-c',code,str(SCRIPTS),str(path),str(i)]) for i in range(3)]
    assert all(child.wait(timeout=30)==0 for child in children)
    assert len(json.loads(path.read_text())['lifetimes'])==12


def test_latest_observation_is_not_averaging_midpoint(tmp_path):
    values=aggregate_measurements([snapshot(tmp_path,[row('prm_lifetime',.001,'2026-09-29T01:00:00-05:00')])])
    assert values[0]['timestamp'] != values[0]['last_observed_at']
    assert aware_time(values[0]['last_observed_at']) == aware_time('2026-09-29T01:00:00-05:00')


def test_six_hour_track_means_preserve_samples_and_raw_history(tmp_path):
    def entry(file, value, error, sample='mixed', timestamp='2026-10-02T12:10:00-05:00'):
        return dict(input_file=file, timestamp=timestamp, lifetime_us=value,
                    error_us=error, sample=sample, method='track', calculation_source='flow_segments', fit_status='ok')
    rows=[entry('a',1000,30),entry('b',1400,40),entry('beam',1700,60,'beam'),
          dict(entry('failed',None,None),fit_status='failed')]
    path=tmp_path/'history.json'
    lifetime.update_json(path,rows)
    data=json.loads(path.read_text())
    assert data['lifetimes']==rows
    mixed=next(r for r in data['plot_lifetimes'] if r['sample']=='mixed')
    assert mixed['lifetime_us']==1200 and mixed['n_measurements']==2
    assert mixed['n_rejected']==1
    assert mixed['propagated_fit_error_us']==25
    assert mixed['between_file_sem_us']==pytest.approx(200)
    assert mixed['error_us']==pytest.approx(200)
    assert len(mixed['members'])==3
    assert mixed['period_start']=='2026-10-02T12:00:00-05:00'
    assert mixed['timestamp']=='2026-10-02T15:00:00-05:00'
    assert mixed['last_observed_at']=='2026-10-02T12:10:00-05:00'
    assert len(lifetime.aggregate_track_lifetimes(rows+[rows[0]]))==2


@pytest.mark.parametrize('date,offset,hours',[('2026-03-08','-06:00',5),('2026-11-01','-05:00',7)])
def test_track_six_hour_windows_follow_local_dst(date,offset,hours):
    rows=lifetime.aggregate_track_lifetimes([dict(timestamp=date+'T00:30:00'+offset,lifetime_us=1000,error_us=None,input_file='flow.h5',calculation_source='flow_segments')])
    assert (aware_time(rows[0]['period_end'])-aware_time(rows[0]['period_start'])).total_seconds()==hours*3600
    assert rows[0]['error_us'] is None


def test_lines_connect_valid_points_across_rejected_and_missing_windows():
    def gas(hour,value):
        return dict(timestamp=f'2026-10-02T{hour+3:02d}:00:00Z',method='gas',
                    period_start=f'2026-10-02T{hour:02d}:00:00Z',
                    period_end=f'2026-10-02T{hour+6:02d}:00:00Z',lifetime_us=value)
    valid=gas(0,100)
    rejected=dict(gas(6,None),fit_status='unavailable')
    last=gas(12,120)
    assert lifetime.connected_values([valid,rejected,last])[1]==[100,120]
    assert lifetime.connected_values([valid,last])[1]==[100,120]
    adjacent=gas(6,110)
    assert lifetime.connected_values([valid,adjacent])[1]==[100,110]



def test_published_json_only_corroborates_and_cannot_seed_plot(tmp_path):
    reference=tmp_path/'reference.json'
    reference.write_text(json.dumps({'lifetimes':[dict(timestamp='2026-10-02T12:00:00Z',lifetime_us=1400,error_us=50)]}))
    output=tmp_path/'history.json'
    empty=lifetime.update_json(output,[],reference_histories=[reference])
    assert empty['plot_lifetimes']==[] and empty['lifetimes']==[]
    assert empty['reference_comparisons'][0]['match_status']=='unmatched'
    derived=dict(timestamp='2026-10-02T07:00:00-05:00',lifetime_us=1200,error_us=40,
        input_file='source.FLOW.hdf5',sample='all_mip',method='track',calculation_source='flow_segments',fit_status='ok')
    data=lifetime.update_json(output,[derived],reference_histories=[reference])
    assert len(data['reference_lifetimes'])==1
    assert data['plot_lifetimes'][0]['lifetime_us']==1200
    comparison=data['reference_comparisons'][0]
    assert comparison['match_status']=='matched' and comparison['difference_us']==-200
    assert data['lifetimes'][0]==derived


def test_refit_reads_segment_values_not_hdf5_lifetime_metadata(tmp_path,monkeypatch):
    import lifetime_funcs as funcs
    import matplotlib.pyplot as plt
    path=tmp_path/'source.FLOW.hdf5'
    a=np.array([(2.,12.,100.,4.),(3.,30.,200.,9.)],dtype=[('dx','f8'),('dQ','f8'),('t','f8'),('nhits','f8')])
    with h5py.File(path,'w') as f:
        f.create_dataset(lifetime.LEGACY_SAMPLE[1],data=a)
        f['analysis/rock_muon_segments'].attrs['electron_lifetime_us']=999999
    def fit(**kwargs):
        np.testing.assert_array_equal(kwargs['dqdx'],[6,10])
        np.testing.assert_array_equal(kwargs['nhits'],[2,3])
        return 1200,40,plt.figure()
    monkeypatch.setattr(funcs,'langau_lifetime',fit)
    result=lifetime.fit_samples(path,tmp_path/'fit.png',timestamp='2026-10-02T12:00:00Z')[0]
    assert result['lifetime_us']==1200 and lifetime.is_flow_track(result)
    assert result['flow_datasets'][0]['shape']==[2]


def test_remote_flow_uses_conditional_ranges_without_full_download(tmp_path):
    import hashlib
    import threading
    from http.server import BaseHTTPRequestHandler,ThreadingHTTPServer
    from flow_input import open_flow
    path=tmp_path/'remote.h5'
    data=np.arange(16,dtype='f8')
    with h5py.File(path,'w') as f:
        f.create_dataset('unrelated_bulk_payload',data=np.zeros(8*1024*1024,dtype='u1'))
        f.create_dataset(lifetime.LEGACY_SAMPLE[1],data=data)
    payload=path.read_bytes();etag='"fixture-etag"';transferred=[]
    class Handler(BaseHTTPRequestHandler):
        def log_message(self,*args): pass
        def do_HEAD(self):
            self.send_response(200);self.send_header('Content-Length',str(len(payload)))
            self.send_header('ETag',etag);self.send_header('Accept-Ranges','bytes');self.end_headers()
        def do_GET(self):
            if self.headers.get('If-Match')!=etag or not self.headers.get('Range'):
                self.send_error(412);return
            start,end=map(int,self.headers['Range'].removeprefix('bytes=').split('-'))
            end=min(end,len(payload)-1);content=payload[start:end+1]
            self.send_response(206);self.send_header('Content-Length',str(len(content)))
            self.send_header('Content-Range',f'bytes {start}-{end}/{len(payload)}')
            self.end_headers();self.wfile.write(content);transferred.append(len(content))
    server=ThreadingHTTPServer(('127.0.0.1',0),Handler)
    thread=threading.Thread(target=server.serve_forever,daemon=True);thread.start()
    try:
        with open_flow(f'http://127.0.0.1:{server.server_port}/remote.h5') as (source,provenance):
            np.testing.assert_array_equal(source[lifetime.LEGACY_SAMPLE[1]][:],data)
        assert 0 < sum(transferred) < len(payload)/2
        assert provenance['source_etag']==etag
        assert sum(transferred) <= provenance['requested_range_bytes'] < len(payload)/2
    finally:
        server.shutdown();server.server_close();thread.join()



def test_annotations_survive_refresh_and_appear_in_interactive_plot(tmp_path):
    path=tmp_path/'history.json';output=tmp_path/'overlay.png'
    event=dict(timestamp='2026-09-17T10:32:00-05:00',label='DF-560 range: 0–10 to 0–1 ppm')
    lifetime.update_json(path,[],annotations=[event])
    lifetime.update_json(path,[],output)
    assert json.loads(path.read_text())['annotations']==[event]
    html=Path(str(output)+'.html').read_text()
    assert 'DF-560 range' in html
    assert Image.open(output).size==(3000,2000)


def test_raw_prm_uses_observation_times_and_survives_refresh(tmp_path):
    """A retry must not double-count readings or turn daily means into raw data."""
    from purity_sources import prm_observations
    source = snapshot(tmp_path, [row('prm_lifetime', .001, '2026-09-29T01:00:00-05:00'),
                                 row('prm_lifetime', .003, '2026-09-29T02:00:00-05:00')])
    frame = load_measurements([source, source])
    raw = prm_observations(frame)
    means = aggregate_with_quality([], measurements=frame)
    path = tmp_path/'history.json'
    lifetime.update_json(path, means, raw_prm_measurements=raw)
    lifetime.update_json(path, means, raw_prm_measurements=raw)
    saved = lifetime.update_json(path, [])
    assert len(saved['raw_prm_measurements']) == 2
    series = lifetime.overlay_groups(saved['lifetimes'], True, saved['raw_prm_measurements'])
    group = series[('prm', 'prm_lifetime')]
    assert [r['lifetime_us'] for r in group['points']] == [1000, 3000]
    assert all(r['error_us'] is None for r in group['points'])
    assert [r['lifetime_us'] for r in group['means']] == [2000]
    assert aware_time(group['means'][0]['timestamp']) == aware_time('2026-09-29T01:30:00-05:00')
    assert aware_time(group['means'][0]['first_observed_at']) == aware_time('2026-09-29T01:00:00-05:00')
    assert aware_time(group['means'][0]['last_observed_at']) == aware_time('2026-09-29T02:00:00-05:00')
    assert all(r['timestamp'] != group['means'][0]['timestamp'] for r in group['points'])
    old = lifetime.overlay_groups(means, True, [])
    assert old[('prm', 'prm_lifetime')]['points'] == []


def test_review_is_persistent_reversible_and_does_not_erase_low_fits(tmp_path):
    """Review by full source identity, never by the measured lifetime value."""
    def fit(name, value):
        return dict(timestamp='2026-10-05T17:00:00-05:00', input_file=name,
                    calculation_source='flow_segments', method='track', sample='all_mip',
                    lifetime_us=value, error_us=20, fit_status='ok')
    rows = [fit('noise/file.FLOW.hdf5', 700), fit('cosmic/file.FLOW.hdf5', 900)]
    policy = {'rules': [dict(input_file=rows[0]['input_file'], action='exclude', reason='Intentional noise study')]}
    path = tmp_path/'history.json'
    lifetime.update_json(path, rows, track_review=policy)
    data = lifetime.update_json(path, rows)  # A refit must retain the review.
    assert data['plot_lifetimes'][0]['lifetime_us'] == 900
    assert data['plot_lifetimes'][0]['n_excluded'] == 1
    assert [r['lifetime_us'] for r in data['lifetimes']] == [700, 900]
    assert all(r['fit_status'] == 'ok' for r in data['lifetimes'])
    series = lifetime.overlay_groups(data['lifetimes'], True)
    assert len(series[('all_mip', '')]['points']) == 2
    assert lifetime.connected_values(series[('all_mip', '')]['means'])[1] == [900]
    restored = lifetime.update_json(path, [], track_review={'rules': []})
    assert restored['plot_lifetimes'][0]['lifetime_us'] == 800
    assert not any(r.get('monitoring_excluded') for r in restored['lifetimes'])
    assert not any('monitoring_excluded' in r for r in rows)  # No caller mutation.


def test_raw_plot_only_connects_means_and_publishes_exclusions(tmp_path, monkeypatch):
    """Inspect actual Plotly traces, not just an HTML file's existence."""
    import plotly.graph_objects as go
    captured = []
    original = go.Figure.write_html
    def capture(self, *args, **kwargs):
        captured.append(self.to_plotly_json())
        return original(self, *args, **kwargs)
    monkeypatch.setattr(go.Figure, 'write_html', capture)
    rows = [dict(timestamp=f'2026-10-05T{hour}:00:00-05:00', input_file=f'file{i}.h5',
                 sample='all_mip', method='track', calculation_source='flow_segments',
                 lifetime_us=value, error_us=30, fit_status='ok')
            for i, (hour, value) in enumerate([('12', 1000), ('13', 1400), ('14', 600)])]
    policy = {'rules': [dict(input_file='file2.h5', action='exclude', reason='Noise study') ]}
    output = tmp_path/'overlay.png'
    lifetime.update_json(tmp_path/'history.json', rows, output, track_review=policy)
    for p in (output, lifetime.raw_plot_path(output)):
        assert Image.open(p).size == (3000, 2000)
        assert Path(str(p)+'.html').exists()
    traces = captured[1]['data']
    assert [t['y'] for t in traces if t['mode'] == 'lines'] == [[1200]]
    points = [t for t in traces if t['mode'] == 'markers']
    assert points[0]['y'] == [1000, 1400]
    excluded = next(t for t in points if 'excluded from mean' in t['name'])
    assert excluded['y'] == [600] and 'Noise study' in excluded['text'][0]
    mean = next(t for t in points if t['name'].endswith('(mean)'))
    assert mean['y'] == [1200]
    assert aware_time(mean['x'][0]) == aware_time('2026-10-05T12:30:00-05:00')
    assert 'Raw PRM readings unavailable' in Path(str(lifetime.raw_plot_path(output))+'.html').read_text()


def test_raw_prm_conflicts_and_ambiguous_reviews_fail_before_publication(tmp_path):
    path = tmp_path/'history.json'
    raw = dict(timestamp='2026-10-01T00:00:00Z', source='PRM', lifetime_us=1000)
    lifetime.update_json(path, [], raw_prm_measurements=[raw])
    before = path.read_bytes()
    with pytest.raises(ValueError, match='conflicting raw PRM'):
        lifetime.update_json(path, [], raw_prm_measurements=[dict(raw, lifetime_us=900)])
    assert path.read_bytes() == before
    entry = dict(input_file='a.h5', calculation_source='flow_segments', method='track')
    rule = dict(input_file='a.h5', action='exclude', reason='study')
    with pytest.raises(ValueError, match='overlapping'):
        lifetime.apply_track_review([entry], {'rules': [rule, rule]})
    with pytest.raises(ValueError, match='reason'):
        lifetime.apply_track_review([entry], {'rules': [dict(rule, reason='')]})


def test_track_audit_shows_singletons_and_leave_one_out_influence():
    from track_audit import audit_tracks
    def entry(name, value, stamp):
        return dict(input_file=name, lifetime_us=value, error_us=20, timestamp=stamp,
                    sample='all_mip', method='track', calculation_source='flow_segments')
    rows = [entry('singleton', 860, '2026-10-03T13:00:00-05:00'),
            entry('a', 740, '2026-10-05T17:00:00-05:00'),
            dict(entry('b', 920, '2026-10-05T17:01:00-05:00'), monitoring_excluded=True)]
    singleton, pair = audit_tracks(rows)['windows']
    assert singleton['contributions'][0]['leave_one_out_mean_us'] is None
    assert pair['unreviewed_mean_us'] == 830
    assert pair['monitoring_mean_us'] == 740
    assert [r['leave_one_out_mean_us'] for r in pair['contributions']] == [920, 740]
    assert [r['mean_change_when_omitted_us'] for r in pair['contributions']] == [90, -90]
