"""Deployment contracts with synthetic inputs; these are not physics evidence.

Real detector replay is documented separately. These tests protect scheduling,
source membership and publication so an operational retry cannot silently turn
one acquisition into two measurements, or replace a good page with a partial one.
"""
import hashlib
import json
import os
import re
from pathlib import Path
import shlex
import subprocess
import sys
import time
from types import SimpleNamespace

import h5py
import numpy as np
from PIL import Image
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO/'actions/lifetime'))
import monitor
import monitor_config
import monitor_plot
import track_shards
import watcher


@pytest.fixture(autouse=True)
def local_locks(monkeypatch):
    monkeypatch.setenv('ARCUBE_NEARLINE_LOCAL_OUTPUT', '1')


def profile(tmp_path, **overrides):
    config = {k:str(tmp_path/k) for k in ('flow_root','packet_root','snapshot_root','state_root','publish_root')}
    config.update(overrides)
    for key in ('flow_root','packet_root','snapshot_root'):
        Path(config[key]).mkdir(parents=True, exist_ok=True)
    path = tmp_path/'monitor config.json'
    path.write_text(json.dumps(config))
    return path, monitor_config.load_config(path)


def acquisition(config, run='stable', stamp='2026_10_05_08_00_00_CDT'):
    flow = Path(config['flow_root'])/run/f'packet-test-{stamp}.FLOW.hdf5'
    native = Path(config['packet_root'])/run/f'packet-test-{stamp}.h5'
    flow.parent.mkdir(parents=True, exist_ok=True)
    native.parent.mkdir(parents=True, exist_ok=True)
    return flow, native


def synthetic_flow(path):
    segments = np.array([(0,1.,10.,3.,50.), (1,1.,20.,3.,100.)],
        dtype=[('rock_segment_id','i4'), ('dx','f8'), ('dQ','f8'), ('nhits','f8'), ('t','f8')])
    with h5py.File(path, 'w') as h:
        h.create_dataset('analysis/rock_muon_segments/data', data=segments)
        h.create_dataset('analysis/rock_muon_tracks/data', data=np.array([(0,),(1,)],dtype=[('event_id','i8')]))
        ref = h.create_dataset('analysis/rock_muon_tracks/ref/analysis/rock_muon_segments/ref', data=[[1,0],[0,1]])
        ref.attrs.update(dset0='/analysis/rock_muon_tracks/data', dset1='/analysis/rock_muon_segments/data')
        h.create_dataset('charge/events/data',data=np.array([(0,1),(1,0)],dtype=[('id','i8'),('n_ext_trigs','i8')]))
        for group in ('geometry_info','lar_info','charge/calib_prompt_hits'):
            h.create_group(group)
        h.create_group('run_info').attrs.update(crs_ticks=.1,e_field=.05,is_mc=False,charge_thresholds='medm',data_packet_type=0)
    return segments


def snapshots(config):
    path = Path(config['snapshot_root'])/'slowcontrols-test.json'
    rows = [dict(timestamp=f'2026-10-05T{hour:02}:00:00-05:00',value=value,quantity='prm_lifetime',unit='us',source='test PRM')
            for hour,value in [(1,1000.),(5,3000.),(9,2000.)]]
    path.write_text(json.dumps(dict(measurements=rows)))
    return path


def make_watcher(config_path, **kwargs):
    options = dict(prog=REPO/'actions/lifetime.sh', paths=[],exts=[],cond=lambda p:True,
        sleep_between_scans=30,min_file_age=0,max_file_age=None,priority=0,dry_run=True,lifetime_config=config_path)
    options.update(kwargs)
    return watcher.Watcher(**options)


def test_dependency_revisions_packet_first_late_and_replaced(tmp_path):
    path, config = profile(tmp_path)
    flow, packet = acquisition(config)
    packet.write_bytes(b'first')
    assert monitor_config.discovery(packet,config) is None
    flow.write_bytes(b'FLOW')
    both = monitor_config.discovery(flow,config)
    assert both == monitor_config.discovery(packet,config)
    packet.unlink()
    only = monitor_config.discovery(flow,config)
    assert only['key'] == both['key'] and only['revision'] != both['revision']
    packet.write_bytes(b'new generation')
    assert monitor_config.discovery(flow,config)['revision'] != both['revision']
    other, _ = acquisition(config,run='different')
    other.write_bytes(b'FLOW')
    assert monitor_config.discovery(other,config)['key'] != both['key']
    assert monitor_config.discovery(flow.with_suffix('.partial'),config) is None


def test_dry_discovery_needs_no_mongo_and_deduplicates_pair(tmp_path, capsys):
    path,config = profile(tmp_path)
    flow,packet = acquisition(config,run='name with spaces')
    flow.write_bytes(b'FLOW'); packet.write_bytes(b'packet')
    snapshots(config)
    w = make_watcher(path)
    w.snarf(); w.snarf()
    commands = [line.removeprefix('Would queue: ') for line in capsys.readouterr().out.splitlines()]
    assert len(commands)==2 and w.db is None
    command = next(c for c in commands if '.FLOW.hdf5' in c)
    assert shlex.split(command)==[str(REPO/'actions/lifetime.sh'),str(flow),'--monitor-config',str(path)]
    # CLI itself must also run in the standard-library-only system Python.
    done = subprocess.run(['/usr/bin/python3',str(REPO/'watcher.py'),str(REPO/'actions/lifetime.sh'),
        '--lifetime-config',str(path),'--once','--dry-run','--min-file-age','0'],capture_output=True,text=True)
    assert done.returncode==0, done.stderr
    assert done.stdout.count('Would queue:')==2


def test_stability_and_freshness_apply_to_companions(tmp_path):
    path,config=profile(tmp_path)
    flow,packet=acquisition(config);flow.write_bytes(b'x');packet.write_bytes(b'y')
    old=time.time()-10000; fresh=time.time()-100
    os.utime(flow,(old,old))
    w=make_watcher(path,min_file_age=30,max_file_age=300)
    assert w.candidate(flow) is None  # Packet is still being written.
    os.utime(packet,(fresh,fresh))
    assert w.candidate(flow) is not None  # Old FLOW, newly completed packet.
    os.utime(packet,(old,old))
    assert w.candidate(flow) is None
    packet.unlink();flow.unlink()
    assert w.candidate(flow) is None


def test_scheduling_spec_deduplicates_across_watcher_restart(tmp_path, monkeypatch):
    path,config=profile(tmp_path);flow,packet=acquisition(config);flow.write_bytes(b'x')
    records=[]
    class Collection:
        def find_one(self,query):
            return next((r for r in records if all(r.get(k)==v for k,v in query.items())),None)
    class Launchpad:
        name='test'
        connection={'test':{'fireworks':Collection()}}
        def add_wf(self,fw):
            fw.fw_id=len(records)+1
            records.append({'spec.'+k:v for k,v in fw.spec.items()})
    monkeypatch.setitem(sys.modules,'fireworks',SimpleNamespace(
        ScriptTask=SimpleNamespace(from_str=lambda script:script),
        Firework=lambda task,name,spec:SimpleNamespace(task=task,name=name,spec=spec)))
    w=make_watcher(path,dry_run=False,launchpad=Launchpad())
    first=w.maybe_make_firework(flow)
    assert first.spec['_purity_key'] and first.spec['_purity_revision']
    assert make_watcher(path,dry_run=False,launchpad=Launchpad()).maybe_make_firework(flow) is None
    packet.write_bytes(b'late')
    assert w.maybe_make_firework(packet) is not None
    assert len(records)==2


@pytest.mark.parametrize('change', ['public_in_state','state_in_public','public_in_input','unexpanded'])
def test_config_rejects_unsafe_or_unresolved_layout(tmp_path,change):
    overrides={
        'public_in_state':dict(publish_root=str(tmp_path/'state_root'/'plots')),
        'state_in_public':dict(state_root=str(tmp_path/'publish_root'/'private')),
        'public_in_input':dict(publish_root=str(tmp_path/'flow_root'/'plots')),
        'unexpanded':dict(calibration_root='${ABSENT_PURITY_TEST_VARIABLE}')}
    with pytest.raises(ValueError):profile(tmp_path,**overrides[change])


def test_track_partition_uses_reference_rows_not_position(tmp_path):
    path=tmp_path/'flow.h5';original=synthetic_flow(path)
    with h5py.File(path,'r') as h:arrays,origin,links=track_shards.stored_segments(h)
    np.testing.assert_array_equal(arrays['cosmic'],original[:1])
    np.testing.assert_array_equal(arrays['beam'],original[1:])
    assert links['source_rows']==dict(beam=[1],cosmic=[0])
    with h5py.File(path,'r+') as h:
        h['analysis/rock_muon_tracks/ref/analysis/rock_muon_segments/ref'][:]=[[1,0],[0,0]]
    with h5py.File(path,'r') as h,pytest.raises(ValueError,match='exactly one'):
        track_shards.stored_segments(h)


@pytest.mark.parametrize('changed_setting', ['calibration', 'selection'])
def test_pool_retries_do_not_duplicate_and_config_changes_do_not_merge(tmp_path,monkeypatch,changed_setting):
    path,config=profile(tmp_path);source,_=acquisition(config);synthetic_flow(source)
    shards=tmp_path/'shards';calls=[]
    def fit(file,plot,timestamp):
        with h5py.File(file,'r') as h:
            calls.append(sum(len(h[p]) for _,p in track_shards.SPLIT_SAMPLES))
        return [dict(timestamp=timestamp,sample='all_mip',fit_status='failed',lifetime_us=None,error_us=None)]
    monkeypatch.setattr(track_shards,'fit_samples',fit)
    meta=dict(timestamp='2026-10-05T08:00:00-05:00',cohort='stable')
    for _ in range(2):track_shards.write_shard(source,shards,'one',meta)
    track_shards.write_shard(source,shards,'two',dict(meta,timestamp='2026-10-05T10:00:00-05:00'))
    rows=track_shards.pool_tracks(shards,tmp_path/'fits')
    assert len(rows)==1 and rows[0]['n_files']==2 and calls==[4]
    assert rows[0]['period_start']=='2026-10-05T06:00:00-05:00'
    assert track_shards.pool_tracks(shards,tmp_path/'fits')==rows and calls==[4]
    with h5py.File(source,'r+') as h:
        if changed_setting == 'calibration':
            h['charge/calib_prompt_hits'].attrs['pedestal_file']='new.json'
        else:
            h['analysis/rock_muon_segments'].attrs['length_cut']=100.
    track_shards.write_shard(source,shards,'two',meta)
    revised=track_shards.pool_tracks(shards,tmp_path/'fits')
    assert len(revised)==2 and sum(r['n_files'] for r in revised)==2
    # Altering an array must fail even if a matching fit result is cached.
    descriptor=json.loads((shards/'one.sample.json').read_text())
    object_path=shards/descriptor['arrays_file']
    with np.load(object_path) as f:arrays={k:f[k] for k in f.files}
    arrays['beam']['dQ']*=2
    np.savez_compressed(object_path,**arrays)
    with pytest.raises(ValueError,match='hash mismatch'):track_shards.pool_tracks(shards,tmp_path/'fits')


def test_split_selection_settings_and_timing_are_checked(tmp_path):
    source=tmp_path/'flow.h5'; original=synthetic_flow(source)
    with h5py.File(source,'r+') as h:
        for (sample,path),timing in zip(track_shards.SPLIT_SAMPLES,['beam','off_beam']):
            h.create_dataset(path,data=original[:1] if sample=='beam' else original[1:])
            h[str(Path(path).parent)].attrs.update(event_timing=timing,length_cut=100.,
                selected_track_count=1 if sample=='beam' else 2)
        arrays,_,_=track_shards.stored_segments(h)
        assert sum(map(len,arrays.values()))==len(original)
        # Counts and old fit results must not split otherwise equal settings.
        before=track_shards.selection_metadata(h)
        group=h['analysis/cosmic_muon_segments']
        group.attrs['selected_track_count']=99
        group.attrs['electron_lifetime_us']=2000.
        assert track_shards.selection_metadata(h)==before
        group.attrs['length_cut']=150.
        with pytest.raises(ValueError,match='different recorded settings'):
            track_shards.stored_segments(h)
        group.attrs['length_cut']=100.
        group.attrs['event_timing']='all'
        with pytest.raises(ValueError,match='unexpected event_timing'):
            track_shards.stored_segments(h)


def test_prm_only_end_to_end_atomic_output_and_reading_time(tmp_path,monkeypatch):
    path,config=profile(tmp_path);snap=snapshots(config)
    payload=monitor.run_monitor(path,snap)
    assert len(payload['lifetimes'])==1
    row=payload['lifetimes'][0]
    assert row['lifetime_us']==2000 and row['error_us']==pytest.approx(1000/np.sqrt(3))
    root=Path(config['publish_root'])
    assert Image.open(root/'elifetime_time_series.png').size==(3000,2000)
    html=(root/'elifetime_time_series.png.html').read_text()
    assert '2026-10-05T05:00:00-05:00' in html  # Actual reading-time centroid.
    assert not re.search(r"<script[^>]+\bsrc\s*=", html)
    current=os.readlink(root/'.lifetime-current')
    previous=(root/'elifetime_time_series.png').read_bytes()
    def fail(*args):raise RuntimeError('renderer failed')
    monkeypatch.setattr(monitor_plot,'render',fail)
    with pytest.raises(RuntimeError,match='renderer failed'):monitor.publish(payload,config)
    assert os.readlink(root/'.lifetime-current')==current
    assert (root/'elifetime_time_series.png').read_bytes()==previous


def test_packet_failure_publishes_controls_and_tracks_but_exits_unsuccessfully(tmp_path,monkeypatch):
    path,config=profile(tmp_path);flow,packet=acquisition(config);synthetic_flow(flow);packet.write_bytes(b'invalid')
    snapshots(config)
    with pytest.raises(RuntimeError,match='pedestal calibration'):
        monitor.run_monitor(path,flow)
    root=Path(config['publish_root'])
    payload=json.loads((root/'elifetime_status.json').read_text())
    assert any(r['sample']=='prm' for r in payload['lifetimes'])
    assert any(r['sample']=='all_mip' for r in payload['lifetimes'])
    assert payload['inputs'][0]['packet_status']=='failed'
    assert payload['input_errors']


def test_failed_packet_marker_is_open_and_external_subset_is_not_drawn(tmp_path,monkeypatch):
    _,config=profile(tmp_path)
    row=dict(timestamp='2026-10-05T08:00:00-05:00',period_start='2026-10-05T06:00:00-05:00',
             period_end='2026-10-05T12:00:00-05:00',source='stable',fit_status='ok',lifetime_us=1000,error_us=100)
    rows=[dict(row,sample='all_mip',method='track_pool'),dict(row,sample='beam',method='track_pool'),
          dict(row,sample='packet',method='packet_pool',fit_status='unavailable',lifetime_us=None,error_us=None,
               candidate_lifetime_us=3300,fit_message='attenuation model disagreement')]
    captured=[]
    original=monitor_plot.go.Figure.to_html
    def record(self,*args,**kwargs):
        captured.extend(self.to_plotly_json()['data'])
        return original(self,*args,**kwargs)
    monkeypatch.setattr(monitor_plot.go.Figure,'to_html',record)
    monitor_plot.render(dict(lifetimes=rows,raw_prm_measurements=[],generated_at=row['timestamp'],
        inputs=[],input_errors=[],uncertainty='test'),tmp_path,config)
    assert {t['name'] for t in captured}=={'All tracks · 6 h','Packets · 6 h'}
    candidate=next(t for t in captured if t['name'].startswith('Packets'))
    assert candidate['marker']['symbol']==['cross-open'] and candidate['error_y']['array']==[None]
    assert rows[-1]['lifetime_us'] is None  # Rendering must not promote failures.
    assert '<td>beam</td>' in (tmp_path/'elifetime_diagnostics.html').read_text()


def test_shell_dispatch_quotes_input_and_config_paths(tmp_path):
    recorder=tmp_path/'python recorder'
    recorder.write_text('#!/usr/bin/env python3\nimport json,sys\nprint(json.dumps(sys.argv[1:]))\n')
    recorder.chmod(0o755)
    env=dict(os.environ,ARCUBE_NEARLINE_PYTHON=str(recorder))
    cmd=['bash',str(REPO/'actions/lifetime.sh'),'input with spaces.FLOW.hdf5','--monitor-config','config with spaces.json']
    result=subprocess.run(cmd,env=env,capture_output=True,text=True,check=True)
    assert json.loads(result.stdout)[1:]==['--monitor-config','config with spaces.json','--input-file','input with spaces.FLOW.hdf5']


def test_configured_cli_rejects_ignored_legacy_options(tmp_path):
    path,_=profile(tmp_path)
    result=subprocess.run([sys.executable,str(REPO/'actions/lifetime/lifetime.py'),
        '--monitor-config',str(path),'--export-slow-controls',str(tmp_path/'export.json')],
        capture_output=True,text=True)
    assert result.returncode==2 and 'supply only an optional --input-file' in result.stderr
    assert not (tmp_path/'export.json').exists()
