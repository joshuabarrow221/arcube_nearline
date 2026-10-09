"""Synthetic contract tests only; none of these objects is detector evidence."""
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'actions/lifetime'))
import packet_pooling as pooling
import lifetime
from purity_sources import aware_time

CONFIG = dict(time_bin_us=10.,drift_max_us=200.,minimum_bin_objects=50,
              maximum_slice_chi2_ndf=3.,fit_tmin_us=20.,fit_tmax_us=180.)


def shard(directory,name,stamp,cohort='stable',coverage=.99,configuration=None):
    frame=pd.DataFrame(dict(q_sum=[10.,11.],time_us=[25.,35.],ext=[0,1],io_group=[1,2]))
    return pooling.store_sample(frame,dict(packet_file=name,timestamp=stamp,cohort=cohort,
        configuration=configuration or CONFIG,pedestal_coverage=coverage,duration_s=60.),directory)


def test_windows_fixed_before_fits_and_preserve_dst():
    assert pooling.window_start('2026-10-05T17:59:00-05:00',6,'America/Chicago').hour==12
    assert pooling.window_start('2026-10-05T17:59:00-05:00',24,'America/Chicago').hour==0
    assert pooling.window_start('2026-10-05T17:59:00-05:00',48,'America/Chicago') == pooling.window_start('2026-10-04T01:00:00-05:00',48,'America/Chicago')
    with pytest.raises(ValueError):pooling.window_start('2026-10-05T00:00:00Z',13,'UTC')
    start=pooling.window_start('2026-03-08T12:00:00-05:00',24,'America/Chicago')
    end=(start.tz_localize(None)+pd.Timedelta(hours=24)).tz_localize('America/Chicago')
    assert (end-start).total_seconds()==23*3600


def test_pooling_deduplicates_copies_and_separates_cohorts_and_configuration(tmp_path,monkeypatch):
    monkeypatch.setattr(pooling,'fit_objects',lambda objects,*args:dict(lifetime_us=1000.,error_us=10.,fit_status='ok'))
    shard(tmp_path,'/first/a.h5','2026-10-05T01:00:00Z')
    shard(tmp_path,'/copy/a.h5','2026-10-05T01:00:00Z')
    shard(tmp_path,'b.h5','2026-10-05T02:00:00Z')
    shard(tmp_path,'c.h5','2026-10-05T02:00:00Z',cohort='noise-study')
    shard(tmp_path,'d.h5','2026-10-05T02:00:00Z',configuration=dict(CONFIG,time_bin_us=20.))
    rows=pooling.pool_samples(tmp_path,24)
    assert len(rows)==3
    assert sorted(r['n_files'] for r in rows)==[1,1,2]
    assert sum(r['n_objects'] for r in rows)==8
    cached=pooling.pool_samples(tmp_path,24)
    # Identical membership reuses the fit, but each inventory has a fresh
    # timestamp to order publication against concurrent extraction/retries.
    assert [{k:v for k,v in r.items() if k!='pool_snapshot_at'} for r in cached] == [
        {k:v for k,v in r.items() if k!='pool_snapshot_at'} for r in rows]


def test_low_coverage_cannot_be_hidden_by_a_good_fit_or_another_file(tmp_path,monkeypatch):
    monkeypatch.setattr(pooling,'fit_objects',lambda *args:dict(lifetime_us=1000.,error_us=10.,fit_status='ok'))
    shard(tmp_path,'a.h5','2026-10-05T01:00:00Z',coverage=.8)
    shard(tmp_path,'b.h5','2026-10-05T02:00:00Z',coverage=1.)
    result=pooling.pool_samples(tmp_path,24)[0]
    assert result['lifetime_us'] is None and result['fit_status']=='unavailable'
    assert 'pedestal coverage' in result['fit_message']


@pytest.mark.parametrize('cached', [False, True])
def test_changed_payload_rejected_before_pooling(tmp_path, monkeypatch, cached):
    monkeypatch.setattr(pooling,'fit_objects',lambda *args:dict(lifetime_us=1000.,error_us=10.,fit_status='ok'))
    meta=shard(tmp_path,'a.h5','2026-10-05T01:00:00Z')
    if cached:
        assert pooling.pool_samples(tmp_path,24)[0]['fit_status']=='ok'
    p=tmp_path/meta['objects_file']
    with np.load(p,allow_pickle=False) as f:a=f['objects']
    a['q_sum']*=2
    np.savez_compressed(p,objects=a)
    with pytest.raises(ValueError,match='hash mismatch'):pooling.pool_samples(tmp_path,24)


def test_attenuation_gates_and_covariance_on_known_curve(monkeypatch):
    import packet_lifetime as reference
    monkeypatch.setattr(reference,'fit_langau',lambda q:dict(mpv=float(np.mean(q)),red_chi2=1.))
    monkeypatch.setattr(reference,'bootstrap_mpv',lambda q,n,seed:(float(np.mean(q)),.05,n))
    # Exact slice centers lie on tau=2 ms. This tests the second fit's units
    # and validity gates independently of the reference slice fitter.
    a=np.zeros(20*60,dtype=[('q_sum','f8'),('time_us','f8')])
    a['time_us']=np.repeat(np.arange(5,200,10),60)
    a['q_sum']=20*np.exp(-a['time_us']/2000)
    good=pooling.fit_objects(a,CONFIG)
    assert good['fit_status']=='ok' and good['lifetime_us']==pytest.approx(2000,rel=1e-5)
    assert good['error_us']>0
    short=pooling.fit_objects(a[a['time_us']<45],CONFIG)
    assert short['lifetime_us'] is None and 'too few' in short['fit_message']
    a['q_sum']+=np.repeat([2.,-2.]*10,60)
    bad=pooling.fit_objects(a,CONFIG)
    assert bad['lifetime_us'] is None and 'model disagreement' in bad['fit_message']
    monkeypatch.setattr(reference,'fit_langau',lambda q:dict(mpv=float('nan'),red_chi2=float('nan')))
    unavailable=pooling.fit_objects(a,CONFIG)
    assert unavailable['lifetime_us'] is None
    json.dumps(unavailable,allow_nan=False)  # Failure diagnostics remain valid JSON.


def test_pool_history_does_not_regress_and_only_selected_period_is_published(tmp_path):
    def r(hours,evaluated,value):
        return dict(method='packet_pool',sample='packet',source=f'cohort/{hours}',window_hours=hours,
                    timestamp='2026-10-05T12:00:00Z',period_start='2026-10-05T00:00:00Z',
                    evaluated_at=evaluated,lifetime_us=value,error_us=10,fit_status='ok')
    p=tmp_path/'history.json'
    lifetime.update_json(p,[r(24,'2026-10-06T00:00:00Z',1000),r(6,'2026-10-06T00:00:00Z',1100)])
    data=lifetime.update_json(p,[r(24,'2026-10-05T00:00:00Z',500)],packet_window_hours=24)
    assert [v['lifetime_us'] for v in data['plot_lifetimes']]==[1000]
    assert len(data['lifetimes'])==2


def test_inventory_retires_changed_configuration_and_ignores_late_worker(tmp_path, monkeypatch):
    """A re-extracted file must not appear twice under old and new settings."""
    monkeypatch.setattr(pooling,'fit_objects',lambda *args:dict(lifetime_us=1000.,error_us=10.,fit_status='ok'))
    samples=tmp_path/'samples'; history=tmp_path/'history.json'
    shard(samples,'a.h5','2026-10-05T01:00:00Z')
    old_inventory={}; old=pooling.pool_samples(samples,snapshot=old_inventory)
    lifetime.update_json(history,old,packet_snapshot=old_inventory)
    shard(samples,'a.h5','2026-10-05T01:00:00Z',configuration=dict(CONFIG,pedestal_identity='new calibration'))
    new_inventory={}; new=pooling.pool_samples(samples,snapshot=new_inventory)
    lifetime.update_json(history,new,packet_snapshot=new_inventory)
    data=lifetime.update_json(history,old,packet_snapshot=old_inventory)
    assert len(data['lifetimes'])==2  # Preserve audit history.
    assert len(data['plot_lifetimes'])==1
    assert data['plot_lifetimes'][0]['source']==new[0]['source']


def test_static_calibration_identity_and_missing_channels(tmp_path):
    """Synthetic geometry/calibration reproduces Q_raw and rejects wrong gain.

    Unknown channels must remain missing instead of getting a default pedestal.
    A changed Q_raw relation fails even when filenames and metadata look right.
    """
    import h5py
    from calibration_pedestal import calibration_map, tile_lookup
    from packet_pedestal import key
    pedfile=tmp_path/'pedestal.json'; pedfile.write_text(json.dumps({'1000300405':{'pedestal_mv':580.}}))
    flowfile=tmp_path/'test.FLOW.hdf5'
    packet_dtype=[(c,'i4') for c in ('io_group','io_channel','chip_id','channel_id','packet_type','dataword')]
    hit_dtype=[(c,'i4') for c in ('io_group','io_channel','chip_id','channel_id')]+[('Q_raw','f8')]
    packets=np.zeros(150,dtype=packet_dtype)
    for c,v in dict(io_group=1,io_channel=1,chip_id=4,channel_id=5).items():packets[c]=v
    packets['dataword']=np.arange(150)%50+30
    hits=np.zeros(150,dtype=hit_dtype)
    for c in ('io_group','io_channel','chip_id','channel_id'):hits[c]=packets[c]
    hits['Q_raw']=(packets['dataword']/256*(1568-478.1)+478.1-580)/4.522
    with h5py.File(flowfile,'w') as f:
        group=f.create_group('geometry_info/tile_id')
        meta=np.zeros(1,dtype=[('min_max_keys','i8',(2,2))]);meta['min_max_keys']=[[[1,2],[1,2]]]
        group.attrs['meta']=meta
        tiles=np.zeros(5,dtype=[('data','i8'),('filled','?')]);tiles[1]=(3,True)
        group['data']=tiles
        assert tile_lookup(group,[1,2,1,0],[1,1,2,1]).tolist()==[3,-1,-1,-1]
        f['charge/packets/data']=packets; f['charge/calib_prompt_hits/data']=hits
        f['charge/calib_prompt_hits/ref/charge/packets/ref']=np.column_stack([np.arange(150)]*2)
        f['charge/calib_prompt_hits'].attrs['pedestal_file']='/calibration/pedestal.json'
    channels=packets[:2].copy();channels['channel_id'][1]=6
    mapping,info=calibration_map(str(flowfile),pedfile,channels)
    address=int(key(1,1,4,5))
    assert mapping=={address:pytest.approx((580-478.1)*256/(1568-478.1))}
    assert info['maximum_Q_raw_residual']<1e-10
    assert info['verified_pairs']==150
    with h5py.File(flowfile,'a') as f:
        wrong=hits.copy();wrong['Q_raw']*=2
        f['charge/calib_prompt_hits/data'][:]=wrong
    with pytest.raises(ValueError,match='identity not verified'):calibration_map(str(flowfile),pedfile,channels)
    with h5py.File(flowfile,'a') as f:f['charge/calib_prompt_hits'].attrs['gain_file']='custom.json'
    with pytest.raises(ValueError,match='non-default'):calibration_map(str(flowfile),pedfile,channels)


def test_identical_objects_in_different_windows_reuse_fit(tmp_path,monkeypatch):
    calls=[]
    def fit(objects,*args):
        calls.append(len(objects))
        return dict(lifetime_us=1000.,error_us=10.,fit_status='ok')
    monkeypatch.setattr(pooling,'fit_objects',fit)
    shard(tmp_path,'a.h5','2026-10-05T01:00:00Z')
    six=pooling.pool_samples(tmp_path,6)
    day=pooling.pool_samples(tmp_path,24)
    assert calls==[2] and six[0]['window_hours']==6 and day[0]['window_hours']==24


def test_pool_refresh_failure_still_publishes_prm(tmp_path):
    import subprocess
    samples=tmp_path/'samples';samples.mkdir()
    (samples/'broken.sample.json').write_text('{broken')
    history=tmp_path/'history.json'
    lifetime.update_json(history,[dict(method='prm',sample='prm',source='test',fit_status='ok',
        timestamp='2026-10-05T12:00:00Z',period_start='2026-10-05T00:00:00Z',
        period_end='2026-10-06T00:00:00Z',lifetime_us=1000.,error_us=None,n_measurements=1)])
    result=subprocess.run([sys.executable,str(Path(lifetime.__file__)),
        '--output_file_json',str(history),'--packet-samples',str(samples)],capture_output=True,text=True)
    assert result.returncode==0,result.stderr
    data=json.loads(history.read_text())
    assert [r['lifetime_us'] for r in data['plot_lifetimes']]==[1000.]
    assert any(r.get('method')=='packet_pool_refresh' and r['fit_status']=='failed' for r in data['lifetimes'])


def candidate_row(hours=24, source='cosmic', value=1800.):
    """Synthetic display fixture, never an input to the detector study."""
    return dict(method='packet_pool',sample='packet',source=source,cohort='cosmic',window_hours=hours,
                timestamp='2026-09-28T12:00:00-05:00',period_start='2026-09-28T00:00:00-05:00',
                period_end='2026-09-29T00:00:00-05:00',fit_status='unavailable',lifetime_us=None,error_us=None,
                candidate_lifetime_us=value,fit_message='insufficient pedestal coverage',minimum_pedestal_coverage=.92,
                alpha_per_ms=1000/1800.,alpha_error_per_ms=.1,covariance=[[.01,0],[0,.01]],n_files=1,n_objects=30000)


def test_candidate_display_does_not_promote_failed_fits_or_invent_values(tmp_path):
    row=candidate_row()
    omitted=[dict(row,candidate_lifetime_us=v) for v in (None,0,-1,float('nan'),float('inf'))]
    omitted += [dict(row,method='track'),dict(row,monitoring_excluded=True),
                dict(row,fit_status='ok',lifetime_us=1800.,error_us=30.)]
    points=lifetime.packet_candidate_points([row,*omitted])
    assert len(points)==1 and points[0]['display_lifetime_us']==1800.
    assert points[0]['display_error_us']==pytest.approx(324.)
    assert not lifetime.valid_lifetime(points[0])
    assert points[0]['lifetime_us'] is None and points[0]['fit_status']=='unavailable'
    assert 'display_lifetime_us' not in row  # No caller/history mutation.
    no_cov=lifetime.packet_candidate_points([dict(row,covariance=None)])[0]
    assert no_cov['display_error_us'] is None
    path=tmp_path/'history.json'
    data=lifetime.update_json(path,[row,candidate_row(48,'cosmic/48',1980.)],packet_window_hours=24)
    assert [r['display_lifetime_us'] for r in data['packet_candidate_points']]==[1800.]
    assert all(not lifetime.valid_lifetime(r) for r in data['plot_lifetimes'])
    assert len(data['lifetimes'])==2


def test_candidate_png_and_html_show_failure_and_conditional_error(tmp_path,monkeypatch):
    from PIL import Image
    import plotly.graph_objects as go
    from matplotlib.axes import Axes
    captured=[];collections=[]
    write=go.Figure.write_html;scatter=Axes.scatter
    def capture(self,*args,**kwargs):
        captured.append(self.to_plotly_json());return write(self,*args,**kwargs)
    def capture_scatter(self,*args,**kwargs):
        result=scatter(self,*args,**kwargs);collections.append(result);return result
    monkeypatch.setattr(go.Figure,'write_html',capture)
    monkeypatch.setattr(Axes,'scatter',capture_scatter)
    row=candidate_row();output=tmp_path/'plot.png'
    # Include a second window with no estimate. It must not become a zero,
    # guessed point, or accepted lifetime merely because candidates are shown.
    missing=dict(row,timestamp='2026-09-29T12:00:00-05:00',period_start='2026-09-29T00:00:00-05:00',candidate_lifetime_us=None)
    lifetime.update_json(tmp_path/'history.json',[row,missing],output)
    for png,figure in zip([output,lifetime.raw_plot_path(output)],captured):
        assert Image.open(png).size==(3000,2000)
        trace=next(t for t in figure['data'] if 'Packet fit candidate' in t.get('name',''))
        assert trace['y']==[1800.] and trace['marker']['symbol']=='hexagon-open'
        assert trace['line']['dash']=='dash' and trace['error_y']['array']==pytest.approx([324.])
        assert 'insufficient pedestal coverage' in trace['text'][0]
        assert '92.00%' in trace['text'][0]
        assert 'Conditional fit-covariance' in trace['text'][0]
    assert len(collections)==2
    assert all(len(c.get_facecolors())==0 for c in collections)
    assert all(c.get_linestyles()[0][1] is not None for c in collections)
