"""Charge-scale semantics and full-pool diagnostic checks, not detector calibration."""
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import h5py
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'actions/lifetime'))
import charge_scale as scale


def receipt():
    times=np.linspace(45.,175.,14)
    amplitude,tau=16.,320.
    return dict(selection='quality-stored',period_start='2026-11-02T18:00:00+00:00',
        period_end='2026-11-03T00:00:00+00:00',n_valid_segments=123,errors=[],
        members=[dict(input_file='/run/cohort/packet.FLOW.hdf5')],
        reconstruction_signature={'quality':{'quality_code_sha256':'abc'}},
        pooled_segments='/run/pool/segments.h5',fit_status='ok',lifetime_us=tau,error_us=6.,
        fit=dict(fit_status='ok',fit_diagnostics=dict(estimator_version='demario-proxy-active-bins-v3-min10',
            candidate_parameters=dict(amplitude=amplitude,lifetime_us=tau),
            attenuation_covariance=[[.04,-.6],[-.6,36.]],
            fitted_bins=list(range(14)),drift_time_us=times.tolist(),
            mpv=(amplitude*np.exp(-times/tau)).tolist(),mpv_error=[.1]*14,chi2_ndf=1.3,ndf=12)))


def test_saved_intercept_uses_proxy_units_and_full_covariance():
    result=scale.charge_scale_from_summary(receipt())
    assert result['status']=='ok'
    assert result['units']=='legacy proxy units'
    assert result['intercept']==16.
    assert result['intercept_error']==pytest.approx(.2)
    assert result['tau_us']==320.
    assert result['correlation']==pytest.approx(-.5)
    assert result['curve']['x'][0]==0.
    assert result['curve']['y'][0]==16.
    assert result['curve']['error'][0]==pytest.approx(.2)
    assert result['fit_range_us']==[45.,175.]
    index=70
    t=result['curve']['x'][index]
    gradient=np.array([np.exp(-t/320),16*np.exp(-t/320)*t/320**2])
    expected=np.sqrt(gradient@np.array([[.04,-.6],[-.6,36.]])@gradient)
    assert result['curve']['error'][index]==pytest.approx(expected)


@pytest.mark.parametrize('problem',['failed','rejected','input_error','singular','missing_covariance','nonpositive','over_limit'])
def test_invalid_saved_results_do_not_publish_intercept(problem):
    data=receipt()
    diag=data['fit']['fit_diagnostics']
    if problem in ('failed','rejected'):data['fit_status']=problem
    elif problem=='input_error':data['errors']=[{'error':'bad input'}]
    elif problem=='singular':diag['attenuation_covariance']=[[1.,1.],[1.,1.]]
    elif problem=='missing_covariance':diag.pop('attenuation_covariance')
    elif problem=='nonpositive':diag['candidate_parameters']['amplitude']=-1
    else:diag['candidate_parameters']['lifetime_us']=20001.
    result=scale.charge_scale_from_summary(data)
    assert result['status']!='ok'
    assert result['intercept'] is None and result['intercept_error'] is None
    assert result['curve']['x']==[]


def segments(count=10):
    array=np.zeros(19*count,dtype=[('dQ','f8'),('dx','f8'),('nhits','f8'),('t','f8')])
    array['dx']=2.
    array['nhits']=5.
    array['t']=np.repeat((scale.TIME_BINS_US[:-1]+scale.TIME_BINS_US[1:])/2/.1,count)
    array['dQ']=80*np.exp(-array['t']*.1/300.)
    return array


def test_direct_fit_uses_all_segments_and_fixed_bins_not_display_sample(monkeypatch):
    expected_times=(scale.TIME_BINS_US[:-1]+scale.TIME_BINS_US[1:])/2
    calls=[]
    def peak(hist,centers):
        index=scale.NOMINAL_BINS[len(calls)]
        calls.append(hist.copy())
        return 40*np.exp(-expected_times[index]/300),.2,8
    monkeypatch.setattr(scale,'_fit_direct_peak',peak)
    result=scale.fit_direct_charge_scale(segments(100))
    assert result['status']=='ok'
    assert result['intercept']==pytest.approx(40,rel=1e-5)
    assert result['tau_us']==pytest.approx(300,rel=1e-5)
    assert result['n_segments']==1900
    assert result['usable_bins']==14
    assert len(calls)==14
    assert all(row['n_segments']==100 for row in result['bin_status'])
    assert result['units']=='FLOW charge units / cm'
    assert result['curve']['x'][0]==0


def test_direct_failures_remain_null(monkeypatch):
    monkeypatch.setattr(scale,'_fit_direct_peak',lambda *_: (_ for _ in ()).throw(ValueError('unidentifiable peak')))
    result=scale.fit_direct_charge_scale(segments())
    assert result['status']=='failed'
    assert result['intercept'] is None and result['curve']['x']==[]
    assert result['usable_bins']==0
    assert len([row for row in result['bin_status'] if row['status']=='peak_fit_failed'])==14


def test_direct_bins_do_not_double_count_edges(monkeypatch):
    monkeypatch.setattr(scale,'_fit_direct_peak',lambda *_: (_ for _ in ()).throw(ValueError('no peak')))
    data=segments(1)[:len(scale.TIME_BINS_US)].copy()
    # Only19availablearrayrows; allinternaledgesplusfirstboundary areunique.
    data['t']=scale.TIME_BINS_US[:-1]/.1
    result=scale.fit_direct_charge_scale(data)
    assert sum(row['n_segments'] for row in result['bin_status'])==len(data)


def test_direct_peak_requires_resolved_four_parameter_covariance(monkeypatch):
    centers=np.linspace(10.,80.,30)
    hist=np.exp(-.5*((centers-40)/12)**2)
    monkeypatch.setattr(scale,'curve_fit',lambda *args,**kwargs:(np.array([40.,4.,5.,1.]),np.ones((4,4))))
    with pytest.raises(ValueError,match='covariance'):
        scale._fit_direct_peak(hist,centers)


def test_direct_peak_rejects_finite_boundary_parameter(monkeypatch):
    centers=np.linspace(10.,80.,30)
    hist=np.exp(-.5*((centers-40)/12)**2)
    monkeypatch.setattr(scale,'curve_fit',lambda *args,**kwargs:(np.array([40.,.01,5.,1.]),np.eye(4)))
    with pytest.raises(ValueError,match='boundary'):
        scale._fit_direct_peak(hist,centers)


def test_direct_peak_recovers_known_langau_peak_with_counting_noise():
    centers=(scale.CHARGE_BINS[:-1]+scale.CHARGE_BINS[1:])/2
    expected=scale.langau_func(centers,60.,8.,10.,1.)
    histogram=np.random.default_rng(20261104).poisson(expected*50000)/50000.
    peak,error,ncenters=scale._fit_direct_peak(histogram,centers)
    assert peak==pytest.approx(60.,abs=.5)
    assert 0 < error < .5
    assert ncenters > 4


def test_monitor_preserves_failures_and_deduplicates_copies(tmp_path):
    good,bad=receipt(),receipt()
    bad['fit_status']='failed'
    bad['period_start']='2026-11-03T00:00:00+00:00'
    bad['period_end']='2026-11-03T06:00:00+00:00'
    bad['pooled_segments']='/run/pool2/segments.h5'
    paths=[]
    for name,data in [('good',good),('copy',good),('bad',bad)]:
        path=tmp_path/(name+'.json');path.write_text(json.dumps(data));paths.append(path)
    records=scale.monitor_records(paths)
    assert len(records)==2
    assert [row['proxy']['status'] for row in records]==['ok','unavailable']
    output=tmp_path/'monitor.html';json_output=tmp_path/'monitor.json'
    scale.write_monitor(records,output,json_output)
    assert '<script src=' not in output.read_text().lower()
    assert 'charge_scale_20261102T180000Z_to_20261103T060000Z' in output.read_text()
    saved=json.loads(json_output.read_text())['records']
    assert saved[1]['proxy']['intercept'] is None


def test_monitor_separates_modes_and_recorded_cohorts(tmp_path):
    paths=[]
    for index,selection in enumerate(('stored','quality-stored','quality-stored')):
        data=receipt();data['selection']=selection
        data['pooled_segments']=f'/pool/{index}/segments.h5'
        if index==2:data['reconstruction_signature']['quality']['quality_code_sha256']='changed'
        path=tmp_path/f'{index}.json';path.write_text(json.dumps(data));paths.append(path)
    records=scale.monitor_records(paths)
    assert len(records)==3
    assert records[1]['cohort']!=records[2]['cohort']
    assert records[0]['selection']!=records[1]['selection']


def test_charge_payload_does_not_claim_direct_fit_without_pool():
    result=scale.charge_scale_payload(receipt())
    assert result['proxy']['status']=='ok'
    assert result['direct']['status']=='unavailable'
    assert result['direct']['intercept'] is None


def test_missing_pool_preserves_saved_proxy_and_records_direct_error(tmp_path):
    path=tmp_path/'summary.json';path.write_text(json.dumps(receipt()))
    record=scale.monitor_records([path],direct=True)[0]
    assert record['proxy']['status']=='ok'
    assert record['direct']['status']=='unavailable'
    assert record['direct']['input_error']


def test_direct_monitor_checks_fitted_pool_hash_before_fitting(tmp_path,monkeypatch):
    data=receipt();array=segments()
    pool=tmp_path/'segments.h5'
    key='analysis/rock_muon_segments/data'
    with h5py.File(pool,'w') as target:
        target.attrs['complete']=True
        target[key]=array
    data['pooled_segments']=str(pool)
    data['n_valid_segments']=len(array)
    data['reconstruction_signature']['reconstruction']={'run_info':{'crs_ticks':.1}}
    data['fit']['flow_datasets']=[{'path':key,'sha256':'not-the-original-pool'}]
    path=tmp_path/'summary.json';path.write_text(json.dumps(data))
    calls=[]
    def fit(rows,ticks):
        calls.append((rows.copy(),ticks))
        return dict(status='failed',intercept=None,fit_message='synthetic diagnostic')
    monkeypatch.setattr(scale,'fit_direct_charge_scale',fit)
    record=scale.monitor_records([path],direct=True)[0]
    assert record['proxy']['status']=='ok'
    assert 'SHA-256' in record['direct']['fit_message']
    assert not calls
    data['fit']['flow_datasets'][0]['sha256']=hashlib.sha256(array.tobytes()).hexdigest()
    path.write_text(json.dumps(data))
    record=scale.monitor_records([path],direct=True)[0]
    assert record['direct']['status']=='failed'
    assert len(calls)==1 and len(calls[0][0])==len(array) and calls[0][1]==.1
