"""Synthetic checks of pooling membership/provenance, not detector evidence."""
from pathlib import Path
import hashlib
import sys

import h5py
import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'actions/lifetime'))
from lifetime import SPLIT_SAMPLES
from track_pooling import build_segment_pool, drift_support


def member(directory, name, tick=10.):
    flow, sidecar = directory/(name+'.h5'), directory/(name+'.segments.h5')
    with h5py.File(flow, 'w') as f:
        f.create_dataset('charge/events/data', data=[0])
        for group in ['geometry_info', 'lar_info', 'charge/calib_prompt_hits']:
            f.create_group(group)
        run = f.create_group('run_info')
        run.attrs.update(crs_ticks=.1, e_field=.05, is_mc=False,
                         charge_thresholds='medm', data_packet_type=0)
    data = np.array([(1, 1., 10., 3., tick)],
        dtype=[('rock_segment_id','i4'),('dx','f8'),('dQ','f8'),('nhits','f8'),('t','f8')])
    hashes = {}
    with h5py.File(sidecar, 'w') as f:
        f.attrs.update(source_file=str(flow.resolve()), complete=True,
                       events_processed=1, selector_sha256='unchanged')
        for sample,path in SPLIT_SAMPLES:
            f.create_dataset(path, data=data if sample=='cosmic' else data[:0])
            f[str(Path(path).parent)].attrs['selected_track_count'] = int(sample=='cosmic')
            hashes[path] = hashlib.sha256(f[path][:].tobytes()).hexdigest()
    return dict(flow_file=str(flow), selection_file=str(sidecar),
                flow_sha256=hashlib.sha256(name.encode()).hexdigest(), segment_sha256=hashes)


def test_concatenation_keeps_file_local_ids_and_empty_categories(tmp_path):
    members = [member(tmp_path,'a'), member(tmp_path,'b',tick=150.)]
    output = tmp_path/'pool.h5'
    build_segment_pool(members, output)
    with h5py.File(output,'r') as f:
        # Repeated file-local IDs must survive; lineage makes them unambiguous.
        np.testing.assert_array_equal(f['analysis/cosmic_muon_segments/data']['rock_segment_id'], [1,1])
        np.testing.assert_array_equal(f['analysis/cosmic_muon_segments/member_index'][:], [0,1])
        np.testing.assert_array_equal(f['analysis/cosmic_muon_segments/data']['t'], [10.,150.])
        assert len(f['analysis/beam_rock_muon_segments/data']) == 0


def test_duplicate_source_cannot_acquire_extra_weight(tmp_path):
    first = member(tmp_path,'a')
    with pytest.raises(ValueError, match='duplicate FLOW'):
        build_segment_pool([first,first],tmp_path/'pool.h5')


@pytest.mark.parametrize('failure', ['incomplete', 'selector', 'calibration', 'tampered'])
def test_bad_provenance_cannot_enter_pool(tmp_path, failure):
    first, second = member(tmp_path,'a'), member(tmp_path,'b')
    with h5py.File(second['selection_file'],'r+') as f:
        if failure=='incomplete':f.attrs['complete']=False
        if failure=='selector':f.attrs['selector_sha256']='changed'
        if failure=='tampered':
            path='analysis/cosmic_muon_segments/data';data=f[path][:]
            data['dQ'] *= 2;f[path][:]=data
    if failure=='calibration':
        with h5py.File(second['flow_file'],'r+') as f:
            f['charge/calib_prompt_hits'].attrs['gain_file']='different.json'
    with pytest.raises(ValueError):
        build_segment_pool([first,second],tmp_path/'pool.h5')


def test_drift_support_retains_required_early_bin(tmp_path):
    first = member(tmp_path,'a')
    with h5py.File(first['selection_file'],'r') as f:
        data=f['analysis/cosmic_muon_segments/data'][:]
    one=drift_support(data)
    five=drift_support(np.repeat(data,5))
    assert len(one)==19 and one[0]['n_segments']==1
    assert not one[0]['histogram_supported']
    assert five[0]['histogram_supported']
    assert five[1]['n_segments']==0 and not five[1]['histogram_supported']
