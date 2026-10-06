"""Read-only FLOW adapter to DeMario's existing RockMuonSelection routines.

Writes only selected segments to a separate HDF5 file, leaving the source FLOW
untouched. It uses the source's event→hit references and stored geometry.
"""
import hashlib
import inspect
from pathlib import Path
import h5py
import numpy as np
from packet_pedestal import take


def select_muons(input_file, output_file, max_events=None):
    from proto_nd_flow.selection.RockMuon_selection import RockMuonSelection
    selector = RockMuonSelection(name='lifetime_selection', classname='RockMuonSelection', data_manager=None)
    selected = {'beam': [], 'cosmic': []}
    counts = {'beam': 0, 'cosmic': 0}
    with h5py.File(input_file, 'r') as source:
        bounds = source['geometry_info'].attrs['lar_detector_bounds']
        events = source['charge/events/data']
        hits = source['charge/calib_prompt_hits/data']
        refs = source['charge/events/ref/charge/calib_prompt_hits/ref']
        regions = source['charge/events/ref/charge/calib_prompt_hits/ref_region']
        n_events = len(events) if max_events is None else min(len(events), max_events)
        for index in range(n_events):
            event = events[index]
            region = regions[index]
            links = refs[int(region['start']):int(region['stop'])]
            indices = links[links[:, 0] == index, 1]
            event_hits = take(hits, np.unique(indices))
            valid = np.isfinite(event_hits['x']) & np.isfinite(event_hits['y']) & np.isfinite(event_hits['z'])
            event_hits = event_hits[valid]
            if len(event_hits) < 3:
                continue
            sample = 'beam' if event['n_ext_trigs'] > 0 else 'cosmic'
            for indices in selector.cluster(event_hits):
                if len(indices) <= 10:
                    continue
                cluster = selector.clean_noise_hits(event_hits[indices])
                if len(cluster) < 3:
                    continue
                muon, *_ = selector.select_muon_track(cluster, bounds)
                if len(muon):
                    selector.track_count += 1
                    segments, _, _ = selector.segments(muon)
                    selected[sample].extend(segments)
                    counts[sample] += 1
            if index % 100 == 0:
                print(f'Selection: {index+1}/{n_events} events, tracks {counts}', flush=True)
        complete = n_events == len(events)
    Path(output_file).parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(output_file, 'w') as out:
        out.attrs['source_file'] = str(Path(input_file).resolve())
        out.attrs['complete'] = complete
        out.attrs['events_processed'] = n_events
        out.attrs['selector_sha256'] = hashlib.sha256(Path(inspect.getfile(RockMuonSelection)).read_bytes()).hexdigest()
        for sample, name in [('beam','beam_rock_muon_segments'),('cosmic','cosmic_muon_segments')]:
            group = out.create_group('analysis/'+name)
            array = np.array([tuple(row) for row in selected[sample]], dtype=selector.rock_muon_segments_dtype)
            group.create_dataset('data', data=array)
            group.attrs['selected_track_count'] = counts[sample]
            group.attrs['selected_segment_count'] = len(array)
    return output_file
