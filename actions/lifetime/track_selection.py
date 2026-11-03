"""Read-only FLOW adapter to DeMario's existing RockMuonSelection routines.

Writes only selected segments to a separate HDF5 file, leaving the source FLOW
untouched. It uses the source's event→hit references and stored geometry.
"""
import hashlib
import inspect
import json
import os
import tempfile
from pathlib import Path
import h5py
import numpy as np
from packet_pedestal import take


def projected_segments(selector, muon_hits, segment_length=2.):
    """Partition each io-group along its PCA coordinate, with unique hit membership.

    The installed legacy selector orders and stops on z. That loses tracks with
    constant z and double-counts hits on inclusive boundaries. This adapter uses
    half-open scalar intervals (last interval closed), including the short tail.
    It does not connect or require alignment between different io-groups.
    """
    segments, memberships = [], []
    for iog in np.unique(muon_hits['io_group']):
        hits = muon_hits[muon_hits['io_group'] == iog]
        if len(hits) < 3:
            continue
        positions = np.column_stack([hits[key] for key in ('x', 'y', 'z')])
        center = positions.mean(axis=0)
        _, _, axes = np.linalg.svd(positions-center, full_matrices=False)
        axis = axes[0]
        if axis[np.argmax(np.abs(axis))] < 0:
            axis = -axis
        projected = (positions-center) @ axis
        low, high = float(projected.min()), float(projected.max())
        if high-low <= 1e-9:
            continue
        # SVD roundoff can turn exactly 10 cm into 10+epsilon, creating
        # a spurious tail containing only the endpoint hit.
        count = max(1, int(np.ceil((high-low)/segment_length-1e-10)))
        bins = np.minimum(((projected-low)/segment_length).astype(int), count-1)
        for index in range(count):
            chosen = hits[bins == index]
            if not len(chosen):
                continue
            start = center + (low+index*segment_length)*axis
            end = center + min(low+(index+1)*segment_length, high)*axis
            dx = float(np.linalg.norm(end-start))
            if dx <= 1e-9:
                continue
            selector.segment_count += 1
            row = np.zeros((), dtype=selector.rock_muon_segments_dtype)
            row['rock_segment_id'] = selector.segment_count
            for dim, a, b in zip(('x', 'y', 'z'), start, end):
                row[dim+'_start'], row[dim+'_end'], row[dim+'_mid'] = a, b, (a+b)/2
            row['dE'], row['dQ'] = chosen['E'].sum(), chosen['Q'].sum()
            row['dx'], row['nhits'], row['io_group'] = dx, len(chosen), iog
            row['t'] = (chosen['t_drift'].min()+chosen['t_drift'].max())/2
            segments.append(row)
            memberships.append(chosen['id'])
    return segments, memberships


def select_muons(input_file, output_file, max_events=None, quality_config=None):
    """Reuse the installed selector, writing a separate segment-only sidecar.

    ``max_events`` is a development limit; the production CLI leaves it unset.
    Record the selector hash and completeness because changing reconstruction
    or processing only a prefix can change the inferred lifetime.

    Both timing categories run the same clustering, cleaning and through-going
    track selection below. "Cosmic-enriched" adds only n_ext_trigs == 0 at the
    event level; it is not a dedicated cosmic-muon identification algorithm.
    lifetime.fit_samples forms All tracks by concatenating the two categories'
    segments before fitting. Cosmic-enriched tracks are included in that fit,
    so their separate curve is a correlated subset comparison.
    """
    from proto_nd_flow.selection.RockMuon_selection import RockMuonSelection
    from stored_quality import build_quality_context, _reference, EVENTS, HITS
    if Path(input_file).resolve() == Path(output_file).resolve():
        raise ValueError('selection output must differ from source FLOW')
    if max_events is not None and max_events < 0:
        raise ValueError('max_events must be nonnegative')
    selector = RockMuonSelection(name='lifetime_selection', classname='RockMuonSelection', data_manager=None)
    selected = {'beam': [], 'cosmic': []}
    counts = {'beam': 0, 'cosmic': 0}
    origins = {'beam': [], 'cosmic': []}
    quality_records = []
    with h5py.File(input_file, 'r') as source:
        quality = build_quality_context(source, quality_config)
        bounds = source['geometry_info'].attrs['lar_detector_bounds']
        events = source['charge/events/data']
        hits = source['charge/calib_prompt_hits/data']
        refs = _reference(source, EVENTS, HITS)
        regions = source['charge/events/ref/charge/calib_prompt_hits/ref_region']
        if len(regions) != len(events):
            raise ValueError('incomplete event-to-hit reference regions')
        n_events = len(events) if max_events is None else min(len(events), max_events)
        for index in range(n_events):
            event = events[index]
            region = regions[index]
            # Explicit H5Flow references are the association authority. Hit
            # rows are not assumed contiguous or aligned with event row IDs.
            start, stop = int(region['start']), int(region['stop'])
            if start < 0 or stop < start or stop > len(refs):
                raise ValueError('invalid event-to-hit reference region')
            links = refs[start:stop]
            indices = links[links[:, 0] == index, 1]
            if (np.any(indices < 0) or np.any(indices >= len(hits))
                    or len(np.unique(indices)) != len(indices)):
                raise ValueError('invalid or duplicate event-to-hit reference')
            event_hits = take(hits, np.sort(indices))
            if len(np.unique(event_hits['id'])) != len(event_hits):
                raise ValueError('duplicate event hit IDs')
            keep, diagnostic = quality.evaluate(event_hits)
            diagnostic.update(event_row=index, event_id=int(event['id']), rejected_segments=0)
            quality_records.append(diagnostic)
            bad_ids = set(map(int, event_hits['id'][~keep]))
            valid = np.isfinite(event_hits['x']) & np.isfinite(event_hits['y']) & np.isfinite(event_hits['z'])
            event_hits = event_hits[valid]
            if len(event_hits) < 3:
                continue
            # Compatibility keys describe trigger timing only. Neither this
            # branch nor the selector queries IFBeam or proves particle origin.
            # Zero external triggers is not a verified beam-off interval; a
            # positive count can also accompany cosmics during commissioning.
            sample = 'beam' if event['n_ext_trigs'] > 0 else 'cosmic'
            for indices in selector.cluster(event_hits):
                if len(indices) <= 10:
                    continue
                cluster = selector.clean_noise_hits(event_hits[indices])
                if len(cluster) < 3:
                    continue
                # Keep the local selector's through-going criteria, supplying
                # stored xyz bounds directly. Its legacy select_muon_track
                # wrapper incorrectly flips those bounds as if they were zyx.
                candidate = selector.clean_noise_hits(cluster)
                if len(candidate) < 3:
                    continue
                length, *_ = selector.length(candidate)
                accepted_track = (length >= selector.length_cut
                    and selector.average_distance(candidate) <= 1.5
                    and selector.close_to_two_faces(np.asarray(bounds).reshape(-1), candidate))
                if accepted_track:
                    selector.track_count += 1
                    segments, memberships = projected_segments(selector, candidate)
                    accepted = 0
                    for segment, member_ids in zip(segments, memberships):
                        # Veto whole segments, preserving measured Q, N and dx.
                        # Subtracting flagged hits would bias the charge density.
                        if any(int(hit_id) in bad_ids for hit_id in member_ids):
                            diagnostic['rejected_segments'] += 1
                            continue
                        selected[sample].append(segment)
                        origins[sample].append((index, int(event['id']), selector.track_count))
                        accepted += 1
                    counts[sample] += int(accepted > 0)
            if index % 100 == 0:
                print(f'Selection: {index+1}/{n_events} events, tracks {counts}', flush=True)
        complete = n_events == len(events)
    Path(output_file).parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(suffix='.h5.tmp', dir=Path(output_file).parent)
    os.close(fd)
    try:
        with h5py.File(temporary, 'w') as out:
            # These sidecars belong outside the watched FLOW tree; otherwise a
            # recursive watcher could enqueue them as new detector input files.
            out.attrs['source_file'] = str(Path(input_file).resolve())
            out.attrs['complete'] = complete
            out.attrs['events_processed'] = n_events
            out.attrs['selector_sha256'] = hashlib.sha256(Path(inspect.getfile(RockMuonSelection)).read_bytes()).hexdigest()
            out.attrs.update(quality.metadata)
            out.attrs['adapter_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
            out.attrs['selection_origin'] = 'fresh local selector; scalar-projection segments; whole-segment quality veto'
            out.create_dataset('quality/events/data', data=np.array([json.dumps(row, sort_keys=True) for row in quality_records], dtype=h5py.string_dtype()))
            for sample, name in [('beam','beam_rock_muon_segments'),('cosmic','cosmic_muon_segments')]:
                group = out.create_group('analysis/'+name)
                array = np.array(selected[sample], dtype=selector.rock_muon_segments_dtype)
                group.create_dataset('data', data=array)
                group.create_dataset('source_event_track', data=np.array(origins[sample], dtype=[('event_row','i8'),('event_id','i8'),('track_id','i8')]))
                group.attrs['selected_track_count'] = counts[sample]
                group.attrs['selected_segment_count'] = len(array)
        os.replace(temporary, output_file)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)
    return output_file
