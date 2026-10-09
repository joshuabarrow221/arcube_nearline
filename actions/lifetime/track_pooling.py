"""Pool complete direct-FLOW segment selections before the historical fit.

This adapter does not average fitted lifetimes, change selection, or relax
DeMario's drift-bin requirements. Callers fix calendar/cohort membership before
fitting. Incompatible reconstruction metadata and duplicate inputs are errors.
The output is a small sidecar, never a mutation of source FLOW files.
"""
import hashlib
import json
from pathlib import Path

import h5py
import numpy as np

from lifetime import SPLIT_SAMPLES, valid_segments


def json_value(value):
    """Make HDF5 attributes comparable without losing array structure."""
    if isinstance(value, dict):
        return {key: json_value(item) for key, item in value.items()}
    if isinstance(value, (np.ndarray, np.generic)):
        return json_value(value.tolist())
    if isinstance(value, (tuple, list)):
        return [json_value(item) for item in value]
    if isinstance(value, bytes):
        return value.decode('utf-8')
    return value


def reconstruction_signature(flow):
    """Record metadata relevant to geometry, calibrated charge, and drift time.

    Equal recorded metadata does not establish stable live hardware settings.
    The caller must also keep acquisition/configuration cohorts separate.
    The stored electron-lifetime resource is provenance, not a fitted input:
    the estimator below consumes segment dQ, dx, nhits and t only.
    """
    groups = ['geometry_info', 'lar_info', 'charge/calib_prompt_hits']
    signature = {name: json_value(dict(flow[name].attrs)) for name in groups}
    signature['run_info'] = {name: json_value(flow['run_info'].attrs[name])
        for name in ['crs_ticks', 'e_field', 'is_mc', 'charge_thresholds', 'data_packet_type']}
    return signature


def build_segment_pool(members, destination):
    """Concatenate each split category exactly once, retaining row lineage.

    Each member supplies flow_file, selection_file, flow_sha256, and the saved
    per-dataset segment hashes from its completed direct-FLOW fit attempt.
    Per-file fit *success* is deliberately irrelevant to membership.
    The beam/cosmic keys partition selected segments by the event's recorded
    external-trigger count. The downstream All tracks fit uses both arrays;
    its cosmic-enriched comparison reuses the zero-trigger array and is not an
    independent data sample. No additional particle-origin classifier is run.
    """
    if not members:
        raise ValueError('empty pool')
    names, hashes = set(), set()
    arrays = {sample: [] for sample, _ in SPLIT_SAMPLES}
    origins = {sample: [] for sample, _ in SPLIT_SAMPLES}
    records, signature, dtype = [], None, None
    for index, member in enumerate(members):
        source = Path(member['flow_file']).resolve()
        # A renamed byte-identical FLOW file must not acquire a second weight.
        if source in names or member['flow_sha256'] in hashes:
            raise ValueError('duplicate FLOW member')
        names.add(source); hashes.add(member['flow_sha256'])
        with h5py.File(source, 'r') as flow, h5py.File(member['selection_file'], 'r') as selected:
            if not bool(selected.attrs.get('complete', False)):
                raise ValueError('incomplete FLOW selection')
            if int(selected.attrs['events_processed']) != len(flow['charge/events/data']):
                raise ValueError('selection does not cover every FLOW event')
            if Path(selected.attrs['source_file']).resolve() != source:
                raise ValueError('selection belongs to another FLOW input')
            current = reconstruction_signature(flow)
            current['selector_sha256'] = str(selected.attrs['selector_sha256'])
            if signature is not None and current != signature:
                raise ValueError('incompatible reconstruction or selector metadata')
            signature = current
            counts = {}
            for sample, path in SPLIT_SAMPLES:
                data = selected[path][:]
                if hashlib.sha256(data.tobytes()).hexdigest() != member['segment_sha256'][path]:
                    raise ValueError('segment content differs from completed selection receipt')
                if dtype is not None and data.dtype != dtype:
                    raise ValueError('incompatible segment dtype')
                dtype = data.dtype
                arrays[sample].append(data)
                origins[sample].append(np.full(len(data), index, dtype='i4'))
                counts[sample] = dict(raw_segments=len(data), valid_segments=len(valid_segments(data)),
                    selected_tracks=int(selected[str(Path(path).parent)].attrs['selected_track_count']))
            records.append(dict(member, events_processed=int(selected.attrs['events_processed']), counts=counts))
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(destination, 'w') as pooled:
        pooled.attrs['pool_members'] = json.dumps(records)
        pooled.attrs['reconstruction_signature'] = json.dumps(signature, sort_keys=True)
        pooled.attrs['complete'] = True
        for sample, path in SPLIT_SAMPLES:
            data = np.concatenate(arrays[sample])
            pooled.create_dataset(path, data=data)
            # Segment IDs are file-local. Retain them unchanged and write a
            # parallel source index, rather than treating repeated IDs as duplicates.
            pooled.create_dataset(str(Path(path).parent)+'/member_index', data=np.concatenate(origins[sample]))
    return dict(selection_file=str(destination), members=records, reconstruction_signature=signature)


def drift_support(segments):
    """Expose all legacy-bin prerequisites, including bins outside attenuation.

    The existing fitter first fits all 19 slice histograms, then uses only its
    historical subset for attenuation. Missing early-bin support still aborts
    that procedure. This diagnostic follows the same closed interval masks;
    it never removes a bad bin or changes the fitter's acceptance.
    """
    segments = valid_segments(segments)
    edges = np.linspace(0, 1960, 20)
    rows = []
    for index, (left, right) in enumerate(zip(edges[:-1], edges[1:])):
        chosen = segments[(segments['t'] >= left) & (segments['t'] <= right)]
        nhits = chosen['nhits']/chosen['dx']
        charge = chosen['dQ']/chosen['dx']
        hit_support = int(np.sum((nhits >= 0) & (nhits <= 30)))
        charge_support = int(np.sum((charge >= 0) & (charge <= 250)))
        rows.append(dict(bin=index, start_us=float(left/10), end_us=float(right/10),
            n_segments=len(chosen), hit_histogram_entries=hit_support,
            charge_histogram_entries=charge_support,
            histogram_supported=bool(len(chosen) >= 5 and hit_support and charge_support)))
    return rows
