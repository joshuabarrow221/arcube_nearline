"""Recover static ADC baselines from linked FLOW Q_raw and native ADC values.

CalibHitBuilder.charge_from_dataword is linear: Q_raw = k*(ADC-P).
Fit that line per physical channel, without guessing gain/vref/vcm constants.
Channels with only one ADC value or non-linear Q_raw are not identifiable and
are omitted. Reconstructed t0, Q, E, and track selection are not used.
"""
import h5py
import numpy as np


def key(io, channel, chip, pixel):
    """Pack the complete physical channel address; never merge across IOs."""
    return (((np.asarray(io, dtype=np.int64)*256 + channel)*256 + chip)*64 + pixel)


def take(dataset, indices):
    """Read arbitrary HDF5 row references, preserving order and repetitions."""
    indices = np.asarray(indices, dtype=np.int64)
    unique, inverse = np.unique(indices, return_inverse=True)
    # Read a contiguous span when references are local; h5py point selection
    # over hundreds of thousands of indices is otherwise very slow.
    if not len(unique):
        return np.empty(0, dtype=dataset.dtype)
    if unique[-1]-unique[0] < max(4*len(unique), 10000):
        return dataset[unique[0]:unique[-1]+1][indices-unique[0]]
    return dataset[unique][inverse]


def flow_ped_map(path, chunk=200_000):
    """Infer static ADC intercepts with chunked per-channel sufficient statistics.

    Only linked Q_raw/ADC pairs from the same channel qualify. The tight linear
    residual check tests the calibration identity, not detector noise: both
    numbers were generated from the same packet. This is not a pedestal fit
    to the observed physics charge distribution.
    """
    totals = {}
    hit_path = 'charge/calib_prompt_hits'
    with h5py.File(path, 'r') as flow:
        hits = flow[hit_path+'/data']
        packets = flow['charge/packets/data']
        refs = flow[hit_path+'/ref/charge/packets/ref']
        for start in range(0, len(refs), chunk):
            links = refs[start:start+chunk]
            h, p = take(hits, links[:, 0]), take(packets, links[:, 1])
            valid = np.isfinite(h['Q_raw']) & (p['packet_type'] == 0)
            for field in ('io_group', 'io_channel', 'chip_id', 'channel_id'):
                valid &= h[field] == p[field]
            h, p = h[valid], p[valid]
            channels = key(h['io_group'], h['io_channel'], h['chip_id'], h['channel_id'])
            unique, inverse = np.unique(channels, return_inverse=True)
            x, y = p['dataword'].astype(float), h['Q_raw'].astype(float)
            sums = np.column_stack([np.bincount(inverse, weights=w, minlength=len(unique))
                                    for w in (np.ones(len(x)), x, y, x*x, x*y, y*y)])
            for channel, values in zip(unique, sums):
                if int(channel) in totals:
                    totals[int(channel)] += values
                else:
                    totals[int(channel)] = values
    result = {}
    for channel, (n, sx, sy, sxx, sxy, syy) in totals.items():
        # Analytic OLS needs ADC variation. Constant-ADC channels cannot
        # independently determine slope and pedestal and must be omitted.
        variance = sxx-sx*sx/n
        if n < 3 or variance <= 1e-8:
            continue
        slope = (sxy-sx*sy/n)/variance
        if not np.isfinite(slope) or slope <= 0:
            continue
        intercept = (sy-slope*sx)/n
        residual = max(0., syy-sy*sy/n - slope*slope*variance)
        if residual/n > 1e-6:
            continue
        pedestal = -intercept/slope
        if 0 <= pedestal <= 255:
            result[channel] = float(pedestal)
    if not result:
        raise ValueError('FLOW references cannot determine any static per-channel ADC pedestal')
    return result
