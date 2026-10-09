"""Use the exact static calibration named by FLOW, with an explicit identity check.

This optional adapter avoids reading all FLOW hits just to recover intercepts.
It does not use nearline panel means. The ADC conversion constants are accepted
only after verification against stratified, explicitly linked Q_raw/ADC pairs
in that particular FLOW file. Unknown calibration channels remain missing.
"""
import hashlib
import json
from pathlib import Path

import numpy as np

from flow_input import open_flow
from packet_pedestal import take, key


def tile_lookup(group, io, channel):
    """Decode the persisted two-key ndlar_flow LUT, including its validity mask.

    ndlar_flow.util.lut uses column-major integer strides and reserves row zero
    for the default. Decode only this documented tile-id layout; reject schema
    changes rather than guessing an address from IO-channel number.
    """
    bounds = group.attrs['meta'][0]['min_max_keys']
    if bounds.shape != (2,2): raise ValueError('unsupported geometry tile LUT')
    data = group['data'][:]
    io, channel = np.asarray(io,dtype=np.int64), np.asarray(channel,dtype=np.int64)
    valid = (io>=bounds[0,0]) & (io<=bounds[0,1]) & (channel>=bounds[1,0]) & (channel<=bounds[1,1])
    indices = 1+io-bounds[0,0]+(channel-bounds[1,0])*(bounds[0,1]-bounds[0,0]+1)
    indices = np.where(valid,indices,0)
    return np.where(valid & data['filled'][indices], data['data'][indices], -1)


def calibration_map(flow_file, pedestal_file, channels):
    """Return per-channel ADC pedestals plus source and verification metadata."""
    payload = Path(pedestal_file).read_bytes(); calibration = json.loads(payload)
    vref, vcm, adc_counts, gain = 1568., 478.1, 256., 4.522
    slope = (vref-vcm)/(adc_counts*gain)
    with open_flow(flow_file) as (flow, provenance):
        attrs = flow['charge/calib_prompt_hits'].attrs
        if Path(str(attrs.get('pedestal_file',''))).name != Path(pedestal_file).name:
            raise ValueError('pedestal filename does not match FLOW calibration provenance')
        if attrs.get('configuration_file','') or attrs.get('gain_file',''):
            raise ValueError('non-default channel configuration/gain needs a dedicated adapter')
        geometry = flow['geometry_info/tile_id']
        def pedestals(p):
            tile = tile_lookup(geometry,p['io_group'],p['io_channel'])
            uid = p['io_group'].astype(np.int64)*1_000_000_000+tile*100_000+p['chip_id'].astype(np.int64)*100+p['channel_id']
            mv = np.array([calibration.get(str(int(u)),{}).get('pedestal_mv',np.nan) if t>=0 else np.nan for u,t in zip(uid,tile)])
            return (mv-vcm)*adc_counts/(vref-vcm)
        refs = flow['charge/calib_prompt_hits/ref/charge/packets/ref']
        verified = 0; maximum_residual = 0.; seen_references = set()
        # Stratify across the file, not only the first event or an attractive
        # drift bin. Verify the same charge-calibration identity everywhere.
        for start in np.linspace(0,max(0,len(refs)-2000),3,dtype=int):
            links = refs[start:start+2000]
            # Small files have overlapping strata. Count each linked row once
            # so duplicate reads cannot satisfy the minimum verification size.
            indices = range(start,start+len(links))
            unique = np.array([i not in seen_references for i in indices],dtype=bool)
            seen_references.update(indices)
            links = links[unique]
            if not len(links): continue
            h = take(flow['charge/calib_prompt_hits/data'],links[:,0])
            p = take(flow['charge/packets/data'],links[:,1])
            pedestal = pedestals(p)
            valid = np.isfinite(pedestal) & np.isfinite(h['Q_raw']) & (p['packet_type']==0)
            for field in ('io_group','io_channel','chip_id','channel_id'):
                valid &= h[field]==p[field]
            residual = np.abs(h['Q_raw'][valid]-slope*(p['dataword'][valid].astype(float)-pedestal[valid]))
            if len(residual): maximum_residual = max(maximum_residual,float(residual.max()))
            verified += int(valid.sum())
        if verified < 100 or maximum_residual > 1e-5:
            raise ValueError(f'static calibration identity not verified: pairs={verified}, max residual={maximum_residual}')
        pedestal = pedestals(channels)
        valid = np.isfinite(pedestal) & (pedestal>=0) & (pedestal<=255)
        address = key(channels['io_group'],channels['io_channel'],channels['chip_id'],channels['channel_id'])
        mapping = {int(k):float(p) for k,p in zip(address[valid],pedestal[valid])}
    return mapping, dict(pedestal_file=str(Path(pedestal_file).resolve()),sha256=hashlib.sha256(payload).hexdigest(),
        verified_pairs=verified, maximum_Q_raw_residual=maximum_residual,flow_provenance=provenance,
        adc_constants=dict(vref_mv=vref,vcm_mv=vcm,adc_counts=adc_counts,gain=gain))
