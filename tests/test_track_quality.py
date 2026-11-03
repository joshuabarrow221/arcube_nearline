"""Synthetic topology/clock tests; these do not calibrate detector cuts."""
from dataclasses import replace
import json
from pathlib import Path
import sys

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'actions/lifetime'))
from track_quality import TrackQualityConfig, quality_hit_mask


DTYPE = [(name, 'f8') for name in ('x', 'y', 'z', 'ts_pps', 't_drift')] + [
    (name, 'i8') for name in ('io_group', 'io_channel', 'chip_id', 'channel_id')]


def hits(n):
    result = np.zeros(n, dtype=DTYPE)
    result['io_group'] = result['io_channel'] = 1
    result['chip_id'] = 11
    result['channel_id'] = np.arange(n) % 64
    result['ts_pps'] = np.arange(n) * 10
    result['x'] = np.arange(n) * .2
    result['y'] = (np.arange(n) % 8) * .44
    result['z'] = (np.arange(n) // 8) * .44 + 1300
    return result


def test_hot_pixel_requires_repetition_rate_and_drift_streak():
    data = hits(20)
    data['channel_id'] = 3
    data['y'] = 1
    data['z'] = 1300
    keep, diagnostic = quality_hit_mask(data)
    assert not keep.any()
    assert diagnostic['hot_pixel_hits'] == 20
    assert diagnostic['coherent_noise_hits'] == 0
    # Repeated ADC hits without a drift streak are not an axis-track candidate.
    data['x'] = 0
    assert quality_hit_mask(data)[0].all()
    data['x'] = np.arange(20) * .2
    data['ts_pps'] = np.arange(20) * 1000
    assert quality_hit_mask(data)[0].all()


def test_horizontal_multi_pixel_track_and_isolated_hits_survive():
    data = hits(20)
    data['y'] = 0
    data['z'] = 1300
    assert quality_hit_mask(data)[0].all()
    assert quality_hit_mask(hits(1))[0].all()


def test_chip_burst_uses_arrival_time_and_crosses_fixed_window_edges():
    data = hits(64)
    data['ts_pps'] = np.linspace(490, 590, 64)
    # Drift time is not the electronics coincidence clock.
    data['t_drift'] = np.arange(64) * 1000
    keep, diagnostic = quality_hit_mask(data)
    assert not keep.any()
    assert diagnostic['burst_groups_flagged']['chip'] == 1
    assert diagnostic['hot_pixel_hits'] == 0
    data['ts_pps'] = np.arange(64) * 1000
    data['t_drift'] = 0
    assert quality_hit_mask(data)[0].all()
    json.dumps(diagnostic, allow_nan=False)


def test_uart_aliases_are_one_physical_pixel_and_chip():
    data = hits(64)
    data['io_channel'] = np.arange(64) % 4 + 1
    assert not quality_hit_mask(data)[0].any()
    data['channel_id'] = np.arange(64) % 8
    data['x'] = 0
    assert quality_hit_mask(data)[0].all()


def test_sustained_bursts_are_all_marked_and_other_module_survives():
    burst = hits(192)
    burst['ts_pps'] = np.repeat([0, 450, 900], 64) + np.tile(np.arange(64), 3)
    burst['y'] = np.tile((np.arange(64) % 8) * .44, 3)
    burst['z'] = np.tile((np.arange(64) // 8) * .44 + 1300, 3)
    clean = hits(10)
    clean['io_group'] = 3
    keep, _ = quality_hit_mask(np.concatenate([burst, clean]))
    assert not keep[:192].any()
    assert keep[192:].all()


def test_tile_and_module_bursts_can_be_detected_without_a_chip_trigger():
    data = hits(128)
    data['chip_id'] = np.arange(128) // 32
    data['channel_id'] = np.arange(128) % 32
    data['ts_pps'] = np.arange(128)
    cfg = replace(TrackQualityConfig(), chip_min_hits=1000)
    keep, diagnostic = quality_hit_mask(data, cfg)
    assert not keep.any()
    assert diagnostic['burst_groups_flagged']['chip'] == 0
    assert diagnostic['burst_groups_flagged']['tile'] == 1
    data = np.tile(data, 4)
    data['io_channel'] = np.repeat([1, 5, 9, 13], 128)
    cfg = replace(cfg, tile_min_hits=1000)
    keep, diagnostic = quality_hit_mask(data, cfg)
    assert not keep.any()
    assert diagnostic['burst_groups_flagged']['module'] == 1


def test_dense_straight_track_is_not_a_spatially_broad_burst():
    data = hits(128)
    data['chip_id'] = np.arange(128) // 32
    data['channel_id'] = np.arange(128) % 32
    data['ts_pps'] = np.arange(128)
    data['y'] = np.arange(128) * .44
    data['z'] = 1300 + data['y'] * .2
    assert quality_hit_mask(data)[0].all()


def test_geometry_ids_override_two_by_two_fallback():
    data = hits(64)
    data['io_channel'] = np.arange(64) // 16 + 1
    # Explicit distinct tiles must not be conflated into the fallback tile.
    assert quality_hit_mask(data, tile_ids=data['io_channel'])[0].all()
    with pytest.raises(ValueError, match='invalid geometry'):
        quality_hit_mask(data, tile_ids=np.full(64, -1))


def test_clock_scale_changes_window_and_optional_rollover_is_explicit():
    data = hits(64)
    data['ts_pps'] = np.arange(64) * 10
    assert not quality_hit_mask(data)[0].any()
    assert quality_hit_mask(data, replace(TrackQualityConfig(), tick_us=10))[0].all()
    period = 10_000_000
    data['ts_pps'] = (np.arange(64) + period - 32) % period
    keep, diagnostic = quality_hit_mask(data, replace(TrackQualityConfig(), rollover_ticks=period))
    assert not keep.any()
    assert diagnostic['timestamps_unwrapped']


def test_empty_invalid_and_missing_field_handling():
    keep, diagnostic = quality_hit_mask(hits(0))
    assert len(keep) == 0 and diagnostic['rejected_hits'] == 0
    data = hits(1)
    data['ts_pps'] = np.nan
    keep, diagnostic = quality_hit_mask(data)
    assert not keep.any() and diagnostic['invalid_hits'] == 1
    with pytest.raises(ValueError, match='requires fields'):
        quality_hit_mask(np.zeros(1, dtype=[('x', 'f8')]))


@pytest.mark.parametrize('kwargs', [dict(window_us=0), dict(tick_us=np.nan),
    dict(hot_min_hits=1), dict(chip_min_channels=1.5), dict(tile_min_chips=True),
    dict(hot_min_axis_cosine=1.01), dict(drift_axis='t'), dict(rollover_ticks=-1)])
def test_invalid_config_rejected(kwargs):
    with pytest.raises(ValueError):
        TrackQualityConfig(**kwargs)
