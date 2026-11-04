"""Synthetic topology/clock tests; these do not calibrate detector cuts."""
from dataclasses import replace
import json
from pathlib import Path
import sys

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'actions/lifetime'))
from track_quality import TrackQualityConfig, quality_hit_mask, _burst_mask


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


def test_offset_clock_groups_do_not_require_one_short_global_time_interval():
    period = 10_000_000
    data = np.concatenate([hits(10) for _ in range(3)])
    for group, offset in zip(range(3), [.1, .45, .8]):
        rows = slice(group*10, (group+1)*10)
        data['io_group'][rows] = group+1
        data['ts_pps'][rows] += period*offset
    keep, diagnostic = quality_hit_mask(data, replace(TrackQualityConfig(), rollover_ticks=period))
    assert keep.all()
    assert diagnostic['ambiguous_clock_groups'] == 0
    assert len(diagnostic['clock_domains']) == 3


def test_broad_clock_group_is_vetoed_without_losing_an_independent_group():
    period = 10_000_000
    ambiguous, clean = hits(4), hits(10)
    ambiguous['ts_pps'] = np.array([0, .35, .7, .95])*period
    clean['io_group'] = 3
    keep, diagnostic = quality_hit_mask(np.concatenate([ambiguous, clean]),
                                       replace(TrackQualityConfig(), rollover_ticks=period))
    assert not keep[:4].any() and keep[4:].all()
    assert diagnostic['ambiguous_clock_groups'] == 1
    assert diagnostic['ambiguous_clock_hits'] == 4
    # Timing ambiguity is not evidence of a particular electronics-noise type.
    assert diagnostic['hot_pixel_hits'] == diagnostic['coherent_noise_hits'] == 0
    assert diagnostic['clock_domains'][0]['status'] == 'ambiguous'
    json.dumps(diagnostic, allow_nan=False)


def test_independent_clock_origins_do_not_manufacture_a_module_burst():
    period = 10_000_000
    data = hits(32)
    data['io_group'] = np.repeat([1, 2], 16)
    data['io_channel'] = np.tile(np.repeat([1, 5], 8), 2)
    data['ts_pps'][:16] = np.arange(16)
    data['ts_pps'][16:] = .4*period + np.arange(16)
    cfg = replace(TrackQualityConfig(), rollover_ticks=period,
                  chip_min_hits=1000, tile_min_hits=1000,
                  module_min_hits=16, module_min_channels=24, module_min_tiles=2,
                  burst_min_transverse_rms_cm=.3)
    keep, diagnostic = quality_hit_mask(data, cfg)
    assert keep.all()
    assert diagnostic['burst_groups_flagged']['module'] == 0
    assert diagnostic['maximum_distinct_channels']['module'] == 16


def test_exact_half_period_is_ambiguous_and_rollover_burst_survives_elsewhere():
    period = 10_000_000
    ambiguous, burst = hits(2), hits(64)
    ambiguous['ts_pps'] = [0, period/2]
    burst['io_group'] = 3
    burst['ts_pps'] = (np.arange(64)+period-32) % period
    keep, diagnostic = quality_hit_mask(np.concatenate([ambiguous, burst]),
                                       replace(TrackQualityConfig(), rollover_ticks=period))
    assert not keep.any()
    assert diagnostic['ambiguous_clock_hits'] == 2
    assert diagnostic['coherent_noise_hits'] == 64
    assert diagnostic['timestamps_unwrapped']


@pytest.mark.parametrize('kwargs', [dict(window_us=0), dict(tick_us=np.nan),
    dict(hot_min_hits=1), dict(chip_min_channels=1.5), dict(tile_min_chips=True),
    dict(hot_min_axis_cosine=1.01), dict(drift_axis='t'), dict(rollover_ticks=-1)])
def test_invalid_config_rejected(kwargs):
    with pytest.raises(ValueError):
        TrackQualityConfig(**kwargs)


@pytest.mark.parametrize('arrival_pattern', ['coincident', 'sliding'])
def test_batched_burst_windows_match_independent_window_enumeration(arrival_pattern):
    """Cross the batch boundary, including many windows sharing one stop."""
    count = 1350
    rng = np.random.default_rng(20261104)
    pixels = rng.integers(0, 64, count)
    coordinates = np.column_stack((np.arange(64) % 8, np.arange(64) // 8)).astype(float)
    times = np.zeros(count) if arrival_pattern == 'coincident' else np.arange(count) * .1
    permutation = rng.permutation(count)
    times, pixels = times[permutation], pixels[permutation]
    groups = np.zeros(count, dtype='i8')
    subunits = pixels // 16
    config = replace(TrackQualityConfig(), window_us=50.)
    expected = np.zeros(count, dtype=bool)
    maximum = 0
    ordered = np.argsort(times, kind='stable')
    for left in range(count):
        window = ordered[left:][times[ordered[left:]] <= times[ordered[left]] + config.window_us]
        unique = np.unique(pixels[window])
        maximum = max(maximum, len(unique))
        if len(window) < 32 or len(unique) < 24 or len(np.unique(subunits[window])) < 3:
            continue
        xy = coordinates[unique]
        mean = xy.mean(axis=0)
        covariance = (xy[:, :, None] * xy[:, None, :]).mean(axis=0) - np.outer(mean, mean)
        if np.linalg.eigvalsh(covariance)[0] >= config.burst_min_transverse_rms_cm ** 2:
            expected[window] = True
    actual, flagged_groups, actual_maximum = _burst_mask(
        times, pixels, groups, subunits, coordinates, 32, 24, 3, config)
    assert np.array_equal(actual, expected)
    assert flagged_groups == int(expected.any())
    assert actual_maximum == maximum


@pytest.mark.parametrize('threshold', [np.nextafter(.5, 0.), .5, np.nextafter(.5, 1.)])
def test_batched_burst_rms_boundary_remains_inclusive(threshold):
    # Exact binary coordinates give minor variance exactly 0.25. Repeated
    # windows exercise both full batches and the partial trailing batch.
    coordinates = np.array([[-1., -.5], [-1., .5], [1., -.5], [1., .5]])
    pixels = np.tile(np.arange(4), 300)
    times = np.zeros(len(pixels))
    config = replace(TrackQualityConfig(), burst_min_transverse_rms_cm=threshold)
    mask, groups, maximum = _burst_mask(times, pixels, np.zeros(len(pixels)),
        pixels, coordinates, 4, 4, 4, config)
    assert mask.all() == (threshold <= .5)
    assert mask.any() == (threshold <= .5)
    assert groups == int(threshold <= .5)
    assert maximum == 4
