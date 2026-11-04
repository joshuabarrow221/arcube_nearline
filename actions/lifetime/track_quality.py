"""Event-local electronics noise evidence for lifetime segment selection.

These are prospective, configurable cuts, not detector-calibrated efficiencies.
Use the mask to veto affected *whole segments*, rather than subtracting rejected
hits from a segment's charge. Evaluate the original event hits before PCA or
track cleaning, which could otherwise hide the surrounding noisy electronics.

Arrival ``ts_pps`` is in charge-clock ticks, independent of reconstructed drift
time. Callers should supply ``tick_us`` from FLOW ``run_info/crs_ticks``. Calls
must contain one event: PPS-relative timestamps cannot identify acquisition time
across unrelated events. Rollover handling is independent for each io_group;
an ambiguous group's hits are rejected without invalidating other groups or
the entire source file. Noise split across event boundaries is not joined.
"""
from collections import deque
from dataclasses import asdict, dataclass

import numpy as np


QUALITY_VERSION = 'event-electronics-v2-iogroup-clock'


class ArrivalTimeError(ValueError):
    """A clock group's PPS phases do not determine a short arrival interval."""


@dataclass(frozen=True)
class TrackQualityConfig:
    """Provisional 2x2 thresholds; lengths cm, arrival windows microseconds.

    A hot-pixel candidate needs repeated activity *and* a drift-axis streak.
    A real track confined to one pixel is not distinguishable by this evidence;
    the cut deliberately sacrifices that ambiguous topology. Burst thresholds
    count physical channels, not repeated UART aliases. Their pixel-plane minor
    RMS guard protects approximately straight tracks with high occupancy.
    """
    window_us: float = 50.0
    tick_us: float = 0.1
    drift_axis: str = 'x'
    rollover_ticks: float | None = None
    hot_min_hits: int = 12
    hot_min_rate_hz: float = 100_000.0
    hot_min_drift_span_cm: float = 2.0
    hot_max_transverse_span_cm: float = 0.5
    hot_min_axis_cosine: float = 0.995
    chip_min_hits: int = 32
    chip_min_channels: int = 24
    tile_min_hits: int = 128
    tile_min_channels: int = 96
    tile_min_chips: int = 4
    module_min_hits: int = 512
    module_min_channels: int = 384
    module_min_tiles: int = 4
    burst_min_transverse_rms_cm: float = 0.5

    def __post_init__(self):
        if self.drift_axis not in ('x', 'y', 'z'):
            raise ValueError('drift_axis must be x, y, or z')
        count_fields = ('hot_min_hits', 'chip_min_hits', 'chip_min_channels',
                        'tile_min_hits', 'tile_min_channels', 'tile_min_chips',
                        'module_min_hits', 'module_min_channels', 'module_min_tiles')
        for field in count_fields:
            value = getattr(self, field)
            if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)) or value < 1:
                raise ValueError(f'{field} must be a positive integer')
        if self.hot_min_hits < 2:
            raise ValueError('hot_min_hits must be at least two')
        for field in ('window_us', 'tick_us', 'hot_min_rate_hz', 'hot_min_drift_span_cm',
                      'hot_max_transverse_span_cm', 'hot_min_axis_cosine',
                      'burst_min_transverse_rms_cm'):
            value = getattr(self, field)
            if not np.isfinite(value) or value <= 0:
                raise ValueError(f'{field} must be finite and positive')
        if self.hot_min_axis_cosine > 1:
            raise ValueError('hot_min_axis_cosine must be at most one')
        if self.rollover_ticks is not None:
            if not np.isfinite(self.rollover_ticks) or self.rollover_ticks <= 0:
                raise ValueError('rollover_ticks must be finite and positive')
            if self.window_us >= self.rollover_ticks * self.tick_us / 2:
                raise ValueError('window_us must be less than half the rollover period')

    def to_dict(self):
        return {key: value.item() if isinstance(value, np.generic) else value
                for key, value in asdict(self).items()}


def _labels(columns):
    """Dense collision-free integer identities from electronics columns."""
    return np.unique(np.column_stack(columns), axis=0, return_inverse=True)[1]


def _groups(labels):
    order = np.argsort(labels, kind='stable')
    edges = np.flatnonzero(np.diff(labels[order])) + 1
    return np.split(order, edges) if len(order) else []


def _arrival_times(ticks, config):
    """Resolve a short interval within one io_group, never across clocks."""
    ticks = ticks.astype('f8', copy=True)
    unwrapped = False
    if len(ticks) and config.rollover_ticks is not None:
        period = config.rollover_ticks
        if np.any((ticks < 0) | (ticks >= period)):
            raise ArrivalTimeError('ts_pps lies outside the specified rollover period')
        if np.ptp(ticks) >= period / 2:
            ordered = np.sort(np.unique(ticks))
            gaps = np.diff(np.r_[ordered, ordered[0] + period])
            gap = int(np.argmax(gaps))
            start = ordered[(gap + 1) % len(ordered)]
            ticks[ticks < start] += period
            if np.ptp(ticks) >= period / 2:
                raise ArrivalTimeError('io_group does not fit within half the rollover period')
            unwrapped = True
    if len(ticks):
        ticks -= ticks.min()
    return ticks * config.tick_us, unwrapped


def _hot_mask(times, xyz, pixel_labels, config):
    rejected = np.zeros(len(times), dtype=bool)
    flagged = 0
    axis = ('x', 'y', 'z').index(config.drift_axis)
    transverse = [i for i in range(3) if i != axis]
    for indices in _groups(pixel_labels):
        if len(indices) < config.hot_min_hits:
            continue
        indices = indices[np.argsort(times[indices], kind='stable')]
        local_time = times[indices]
        positions = xyz[indices]
        end = np.searchsorted(local_time, local_time + config.window_us, side='right')
        coverage = np.zeros(len(indices) + 1, dtype=np.int64)
        minima, maxima = [deque() for _ in range(3)], [deque() for _ in range(3)]
        added = 0
        for left, right in enumerate(end):
            # A continuously firing pixel can put most of an event in one
            # window. Monotone queues avoid rescanning that window per hit.
            while added < right:
                for coordinate in range(3):
                    while minima[coordinate] and positions[minima[coordinate][-1], coordinate] >= positions[added, coordinate]:
                        minima[coordinate].pop()
                    while maxima[coordinate] and positions[maxima[coordinate][-1], coordinate] <= positions[added, coordinate]:
                        maxima[coordinate].pop()
                    minima[coordinate].append(added)
                    maxima[coordinate].append(added)
                added += 1
            for coordinate in range(3):
                while minima[coordinate][0] < left:
                    minima[coordinate].popleft()
                while maxima[coordinate][0] < left:
                    maxima[coordinate].popleft()
            count = right - left
            if count < config.hot_min_hits:
                continue
            span_us = max(local_time[right - 1] - local_time[left], config.tick_us)
            if (count - 1) * 1e6 / span_us < config.hot_min_rate_hz:
                continue
            extent = np.array([positions[maxima[c][0], c] - positions[minima[c][0], c]
                               for c in range(3)])
            drift = extent[axis]
            transverse_span = np.linalg.norm(extent[transverse])
            norm = np.linalg.norm(extent)
            if (drift >= config.hot_min_drift_span_cm
                    and transverse_span <= config.hot_max_transverse_span_cm
                    and norm > 0 and drift / norm >= config.hot_min_axis_cosine):
                coverage[left] += 1
                coverage[right] -= 1
        local_bad = np.cumsum(coverage[:-1]) > 0
        if local_bad.any():
            flagged += 1
            rejected[indices[local_bad]] = True
    return rejected, flagged


def _burst_mask(times, pixel_labels, group_labels, subunit_labels,
                pixel_coordinates, min_hits, min_channels, min_subunits, config):
    """Union every qualifying sliding window, including bin-edge crossings.

    Distinct-channel coordinate moments avoid letting one repeating pixel
    dominate the spatial test. Counts and moments update in linear time after
    sorting. A difference array marks union membership without a quadratic
    number of slice writes for continuously firing electronics. Eigenvalues
    are evaluated in bounded batches, keeping the original moment update
    order and covariance arithmetic without one LAPACK dispatch per hit.
    """
    rejected = np.zeros(len(times), dtype=bool)
    flagged_groups, maximum_channels = 0, 0
    pixel_outer = pixel_coordinates[:, :, None] * pixel_coordinates[:, None, :]
    batch_size = 1024
    sums = np.empty((batch_size, 2), dtype='f8')
    outers = np.empty((batch_size, 2, 2), dtype='f8')
    sizes = np.empty(batch_size, dtype='i8')
    starts = np.empty(batch_size, dtype='i8')
    stops = np.empty(batch_size, dtype='i8')

    def mark_windows(count, coverage):
        mean = sums[:count] / sizes[:count, None]
        covariance = (outers[:count] / sizes[:count, None, None]
                      - mean[:, :, None] * mean[:, None, :])
        passing = np.linalg.eigvalsh(covariance)[:, 0] >= config.burst_min_transverse_rms_cm ** 2
        # Starts are unique, but many nested windows can have the same stop.
        coverage[starts[:count][passing]] += 1
        np.add.at(coverage, stops[:count][passing], -1)

    for indices in _groups(group_labels):
        if len(indices) < min_hits:
            continue
        indices = indices[np.argsort(times[indices], kind='stable')]
        local_time = times[indices]
        pixels = pixel_labels[indices]
        subunits = subunit_labels[indices]
        pixel_counts, subunit_counts = {}, {}
        coord_sum, coord_outer = np.zeros(2), np.zeros((2, 2))
        coverage = np.zeros(len(indices) + 1, dtype=np.int64)
        right, pending = 0, 0
        for left in range(len(indices)):
            while right < len(indices) and local_time[right] - local_time[left] <= config.window_us:
                pixel, subunit = int(pixels[right]), int(subunits[right])
                if pixel not in pixel_counts:
                    coordinate = pixel_coordinates[pixel]
                    coord_sum += coordinate
                    coord_outer += pixel_outer[pixel]
                pixel_counts[pixel] = pixel_counts.get(pixel, 0) + 1
                subunit_counts[subunit] = subunit_counts.get(subunit, 0) + 1
                right += 1
            channels = len(pixel_counts)
            maximum_channels = max(maximum_channels, channels)
            if right - left >= min_hits and channels >= min_channels and len(subunit_counts) >= min_subunits:
                sums[pending] = coord_sum
                outers[pending] = coord_outer
                sizes[pending] = channels
                starts[pending], stops[pending] = left, right
                pending += 1
                if pending == batch_size:
                    mark_windows(pending, coverage)
                    pending = 0
            pixel, subunit = int(pixels[left]), int(subunits[left])
            pixel_counts[pixel] -= 1
            if pixel_counts[pixel] == 0:
                del pixel_counts[pixel]
                coordinate = pixel_coordinates[pixel]
                coord_sum -= coordinate
                coord_outer -= pixel_outer[pixel]
            subunit_counts[subunit] -= 1
            if subunit_counts[subunit] == 0:
                del subunit_counts[subunit]
        if pending:
            mark_windows(pending, coverage)
        local_bad = np.cumsum(coverage[:-1]) > 0
        if local_bad.any():
            flagged_groups += 1
            rejected[indices[local_bad]] = True
    return rejected, flagged_groups, maximum_channels


def quality_hit_mask(hits, config=None, *, tile_ids=None, module_ids=None):
    """Return an input-aligned keep mask and JSON-safe noise evidence.

    Prefer geometry-derived ``tile_ids`` and ``module_ids`` aligned with hits.
    Without them the explicit 2x2 mapping is tile=(io_channel-1)//4+1 and
    module=(io_group-1)//2+1. The tile mapping merges alternate UART routes to
    the same physical channel. It must not be assumed for another detector.
    Invalid required hit fields fail closed; isolated valid hits are retained.
    An ambiguous PPS interval vetoes only that io_group's hits. All coincidence
    scopes, including the module test, stay within the same io_group.
    """
    config = TrackQualityConfig() if config is None else config
    if not isinstance(config, TrackQualityConfig):
        raise TypeError('config must be TrackQualityConfig')
    required = ('x', 'y', 'z', 'ts_pps', 'io_group', 'io_channel', 'chip_id', 'channel_id')
    missing = set(required) - set(hits.dtype.names or ())
    if missing:
        raise ValueError('hit quality requires fields: ' + ', '.join(sorted(missing)))
    if hits.ndim != 1:
        raise ValueError('quality_hit_mask expects a one-dimensional event hit array')
    valid = np.ones(len(hits), dtype=bool)
    for field in required:
        valid &= np.isfinite(hits[field])
    valid &= (hits['io_group'] > 0) & (hits['io_channel'] > 0)
    valid &= (hits['chip_id'] >= 0) & (hits['channel_id'] >= 0)
    mapping = {}
    identifiers = []
    for supplied, fallback, name in (
        (tile_ids, (hits['io_channel'].astype('i8') - 1) // 4 + 1, 'tile'),
        (module_ids, (hits['io_group'].astype('i8') - 1) // 2 + 1, 'module'),
    ):
        values = np.asarray(fallback if supplied is None else supplied)
        if values.shape != (len(hits),):
            raise ValueError(f'{name}_ids must be aligned with the event hits')
        if np.any(~np.isfinite(values[valid])) or np.any(values[valid] < 0):
            raise ValueError(f'{name}_ids contains missing or invalid geometry')
        if np.any(values[valid] != values[valid].astype('i8')):
            raise ValueError(f'{name}_ids must contain integer geometry identities')
        identifiers.append(values[valid].astype('i8'))
        mapping[name] = 'caller-provided geometry' if supplied is not None else '2x2 four UARTs per tile' if name == 'tile' else '2x2 paired io_groups per module'
    # FLOW corrects and unrolls charge clocks per io_group. A long span of
    # phases across groups does not imply that each group's local timing is
    # ambiguous. Do not require one short interval for the whole event.
    candidate = hits[valid]
    clock_valid = np.ones(len(candidate), dtype=bool)
    times = np.zeros(len(candidate), dtype='f8')
    clocks, unwrapped = [], False
    for indices in _groups(candidate['io_group']):
        ticks = candidate['ts_pps'][indices]
        clock = dict(io_group=int(candidate['io_group'][indices[0]]),
                     n_hits=len(indices), phase_min_ticks=float(ticks.min()),
                     phase_max_ticks=float(ticks.max()))
        try:
            local_times, local_unwrapped = _arrival_times(ticks, config)
        except ArrivalTimeError as exc:
            clock_valid[indices] = False
            clock.update(status='ambiguous', reason=str(exc))
        else:
            times[indices] = local_times
            unwrapped |= local_unwrapped
            clock.update(status='ok', span_us=float(np.ptp(local_times)),
                         timestamps_unwrapped=local_unwrapped)
        clocks.append(clock)
    eligible = valid.copy()
    eligible[valid] = clock_valid
    tiles, modules = (values[clock_valid] for values in identifiers)
    selected = hits[eligible]
    times = times[clock_valid]
    xyz = np.column_stack([selected[field] for field in ('x', 'y', 'z')])
    groups = selected['io_group'].astype('i8')
    chips = selected['chip_id'].astype('i8')
    pixels = _labels([groups, tiles, chips, selected['channel_id']])
    chip_groups = _labels([groups, tiles, chips])
    tile_groups = _labels([groups, tiles])
    # Independent recentering cannot establish cross-io_group coincidence.
    # Evaluate module breadth within each clock domain, including its tiles,
    # rather than inventing simultaneity between unrelated PPS phases.
    module_groups = _labels([modules, groups])
    transverse = [i for i, axis in enumerate(('x', 'y', 'z')) if axis != config.drift_axis]
    # Center the pixel-plane coordinates to avoid cancellation in translated
    # detector coordinates (e.g. z ~1300 cm).
    _, representatives = np.unique(pixels, return_index=True)
    coordinates = xyz[representatives][:, transverse]
    if len(coordinates):
        coordinates = coordinates - coordinates.mean(axis=0)
    hot, hot_pixels = _hot_mask(times, xyz, pixels, config)
    burst = np.zeros(len(selected), dtype=bool)
    counts, maxima = {}, {}
    for name, grouped, subunits, min_hits, min_channels, min_subunits in (
        ('chip', chip_groups, chip_groups, config.chip_min_hits, config.chip_min_channels, 1),
        ('tile', tile_groups, chip_groups, config.tile_min_hits, config.tile_min_channels, config.tile_min_chips),
        ('module', module_groups, tile_groups, config.module_min_hits, config.module_min_channels, config.module_min_tiles),
    ):
        mask, count, maximum = _burst_mask(times, pixels, grouped, subunits,
            coordinates, min_hits, min_channels, min_subunits, config)
        burst |= mask
        counts[name] = count
        maxima[name] = maximum
    keep = eligible.copy()
    keep[eligible] = ~(hot | burst)
    diagnostics = dict(
        version=QUALITY_VERSION, thresholds_status='provisional; not efficiency-calibrated',
        config=config.to_dict(), input_hits=len(hits), kept_hits=int(keep.sum()),
        rejected_hits=int((~keep).sum()), invalid_hits=int((~valid).sum()),
        ambiguous_clock_hits=int((~clock_valid).sum()),
        ambiguous_clock_groups=sum(clock['status'] == 'ambiguous' for clock in clocks),
        clock_domains=clocks,
        hot_pixel_hits=int(hot.sum()), coherent_noise_hits=int(burst.sum()),
        hot_pixels_flagged=hot_pixels, burst_groups_flagged=counts,
        maximum_distinct_channels= maxima, geometry_mapping=mapping,
        timestamp_field='ts_pps', timestamps_unwrapped=unwrapped,
        coincidence_scope='within io_group; module breadth evaluated separately for each io_group',
        scope='one event; no cross-event noise history',
    )
    return keep, diagnostics
