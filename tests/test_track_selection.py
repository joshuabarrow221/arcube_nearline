"""Geometric regression checks, not a validation of detector cut efficiency."""
from pathlib import Path
from types import SimpleNamespace
import sys

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'actions/lifetime'))
from track_selection import projected_segments


HITS = np.dtype([(name, 'f8') for name in ('x','y','z','Q','E','t_drift')]
                + [('id','i8'),('io_group','i4')])
SEGMENTS = np.dtype([('rock_segment_id','i4')]
    + [(axis+'_'+suffix, 'f8') for axis in 'xyz' for suffix in ('start','end','mid')]
    + [(name,'f8') for name in ('dE','dQ','dx','t')] + [('nhits','i4'),('io_group','i4')])


def selector():
    return SimpleNamespace(segment_count=-1, rock_muon_segments_dtype=SEGMENTS)


@pytest.mark.parametrize('axis', ['x','y','z'])
def test_all_orientations_conserve_charge_and_unique_hit_membership(axis):
    hits = np.zeros(241, dtype=HITS)
    hits[axis] = np.linspace(-60, 60, len(hits))
    hits['Q'], hits['E'], hits['t_drift'] = 3., 1., np.arange(len(hits))
    hits['id'], hits['io_group'] = np.arange(len(hits)), 1
    rows, ids = projected_segments(selector(), hits)
    assert len(rows) == 60
    assert sum(row['nhits'] for row in rows) == len(hits)
    assert sum(row['dQ'] for row in rows) == hits['Q'].sum()
    np.testing.assert_array_equal(np.sort(np.concatenate(ids)), hits['id'])
    assert all(row['dx'] == pytest.approx(2.) for row in rows)
    assert min(row[axis+'_start'] for row in rows) >= -60
    assert max(row[axis+'_end'] for row in rows) <= 60


def test_disconnected_modules_are_segmented_locally_without_bridging():
    hits = np.zeros(42, dtype=HITS)
    hits['id'] = np.arange(len(hits))
    hits['io_group'] = np.repeat([1,3], 21)
    hits['z'] = np.tile(np.linspace(0,10,21), 2)
    hits['x'][21:] = 90
    hits['Q'] = 1
    rows, members = projected_segments(selector(), hits)
    assert len(rows) == 10
    for row, ids in zip(rows, members):
        assert np.unique(hits['io_group'][ids]).tolist() == [row['io_group']]
        assert row['x_start'] == row['x_end']
    assert sum(row['dQ'] for row in rows) == 42


def test_short_tail_has_actual_length_and_no_empty_zero_length_segment():
    hits = np.zeros(8,dtype=HITS)
    hits['z'], hits['id'], hits['io_group'] = np.linspace(0,3.5,8), np.arange(8), 1
    rows, members = projected_segments(selector(), hits)
    assert [row['dx'] for row in rows] == [2., 1.5]
    assert len(np.concatenate(members)) == 8
