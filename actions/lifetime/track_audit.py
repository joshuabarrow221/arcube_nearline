#!/usr/bin/env python3
"""Explain every track-window mean using its retained direct-FLOW fits.

This is a diagnostic reader, not a refitter or an outlier-removal algorithm.
Leave-one-file-out means quantify sensitivity to a particular file; they do
not establish whether that file or the remaining files measure true purity.
"""
import argparse
import json
from pathlib import Path

import numpy as np

from lifetime import aggregate_track_lifetimes, is_flow_track, record_key, valid_lifetime
from purity_sources import atomic_json


def audit_tracks(entries, timezone_name='America/Chicago'):
    """Report unreviewed and reviewed summaries with per-file influence.

    Construct window membership with the same calendar rules as publication.
    The unreviewed view clears only the monitoring exclusion flag; failed and
    reference-only rows still cannot participate. Thus it answers what changed
    because of review, without resurrecting unusable fits.
    """
    originals = {record_key(r): r for r in entries if is_flow_track(r)}
    unreviewed = [dict(r, monitoring_excluded=False) for r in originals.values()]
    reviewed = {(r['sample'], r['period_start']): r
                for r in aggregate_track_lifetimes(list(originals.values()), timezone_name)}
    windows = []
    for summary in aggregate_track_lifetimes(unreviewed, timezone_name):
        members = [originals[('track', m['input_file'], summary['sample'])] for m in summary['members']]
        valid = [r for r in members if valid_lifetime(r)]
        values = np.array([r['lifetime_us'] for r in valid], dtype=float)
        nominal = reviewed[(summary['sample'], summary['period_start'])]
        contributions = []
        for i, row in enumerate(valid):
            others = np.delete(values, i)
            reduced = float(others.mean()) if len(others) else None
            contributions.append(dict(input_file=row['input_file'], timestamp=row['timestamp'],
                lifetime_us=row['lifetime_us'], error_us=row.get('error_us'),
                n_segments=row.get('n_segments'), run_context=row.get('run_context'),
                chi2_ndf=row.get('fit_diagnostics', {}).get('chi2_ndf'),
                excluded_from_monitoring=row.get('monitoring_excluded', False),
                review_reason=row.get('monitoring_review_reason'),
                unreviewed_weight=1/len(valid), leave_one_out_mean_us=reduced,
                mean_change_when_omitted_us=reduced-float(values.mean()) if reduced is not None else None))
        windows.append(dict(sample=summary['sample'], period_start=summary['period_start'],
            period_end=summary['period_end'], n_valid_files=len(valid),
            n_failed_files=len(members)-len(valid), unreviewed_mean_us=summary['lifetime_us'],
            unreviewed_error_us=summary['error_us'],
            minimum_us=float(values.min()) if len(values) else None,
            maximum_us=float(values.max()) if len(values) else None,
            median_us=float(np.median(values)) if len(values) else None,
            monitoring_mean_us=nominal['lifetime_us'], monitoring_error_us=nominal['error_us'],
            n_excluded=nominal['n_excluded'], contributions=contributions))
    return dict(timezone=timezone_name, windows=windows,
                interpretation='Equal file weights; leave-one-out sensitivity is descriptive, not an exclusion criterion.')


def main():
    """Write an inspectable strict-JSON audit without changing the input history."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('history')
    parser.add_argument('--output', required=True)
    parser.add_argument('--timezone', default='America/Chicago')
    args = parser.parse_args()
    history = json.loads(Path(args.history).read_text())
    atomic_json(args.output, audit_tracks(history['lifetimes'], args.timezone))


if __name__ == '__main__':
    main()
