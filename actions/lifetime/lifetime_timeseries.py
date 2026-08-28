#!/usr/bin/env python3
import argparse
import json
from collections import defaultdict
from datetime import datetime

import matplotlib
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import plotly.graph_objects as go


def load_entries(input_file):
    with open(input_file) as stream:
        raw_entries = json.load(stream)['lifetimes']

    entries_by_sample = defaultdict(list)
    for entry in raw_entries:
        if entry.get('fit_status', 'ok') != 'ok':
            continue
        if entry.get('lifetime_us') is None or entry.get('error_us') is None:
            continue
        sample = entry.get('sample', 'mixed')
        entries_by_sample[sample].append(
            (
                datetime.fromisoformat(entry['timestamp']),
                entry['lifetime_us'],
                entry['error_us'],
            )
        )

    for entries in entries_by_sample.values():
        entries.sort(key=lambda item: item[0])
    return dict(entries_by_sample)


def draw_matplotlib(entries_by_sample, output_file, last_n=None):
    fig, ax = plt.subplots(figsize=(10, 5))
    for sample, entries in sorted(entries_by_sample.items()):
        selected = entries if last_n is None else entries[-last_n:]
        if not selected:
            continue
        timestamps, lifetimes, errors = zip(*selected)
        ax.errorbar(
            timestamps,
            lifetimes,
            yerr=errors,
            fmt='o-',
            capsize=1,
            label=sample,
        )

    ax.set_xlabel('Timestamp [America/Chicago]')
    ax.set_ylabel('Electron lifetime (µs)')
    ax.set_ylim(0, 10_000)
    title = 'Electron lifetime time series'
    if last_n is not None:
        title += f' (last {last_n} points per sample)'
    ax.set_title(title)
    ax.grid(True)
    ax.legend(title='Muon sample')
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%m/%d/%Y\n%H:%M %Z'))
    fig.tight_layout()
    fig.savefig(output_file)
    plt.close(fig)


def draw_plotly(entries_by_sample, output_file):
    fig = go.Figure()
    for sample, entries in sorted(entries_by_sample.items()):
        timestamps, lifetimes, errors = zip(*entries)
        fig.add_trace(
            go.Scatter(
                x=timestamps,
                y=lifetimes,
                error_y={'type': 'data', 'array': errors, 'visible': True},
                mode='markers',
                name=sample,
            )
        )

    fig.update_layout(
        title='Electron lifetime time series',
        xaxis_title='Timestamp [America/Chicago]',
        yaxis_title='Electron lifetime (µs)',
        yaxis_range=[0, 10_000],
        template='plotly_white',
        width=1000,
        height=500,
        legend_title='Muon sample',
    )
    fig.update_xaxes(tickformat='%m/%d/%Y<br>%H:%M %Z')
    # Keep the historical URL: output_file is the .png path and .html is appended.
    fig.write_html(output_file + '.html')


def main(input_file, output_file):
    entries_by_sample = load_entries(input_file)
    if not entries_by_sample:
        raise ValueError('no successful lifetime entries to plot')

    matplotlib.rcParams['timezone'] = 'America/Chicago'
    draw_matplotlib(entries_by_sample, output_file)
    draw_matplotlib(entries_by_sample, output_file + '_last.png', last_n=50)
    draw_plotly(entries_by_sample, output_file)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--input_file',
        required=True,
        help='Path to the electron-lifetime JSON file',
    )
    parser.add_argument(
        '--output_file',
        required=True,
        help='Output PNG path; an HTML plot is written by appending .html',
    )
    args = parser.parse_args()
    print(args)
    main(**vars(args))
