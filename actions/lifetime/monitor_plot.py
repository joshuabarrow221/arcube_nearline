"""Compact shifter plot and a separate complete fit-status page.

The public chart preserves the existing PNG/.png.html URL convention and a
simple Plotly-style presentation. All tracks includes the cosmic-enriched
subset, so their two curves are correlated. Both use the same candidate cuts;
cosmic-enriched adds the event condition n_ext_trigs == 0. The externally
triggered subset remains in diagnostics only. Gray coverage boxes,
long annotations and per-file paths do not obscure the measurement points.
"""
from collections import defaultdict
from html import escape
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from zoneinfo import ZoneInfo

from lifetime import valid_lifetime, packet_candidate_points, TRACK_SAMPLE_DESCRIPTIONS
from purity_sources import aware_time

# Fixed colors remain consistent between PNG, interactive traces and legends.
STYLES = {'prm': ('PRM · daily', '#2878B5', 'o', 'circle'),
          'all_mip': ('All tracks · 6 h', '#EF553B', 'D', 'diamond'),
          'cosmic': ('Cosmic-enriched tracks · 6 h', '#AB63FA', 'v', 'triangle-down'),
          'packet': ('Packets · 6 h', '#19A7A7', 'P', 'cross'),
          'gas': ('Gas O₂ + H₂O · 6 h', '#819B28', 's', 'square'),
          'gas_o2': ('Gas O₂ only · 6 h', '#819B28', 'x', 'x')}


def display_value(row):
    """A candidate remains failed; display copies never promote fit status."""
    if valid_lifetime(row): return row['lifetime_us'], row.get('error_us'), False
    candidates = packet_candidate_points([dict(row, monitoring_excluded=False)])
    if candidates:
        r = candidates[0]; return r['display_lifetime_us'], r['display_error_us'], True
    return None, None, False


def observed_time(row):
    """Use mean actual input times where known; never imply a full-day exposure."""
    if row.get('members'):
        seconds = [aware_time(m['timestamp']).timestamp() for m in row['members']]
        return pd.Timestamp(float(np.mean(seconds)), unit='s', tz='UTC').round('us')
    return aware_time(row['timestamp'])


def render(payload, directory, config):
    directory = Path(directory); tz = ZoneInfo(config['timezone'])
    rows = payload['lifetimes']
    raw_prm = payload.get('raw_prm_measurements', [])
    # Locate daily PRM means at their actual measurement-time centroid. Keep
    # the original aggregation window in JSON and hover for precise provenance.
    positioned = []
    for original in rows:
        r = dict(original)
        if r.get('sample') == 'prm':
            raw = [p for p in raw_prm if p['source'] == r['source'] and
                   aware_time(r['period_start']) <= aware_time(p['timestamp']) < aware_time(r['period_end'])]
            if raw:
                r['timestamp'] = pd.Timestamp(np.mean([aware_time(p['timestamp']).timestamp() for p in raw]), unit='s', tz='UTC').round('us').isoformat()
        positioned.append(r)
    newest = max((aware_time(r['timestamp']) for r in positioned), default=aware_time(payload['generated_at']))
    cutoff = newest-pd.Timedelta(days=config['display_days'])
    groups = defaultdict(list)
    omitted = []
    for r in positioned:
        sample = r.get('sample')
        if sample not in STYLES or r.get('monitoring_excluded'): continue
        if sample in ('gas', 'gas_o2') and sample not in config.get('gas_samples', ['gas']): continue
        if aware_time(r['timestamp']) < cutoff: continue
        v,e,candidate = display_value(r)
        if v is None: continue
        if v > config.get('maximum_display_us', 10000):
            omitted.append(r); continue
        # Connect the displayed means of a method continuously. This is only
        # presentation: independent cohorts/configurations were fitted in
        # separate pools and their identities remain in hover/diagnostics.
        # PRM and gas instruments must never be connected to another sensor.
        key = (sample, r.get('source', '') if sample in ('prm', 'gas', 'gas_o2') else '')
        groups[key].append((r,v,e,candidate))
    fig, ax = plt.subplots(figsize=(15,10), dpi=200)
    fig.subplots_adjust(left=.09,right=.965,bottom=.12,top=.85)
    interactive = go.Figure(); legend_seen = set()
    for (sample,source), points in sorted(groups.items()):
        points.sort(key=lambda p:observed_time(p[0]))
        name,color,marker,symbol = STYLES[sample]
        if sample in ('gas','gas_o2'):
            tag = '1874' if '1874' in source else '1890' if '1890' in source else source
            name += f' ({tag})'
            if '1874' in source: color = '#8D73AF'
        showlegend = name not in legend_seen; legend_seen.add(name)
        times = [observed_time(r).to_pydatetime().astimezone(tz) for r,v,e,c in points]
        values = [v for r,v,e,c in points]
        # Lines connect available measurements, not fabricated/interpolated data.
        ax.plot(times,values,color=color,lw=1.2,alpha=.85)
        for time,(r,v,e,candidate) in zip(times,points):
            ax.errorbar(time,v,yerr=e,fmt=marker,ms=6.5,color=color,capsize=2,lw=1,
                        mfc='none' if candidate else color, label=name if showlegend else '_nolegend_')
            showlegend=False
        hover=[]
        for r,v,e,c in points:
            hover.append(escape(str(r.get('source','')))+'<br>'+escape(str(r.get('period_start','')))+' — '+escape(str(r.get('period_end','')))+
                '<br>Files: '+str(r.get('n_files','—'))+'; observations: '+escape(str(r.get('n_segments',r.get('n_objects',r.get('n_measurements','—')))))+
                '<br>'+escape(str(r.get('uncertainty','')))+'<br>'+escape(str(r.get('fit_message','')))+
                ('<br>'+escape(TRACK_SAMPLE_DESCRIPTIONS[sample]) if sample in TRACK_SAMPLE_DESCRIPTIONS else ''))
        interactive.add_trace(go.Scatter(x=[t.isoformat() for t in times],y=values,mode='lines+markers',
            name=name,legendgroup=name,showlegend=not any(t.name==name for t in interactive.data),
            line=dict(color=color,width=1.5),marker=dict(color=color,size=8,symbol=[symbol+'-open' if c and symbol!='x' else symbol for r,v,e,c in points]),
            error_y=dict(type='data',array=[e for r,v,e,c in points],visible=True),text=hover,
            hovertemplate='%{x}<br>%{y:.1f} µs<br>%{text}<extra>%{fullData.name}</extra>'))
    # Place the title in figure coordinates, above the two-row legend; an axes
    # title with padding would collide with that legend on the exported PNG.
    fig.suptitle('2×2 Electron Lifetime',fontsize=22,y=.975)
    ax.set_ylabel('Electron lifetime [µs]',fontsize=16); ax.set_xlabel('Date · '+config['timezone'],fontsize=14)
    ax.tick_params(labelsize=12);ax.grid(color='#EBF0F8',lw=.8);ax.set_axisbelow(True)
    # Use the available vertical space while keeping the fixed diagnostic cap.
    # Error bars participate in scaling; this never changes fit acceptance or
    # removes a point merely to make the chart look more consistent.
    upper = [v+(e or 0) for points in groups.values() for r,v,e,c in points]
    y_top = min(config['maximum_display_us'], max(1000, 1000*np.ceil(max(upper)*1.08/1000))) if upper else config['maximum_display_us']
    ax.set_ylim(bottom=0,top=y_top)
    if not groups:
        ax.text(.5,.5,'No usable lifetime measurements yet',transform=ax.transAxes,ha='center',fontsize=16)
        ax.set_xlim(cutoff.to_pydatetime(), newest.to_pydatetime())
    ax.xaxis.set_major_locator(mdates.AutoDateLocator(tz=tz));ax.xaxis.set_major_formatter(mdates.ConciseDateFormatter(ax.xaxis.get_major_locator(),tz=tz))
    if groups:ax.legend(loc='lower center',bbox_to_anchor=(.5,1.015),ncol=3,frameon=False,fontsize=11)
    # One short caveat preserves the candidate meaning without the previous
    # shaded boxes or a large in-plot status narrative.
    fig.text(.09,.04,'Open packet markers: checks failed · Error definitions and coverage in diagnostics',fontsize=10,color='#555555')
    fig.savefig(directory/'elifetime_time_series.png',dpi=200);plt.close(fig)
    interactive.update_layout(template='plotly_white',title='2×2 Electron Lifetime',height=700,
        xaxis_title='Date · '+config['timezone'],yaxis_title='Electron lifetime [µs]',
        yaxis=dict(range=[0,y_top]),
        legend=dict(orientation='h',y=1.14,x=0),margin=dict(t=140),hovermode='closest')
    # Self-contained Plotly avoids a CDN dependency on the shifter website.
    chart=interactive.to_html(full_html=False,include_plotlyjs=True,config={'responsive':True,'displaylogo':False})
    page='<html><head><meta charset="utf-8"><title>2×2 Electron Lifetime</title></head><body style="font-family:Arial;margin:20px">'+chart
    page+='<p>PRM: daily means. Tracks and packets: 6 h pooled fits. Gases: 6 h equivalents. Open packet markers fail checks. '
    page+='<a href="elifetime_diagnostics.html">Fit diagnostics and coverage</a> · <a href="elifetime_time_series.png">PNG</a></p></body></html>'
    (directory/'elifetime_time_series.png.html').write_text(page)
    # Full failures and redundant external-trigger fits stay available off the
    # money plot. Escape every source string; file metadata is untrusted text.
    body='<h1>Lifetime diagnostics</h1><p>Generated '+escape(payload['generated_at'])+'</p>'
    body+='<p>Means/pools do not imply continuous exposure. Trigger timing does not prove beam origin. Gas conversion and packet selection remain provisional.</p>'
    # Put the full definition beside diagnostic results, keeping the main PNG
    # uncluttered. These descriptions explain existing selections; they never
    # change pool membership, fit acceptance, numerical values or error bars.
    body+='<h2>Track sample definitions</h2><dl>'
    for sample,label in [('all_mip','All tracks'),('cosmic','Cosmic-enriched tracks'),('beam','Externally triggered tracks (diagnostics only)')]:
        body+='<dt><strong>'+label+'</strong></dt><dd>'+escape(TRACK_SAMPLE_DESCRIPTIONS[sample])+'</dd>'
    body+='</dl><p>All tracks is fitted to the union of the two timing categories, not an average of their fitted lifetimes. '
    body+='The cosmic-enriched and All tracks estimates share segments and are correlated; agreement is not independent corroboration. '
    body+='Different membership can change drift coverage, statistical precision and fitted lifetime. '
    body+='If every selected event has zero external triggers, the two inputs coincide. If none does, the cosmic-enriched fit is unavailable, not zero. '
    body+='Legacy inputs without timing associations support only All tracks for the selection actually stored upstream.</p>'
    body+='<p>'+escape(payload['uncertainty'])+'</p><p><a href="elifetime_status.json">Full provenance/status JSON</a> · <a href="elifetime_time_series.png.html">Main plot</a></p>'
    if payload['input_errors']: body+='<h2>Input processing errors</h2><pre>'+escape(json.dumps(payload['input_errors'],indent=2))+'</pre>'
    if omitted:body+='<p>'+str(len(omitted))+' estimates exceed the configured display range; their values are retained below.</p>'
    body+='<table><tr><th>Period</th><th>Method</th><th>Source/cohort</th><th>Files</th><th>Lifetime [µs]</th><th>Error [µs]</th><th>Status</th><th>χ²/ndf</th><th>Reason</th></tr>'
    for r in sorted(rows,key=lambda r:(r.get('timestamp',''),r.get('sample',''))):
        v,e,c=display_value(r)
        values=[r.get('period_start',r.get('timestamp')),r.get('sample'),r.get('source',r.get('cohort','')),
                r.get('n_files',''),v,e,('diagnostic cohort; ' if r.get('monitoring_excluded') else '')+r.get('fit_status',''),
                r.get('attenuation_chi2_ndf',r.get('fit_diagnostics',{}).get('chi2_ndf','')),r.get('fit_message','')]
        body+='<tr>'+''.join('<td>'+escape('' if x is None else str(x))+'</td>' for x in values)+'</tr>'
    body+='</table><h2>Input availability</h2><pre>'+escape(json.dumps(payload['inputs'],indent=2))+'</pre>'
    (directory/'elifetime_diagnostics.html').write_text('<html><head><meta charset="utf-8"><style>body{font-family:Arial;margin:24px}table{border-collapse:collapse}td,th{padding:8px;border:1px solid #ddd;text-align:left}pre{white-space:pre-wrap}</style></head><body>'+body+'</body></html>')
