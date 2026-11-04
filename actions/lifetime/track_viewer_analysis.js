/* Saved quality evidence and full-pool charge-scale diagnostics. No browser fits. */
(() => {
  'use strict';
  const finite = value => typeof value === 'number' && Number.isFinite(value);
  const fmt = (value, digits=3) => finite(value) ? value.toLocaleString('en-US', {maximumFractionDigits:digits}) : 'not recorded';
  const el = (tag, text, className) => { const node=document.createElement(tag); if(text!==undefined)node.textContent=text; if(className)node.className=className; return node; };
  const style=el('style');
  style.textContent='#quality-explanation{padding:14px 16px;margin:14px 0;background:#fff;border:1px solid #d7e0e9;border-radius:8px}#quality-explanation p+p{margin-top:7px}#quality-explanation table{border-collapse:collapse;margin-top:9px;font-size:12px;width:100%}#quality-explanation th,#quality-explanation td{padding:6px 10px;text-align:left;border-bottom:1px solid #e8edf2}#charge-diagnostics{margin-top:16px}#charge-diagnostics .charge-controls{display:flex;gap:8px;flex-wrap:wrap;align-items:center}#charge-diagnostics .charge-facts{padding:12px 16px;display:flex;gap:10px 25px;flex-wrap:wrap}#charge-diagnostics .charge-facts strong{display:block}#charge-plot{height:370px;width:100%}#charge-model-note{padding:0 16px 12px;color:#596c7e;font-size:12px}';
  document.head.append(style);
  const quality=el('section',undefined,''); quality.id='quality-explanation';
  document.querySelector('.controls').before(quality);
  const section=el('section',undefined,'plot-panel');section.id='charge-diagnostics';
  const head=el('div',undefined,'plot-head');head.append(el('h2','Full-pool charge scale at zero drift'));
  const controls=el('div',undefined,'charge-controls');
  const buttons={};
  for(const [key,title] of [['direct','Direct dQ/ds'],['proxy','Saved lifetime proxy']]) {
    const button=el('button',title);button.type='button';button.setAttribute('aria-pressed',String(key==='direct'));
    button.addEventListener('click',()=>{mode=key;renderCharge();});buttons[key]=button;controls.append(button);
  }
  head.append(controls);
  const facts=el('div',undefined,'charge-facts');facts.setAttribute('aria-live','polite');
  const note=el('p');note.id='charge-model-note';
  const plot=el('div');plot.id='charge-plot';plot.setAttribute('aria-label','Full-pool drift-bin MPVs and exponential extrapolation to zero drift');
  const foot=el('div','These fits use the full pooled input, not the 10–20 displayed tracks. The dashed interval down to t = 0 is an extrapolation. Changing the track page, focus, or rotation does not refit this curve.','plot-foot');
  section.append(head,facts,note,plot,foot);
  document.querySelector('.profile-panel').after(section);
  let data=null,mode='direct';

  function renderQuality() {
    quality.replaceChildren();
    quality.append(el('h2','What “include rejected context” means'));
    quality.append(el('p','This display switch reveals exported hits and pieces of rejected segments belonging to the exported parent tracks. Hits may be subsampled for display. It does not run cuts or change the pool. One failed member hit vetoes its entire segment; a rejected marker is not a truth label for instrumental noise.','small'));
    const evidence=data.quality_evidence;
    if(!evidence?.available) {
      quality.append(el('p',data.summary?.selection==='stored' ? 'Stored mode applies no new hot-pixel or coincidence veto. The upstream selection and numerical validity checks still apply.' : 'This export has no saved aggregate event evidence. Re-export with the current viewer to include it.','muted small'));
      return;
    }
    const counts=evidence.counters || {}, groups=evidence.burst_groups_flagged || {};
    quality.append(el('p',`${fmt(evidence.n_quality_rejected,0)} quality-vetoed segments out of ${fmt(evidence.n_source_segments,0)} original segments. Recorded event flags: ${fmt(counts.hot_pixel_hits,0)} hot-pixel hits; ${fmt(counts.coherent_noise_hits,0)} coherent-window hits. Qualifying chip/tile/module group counts: ${fmt(groups.chip,0)} / ${fmt(groups.tile,0)} / ${fmt(groups.module,0)}.`));
    quality.append(el('p',`Evidence covers ${fmt(evidence.n_events,0)} unique parent events across the receipt, including hits outside the displayed tracks. Policy: ${(evidence.versions || []).join(', ')}. Chip/tile/module flags can overlap; do not add them as disjoint populations.`, 'muted small'));
    const details=el('details');details.append(el('summary','Recorded thresholds and other vetoes'));
    let config=data.summary?.quality?.quality_config;
    try{if(typeof config==='string')config=JSON.parse(config);}catch{config=null;}
    if(config) {
      details.append(el('p',`Coincidence uses sliding ${fmt(config.window_us)} µs arrival-time windows, not just how horizontal a reconstructed object looks. The distinct-pixel transverse minor RMS must also be at least ${fmt(config.burst_min_transverse_rms_cm)} cm. These thresholds are provisional and have no measured signal-efficiency calibration.`,'small'));
      const table=el('table'),header=el('tr');for(const title of ['Scope','Minimum hits','Distinct channels','Additional condition'])header.append(el('th',title));table.append(header);
      for(const [title,hits,channels,extra] of [
        ['Chip',config.chip_min_hits,config.chip_min_channels,'Transverse spread guard'],
        ['Tile',config.tile_min_hits,config.tile_min_channels,`${fmt(config.tile_min_chips,0)} chips`],
        ['Module/readout scope',config.module_min_hits,config.module_min_channels,`${fmt(config.module_min_tiles,0)} tiles`]
      ]){const row=el('tr');[title,fmt(hits,0),fmt(channels,0),extra].forEach(value=>row.append(el('td',value)));table.append(row);}details.append(table);
      details.append(el('p',`Hot pixels are a different test: at least ${fmt(config.hot_min_hits,0)} hits on one physical pixel, rate ≥ ${fmt(config.hot_min_rate_hz/1000)} kHz, drift span ≥ ${fmt(config.hot_min_drift_span_cm)} cm, transverse span ≤ ${fmt(config.hot_max_transverse_span_cm)} cm, and drift-axis cosine ≥ ${fmt(config.hot_min_axis_cosine,4)}.`,'small'));
    }
    details.append(el('p',`Other recorded event counts: ambiguous-clock hits ${fmt(counts.ambiguous_clock_hits,0)}; invalid hit fields ${fmt(counts.invalid_hits,0)}; outside-module hits ${fmt(counts.n_outside_module_hits,0)}; unmapped hits ${fmt(counts.n_unmapped_hits,0)}. Missing counters mean not recorded, not zero.`,'small'));
    details.append(el('p','Saved event counters cannot identify the exact trigger for an individual segment. A real track or shower can overlap a flagged hardware/time window, or itself satisfy the occupancy test. The appearance alone cannot establish whether the veto was physically correct.','small'));
    quality.append(details);
  }

  async function renderCharge() {
    if(!data)return;
    for(const [key,button] of Object.entries(buttons))button.setAttribute('aria-pressed',String(key===mode));
    const result=data.charge_scale?.[mode];facts.replaceChildren();
    function fact(title,value){const box=el('div');box.append(el('span',title,'muted small'),el('strong',value));facts.append(box);}
    const accepted=result?.status==='ok';
    fact('Fit status',result?.status || 'not exported');
    fact(mode==='direct' ? 'dQ/ds at t = 0' : 'Proxy at t = 0',accepted ? `${fmt(result.intercept)} ± ${fmt(result.intercept_error)}` : 'Unavailable');
    fact('Lifetime in this model',accepted ? `${fmt(result.tau_us)} ± ${fmt(result.tau_error_us)} µs` : 'Unavailable');
    fact('χ² / ndf',finite(result?.chi2_ndf) ? fmt(result.chi2_ndf) : 'Unavailable');
    fact('Intercept–lifetime correlation',finite(result?.correlation) ? fmt(result.correlation) : 'Unavailable');
    const description=mode==='direct' ? 'Independent diagnostic: Landau–Gaussian MPVs of full-pool dQ/ds, then A exp(−t/τ). Here ds is the stored segment path length dx; it is not the detector X coordinate. This does not replace the saved lifetime estimator.' : 'Existing lifetime estimator: MPV of a legacy histogram-convolution proxy, then A exp(−t/τ). Its intercept is in proxy units, not raw charge per path length.';
    note.textContent=description+' '+(result?.fit_message || (result ? (result.uncertainty || 'Uncertainties are conditional fit-covariance errors; detector/selection systematics and within-track correlations are not included.') : 'Re-export this sample with the current track viewer to include charge-scale diagnostics.'));
    const traces=[];
    const points=result?.fitpoints;
    if(points?.x?.length)traces.push({type:'scatter',mode:'markers',x:points.x,y:points.y,error_y:{type:'data',array:points.error,visible:true},marker:{size:8,color:'#146b80'},name:'Full-pool drift-bin MPVs',hovertemplate:'t=%{x:.2f} µs<br>MPV=%{y:.4g}<extra></extra>'});
    const curve=result?.curve,edge=result?.fit_range_us?.[0];
    if(accepted && curve?.x?.length) {
      if(curve.error?.length===curve.x.length) {
        traces.push({type:'scatter',mode:'lines',x:curve.x,y:curve.y.map((v,i)=>v+curve.error[i]),line:{width:0},hoverinfo:'skip',showlegend:false});
        traces.push({type:'scatter',mode:'lines',x:curve.x,y:curve.y.map((v,i)=>v-curve.error[i]),line:{width:0},fill:'tonexty',fillcolor:'rgba(20,107,128,.12)',hoverinfo:'skip',name:'Conditional 1σ mean curve'});
      }
      for(const extrapolate of [true,false]) {
        const keep=curve.x.map((value,index)=>({value,index})).filter(({value})=>extrapolate ? value<=edge : value>=edge);
        traces.push({type:'scatter',mode:'lines',x:keep.map(v=>v.value),y:keep.map(v=>curve.y[v.index]),line:{color:'#d57522',width:2,dash:extrapolate?'dash':'solid'},name:extrapolate?'Extrapolation to zero drift':'Exponential fit'});
      }
      traces.push({type:'scatter',mode:'markers',x:[0],y:[result.intercept],error_y:{type:'data',array:[result.intercept_error],visible:true},marker:{symbol:'diamond',size:10,color:'#d57522'},name:'Zero-drift intercept'});
    }
    const safe=value=>String(value||'unknown').replace(/[^A-Za-z0-9_+\-]/g,'');
    await Plotly.react(plot,traces,{margin:{l:75,r:25,t:15,b:55},paper_bgcolor:'#fff',plot_bgcolor:'#fff',font:{family:'system-ui, sans-serif',size:12,color:'#354e63'},xaxis:{title:{text:'Drift time [µs]'},rangemode:'tozero',gridcolor:'#e5ebf1'},yaxis:{title:{text:mode==='direct'?'dQ/ds MPV [FLOW charge units / cm]':'Legacy proxy MPV [proxy units]'},gridcolor:'#e5ebf1'},legend:{orientation:'h',y:-.25},uirevision:mode+'-'+data.summary?.receipt,annotations:traces.length?[]:[{text:'No accepted fit or usable MPV points available',xref:'paper',yref:'paper',x:.5,y:.5,showarrow:false}]},{responsive:true,displaylogo:false,toImageButtonOptions:{format:'png',scale:2,filename:'charge_scale__'+safe(data.summary?.period_start)+'_to_'+safe(data.summary?.period_end)+'__'+safe(data.summary?.selection)+'__'+mode}});
  }
  function load(value){data=value;renderQuality();renderCharge().catch(error=>{note.textContent='Charge display failed: '+error.message;});}
  document.addEventListener('trackviewer:load',event=>load(event.detail));
  load(JSON.parse(document.getElementById('viewer-data').textContent));
})();
