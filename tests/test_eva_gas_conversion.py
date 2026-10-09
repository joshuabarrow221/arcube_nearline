"""Synthetic checks for units, source isolation, and refresh reproducibility."""
import json
from pathlib import Path
import sys
import numpy as np
import pandas as pd
import pytest
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'actions/lifetime'))
from purity_sources import aggregate_measurements, gas_conversion_rates
from lifetime import resolve_gas_conversion_config, update_json
CONFIG = {'model':'eva_20261006','field_strength_vpercm':212.125}


def test_field_rates_match_upstream_without_loading_workbook(monkeypatch):
    import plot_lifetime_history as eva
    monkeypatch.setattr(eva, 'load_data', lambda *a: pytest.fail('Workbook must not supply data'))
    k_o2,k_water,meta = gas_conversion_rates(CONFIG)
    for oxygen,water in [(1.26,32.),(15.,48.),(0.,2.)]:
        assert 1/(k_o2*oxygen+k_water*water) == pytest.approx(eva.lifetime_gas_analyser_us(oxygen,water,212.125))
    assert k_water == pytest.approx(0.000093)
    assert meta['field_strength_vpercm'] == 212.125


def test_mean_concentrations_and_prm_are_not_replaced_by_workbook():
    rows=[]
    for t,o,w in [('2026-10-06T00:01:00-05:00',1.,20.),('2026-10-06T05:59:00-05:00',3.,40.)]:
        rows.extend([dict(timestamp=pd.Timestamp(t),quantity=q,source=q,value=v) for q,v in [('o2',o),('h2o',w),('prm_lifetime',1000*o)]])
    frame=pd.DataFrame(rows)
    policy={'minimum_span_hours':0,'max_gap_minutes':None}
    old={x['method']:x for x in aggregate_measurements([],measurements=frame,quality_config=policy)}
    new={x['method']:x for x in aggregate_measurements([],measurements=frame,quality_config=policy,conversion_config=CONFIG)}
    ko,kw,_=gas_conversion_rates(CONFIG)
    assert new['gas']['lifetime_us']==pytest.approx(1/(ko*2+kw*30))
    assert new['gas_o2']['lifetime_us']==pytest.approx(1/(ko*2))
    assert new['prm']==old['prm']
    assert new['gas']['concentration_ppb']==old['gas']['concentration_ppb']


def test_policy_survives_history_only_and_next_database_refresh(tmp_path):
    history=tmp_path/'history.json'
    update_json(history,[],gas_conversion_config=CONFIG)
    update_json(history,[])
    assert resolve_gas_conversion_config(None,history)==CONFIG


@pytest.mark.parametrize('config',[
    dict(CONFIG,field_strength_vpercm=0),dict(CONFIG,field_strength_vpercm=float('nan')),
    dict(CONFIG,offset_ppb=12.4),dict(CONFIG,model='unknown')])
def test_invalid_or_unmapped_calibration_is_explicit(config):
    with pytest.raises(ValueError):gas_conversion_rates(config)
