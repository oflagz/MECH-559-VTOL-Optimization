"""Regression checks for independent solver, Pareto, and interval outputs."""
import csv
import json
from pathlib import Path
import numpy as np
import pytest
from analysis.advanced_validation import dominance_mask

ROOT=Path(__file__).resolve().parents[1]
TAB=ROOT/'results/tables'


def test_true_pareto_dominance_with_ties():
                                                                               
    x=np.array([[1,5],[2,4],[3,3],[2,6],[1,5],[3,5]],float)
    mask=dominance_mask(x)
    assert mask.tolist()==[True,True,True,False,True,False]


def test_dominance_rejects_nonfinite_data():
    with pytest.raises(ValueError): dominance_mask([[1,float('nan')]])


def test_multiobjective_archive_is_feasible_and_nondominated_points_are_not_dominated():
    rows=list(csv.DictReader((TAB/'pareto_qualified_archive.csv').open()))
    assert len(rows)>10
    obj=np.array([[float(x['hover_wh']),float(x['cruise_wh'])] for x in rows])
    mask=dominance_mask(obj)
    assert mask.tolist()==[x['nondominated']=='True' for x in rows]
    assert all(float(x['total_wh'])==pytest.approx(float(x['hover_wh'])+float(x['cruise_wh']),abs=1e-6) for x in rows)


def test_interval_corner_count_and_no_stochastic_claims():
    rows=list(csv.DictReader((TAB/'interval_uncertainty_corners.csv').open()))
    corners=[v for v in rows if v['scenario']!='nominal']
    assert len(corners)==16
    assert all(x['trim_found']=='True' for x in corners)
    assert any(x['all_constraints_satisfied']=='False' for x in corners)
    v=json.loads((TAB/'advanced_validation_summary.json').read_text())
    assert v['uncertainty']['corner_count']==len(corners)
    assert v['uncertainty']['qualified_count']==sum(x['all_constraints_satisfied']=='True' for x in corners)
    assert 'NOT confidence intervals' in v['interpretation']['uncertainty']


def test_derivative_free_comparison_never_reports_infeasible_as_feasible():
    rows=list(csv.DictReader((TAB/'derivative_free_comparison.csv').open()))
    assert len(rows)==6
    for r in rows:
        if r['feasible']=='True':
            assert float(r['max_ineq_residual']) <= 1.01e-5
            assert float(r['max_trim_residual']) < 2.01e-4
    r=json.loads((TAB/'advanced_validation_summary.json').read_text())
    assert r['dfo_qualified']==sum(x['feasible']=='True' for x in rows)
