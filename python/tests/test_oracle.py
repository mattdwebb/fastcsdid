from pathlib import Path
import numpy as np
import pandas as pd
import pytest
from fastcsdidpy import fit, rif, influence_function

ROOT = Path(__file__).resolve().parents[1] / "validation" / "frozen"


@pytest.fixture
def data():
    return pd.read_csv(ROOT / "fixture.csv")


@pytest.mark.parametrize("method", ["none", "reg", "ipw", "dr"])
@pytest.mark.parametrize("control", ["never", "not_yet"])
@pytest.mark.parametrize("weighted", [0, 1])
def test_frozen_r(data, method, control, weighted):
    name = f"{method}_{control}_{weighted}"
    r = fit(data, "y", "time", gvar="g", cluster="cluster", method=method,
            covariates=None if method == "none" else ["x", "z"],
            control_group=control, weights="w" if weighted else None, inference="jackknife")
    def read(suffix):
        return pd.read_csv(ROOT / f"{name}_{suffix}.csv")
    np.testing.assert_allclose(r["group_time"].estimate, read("att").estimate, atol=2e-9, rtol=0)
    np.testing.assert_allclose(r["aggregation"].estimate, read("agg").estimate, atol=2e-9, rtol=0)
    np.testing.assert_allclose(r["influence_att"], read("if"), atol=2e-8, rtol=0)
    np.testing.assert_allclose(r["influence_aggregate"], read("aif"), atol=2e-8, rtol=0)
    np.testing.assert_allclose(r["covariance"], read("cov"), atol=2e-9, rtol=0)
    np.testing.assert_allclose(r["jackknife"]["deleted_estimates"], read("jk").deleted, atol=2e-8, rtol=0)
    np.testing.assert_allclose(rif(r).mean(axis=0), r["group_time"].estimate, atol=1e-10)
    np.testing.assert_allclose(influence_function(r, cluster_sum=True).sum(axis=0), 0, atol=1e-9)


def test_official_drdid():
    checks = pd.read_csv(ROOT / "external_drdid.csv")
    assert len(checks) == 48
    assert checks.att_gap.max() < 2e-8
    assert checks.if_gap.max() < 2e-7


def test_validation_and_modes(data):
    with pytest.raises(NotImplementedError, match="panel"):
        fit(data,"y","time",unit="cluster",gvar="g",cluster="cluster")
    with pytest.raises(ValueError, match="ignore"):
        fit(data,"y","time",gvar="g",cluster="cluster",covariates=["x"],method="none")
    with pytest.raises(TypeError, match="Unimplemented"):
        fit(data,"y","time",gvar="g",cluster="cluster",invented=True)
    with pytest.raises(NotImplementedError, match="anticipation"):
        fit(data,"y","time",gvar="g",cluster="cluster",anticipation=1)
    d = data.copy()
    d.loc[0,"g"] = 5
    with pytest.raises(ValueError, match="constant"):
        fit(d,"y","time",gvar="g",cluster="cluster")


def test_binary_and_default_mapping(data):
    d = data.assign(D=((data.g>0)&(data.time>=data.g)).astype(int))
    common = dict(cluster="cluster",covariates=["x","z"])
    a = fit(d,"y","time",gvar="g",**common)
    b = fit(d,"y","time",treatment="D",method="dr",**common)
    assert a["estimate"] == pytest.approx(b["estimate"])
    d.loc[(d.cluster==1)&(d.time==4),"D"] = 0
    with pytest.raises(ValueError, match="absorbing"):
        fit(d,"y","time",treatment="D",**common)






def test_fixed_support_failed_deletion_not_dropped(data):
    d = data[~data.cluster.isin([2,3,4])]
    r = fit(d,"y","time",gvar="g",cluster="cluster",method="none",inference="jackknife",
            jackknife_support="fixed")
    assert r["jackknife"]["status"] == "IDENTIFICATION_FAILURE"
    assert np.isnan(r["jackknife"]["p"])


def test_window_and_unsupported_cell(data):
    r = fit(data,"y","time",gvar="g",cluster="cluster",event_window=(0,0))
    assert (r["group_time"].event == 0).all()
    d = data[~((data.g==3)&(data.time==4))]
    r = fit(d,"y","time",gvar="g",cluster="cluster",method="none")
    assert r["status"] == "IDENTIFICATION_FAILURE"
    assert np.isfinite(r["aggregation"].set_index("label").loc["calendar:3","estimate"])
