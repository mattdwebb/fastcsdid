from concurrent.futures import ProcessPoolExecutor
import multiprocessing
from pathlib import Path
import numpy as np
import pandas as pd
import pytest
from fastcsdidpy import fit
from fastcsdidpy.probability import student_two_sided

ROOT = Path(__file__).resolve().parents[1] / "validation" / "frozen"


def test_student_tail_r_oracle():
    oracle = pd.read_csv(ROOT/"student_t.csv")
    observed = [student_two_sided(r.statistic,r.df) for r in oracle.itertuples()]
    np.testing.assert_allclose(observed,oracle.p,atol=3e-13,rtol=3e-12)




def _worker(method):
    d = pd.read_csv(ROOT/"fixture.csv")
    r = fit(d,"y","time",gvar="g",cluster="cluster",method=method,covariates=["x","z"])
    return r["estimate"],r["native_se"]


def test_spawn_isolation():
    with ProcessPoolExecutor(max_workers=2,mp_context=multiprocessing.get_context("spawn")) as executor:
        results = list(executor.map(_worker,["reg","dr"]))
    np.testing.assert_allclose(results,[_worker("reg"),_worker("dr")],atol=1e-12)


@pytest.mark.parametrize("method", ["none", "reg", "ipw", "dr"])
def test_fixed_design_linearity(method):
    d = pd.read_csv(ROOT/"fixture.csv")
    kw = dict(gvar="g",cluster="cluster",method=method,covariates=None if method=="none" else ["x","z"])
    d["other"] = np.sin(d.y*2.3)
    d["linear"] = 1.4*d.y-.7*d.other
    a,b,c = [fit(d,col,"time",**kw)["estimate"] for col in ("y","other","linear")]
    assert c == pytest.approx(1.4*a-.7*b,abs=2e-12)


def test_zero_weights_and_singularity():
    d = pd.read_csv(ROOT/"fixture.csv")
    d.loc[d.cluster==1,"w"] = 0
    r = fit(d,"y","time",gvar="g",cluster="cluster",weights="w",covariates=["x","z"])
    assert r["status"] == "pass"
    d["duplicate"] = d.x
    r = fit(d,"y","time",gvar="g",cluster="cluster",covariates=["x","duplicate"])
    assert r["status"] == "IDENTIFICATION_FAILURE"
    assert np.isnan(r["estimate"])


def test_propensity_invariant_to_replicated_observations():
    """Summed-objective roundoff must not prevent a large-sample logit fit."""
    from fastcsdidpy.core import _logit
    rng = np.random.default_rng(20261020)
    X = np.column_stack([np.ones(160), rng.normal(size=160), rng.binomial(1, .4, 160)])
    D = rng.binomial(1, 1 / (1 + np.exp(-X @ np.array([.1, .2, -.3]))))
    w = np.ones(160)
    p, _ = _logit(X, D, w)
    large, _ = _logit(np.repeat(X, 500, axis=0), np.repeat(D, 500), np.repeat(w, 500))
    np.testing.assert_allclose(large[::500], p, atol=2e-10, rtol=0)
