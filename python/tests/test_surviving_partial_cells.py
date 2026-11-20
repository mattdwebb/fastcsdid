from pathlib import Path
import numpy as np
import pandas as pd
import pytest
from fastcsdidpy import fit

@pytest.mark.parametrize("method",["none","reg","ipw","dr"])
def test_changing_support_recomputes_surviving_att(method):
    d=pd.read_csv(Path(__file__).parent/"fixtures/partial_deletion.csv")
    options=dict(gvar="g",cluster="cluster",weights="w",method=method,
                 covariates=None if method=="none" else ["x","z"],data_mode="rcs")
    r=fit(d,"y","time",inference="jackknife",**options)
    jk=r["jackknife"]
    assert jk["status"]=="pass" and len(jk["deleted_estimates"])==13
    assert np.isfinite(jk["deleted_estimates"]).all() and jk["df"]==12
    assert jk["support_changed"][list(jk["delete_clusters"]).index(13)]
    broken=fit(d[d.cluster!=13],"y","time",inference="point",**options)
    assert broken["status"]=="IDENTIFICATION_FAILURE"
    fixed=fit(d,"y","time",inference="jackknife",jackknife_support="fixed",**options)
    assert fixed["jackknife"]["status"]=="IDENTIFICATION_FAILURE"
