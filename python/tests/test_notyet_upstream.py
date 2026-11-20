from pathlib import Path
import pandas as pd
import numpy as np
from fastcsdidpy import fit

def test_notyet_dripw_matches_verified_upstream_analytical_se():
    root=Path(__file__).parent/"fixtures"
    data=pd.read_csv(root/"notyet_upstream_data.csv")
    reference=pd.read_csv(root/"notyet_upstream_reference.csv")
    result=fit(data,y="y",time="time",gvar="g",cluster="cluster",covariates=["x","z"],data_mode="rcs",control_group="not_yet",method="dripw")
    np.testing.assert_allclose(result["estimate"],reference.estimate,atol=1e-8,rtol=0)
    np.testing.assert_allclose(result["native_se"],reference.native_se,atol=1e-8,rtol=0)
