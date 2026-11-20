from pathlib import Path
import json
import numpy as np
import pandas as pd
from fastcsdidpy import fit

root=Path(__file__).resolve().parent
data=pd.read_csv(root/"data.csv")
result=fit(data,y="y",time="time",gvar="g",cluster="cluster",covariates=["x","z"],method="dripw",control_group="not_yet",data_mode="rcs",inference="native")
r=pd.read_csv(root/"R_results.csv").drop(columns="version")
stata=pd.read_csv(root/"Stata_results.csv")
out=pd.concat([r,stata,pd.DataFrame([dict(engine="fast_Python",estimate=result["estimate"],native_se=result["native_se"])])],ignore_index=True)
for col in ["estimate","native_se"]:
    np.testing.assert_allclose(out[col],out[col].iloc[0],rtol=0,atol=1e-8)
out.to_csv(root/"comparison.csv",index=False,float_format="%.17g")
summary=dict(status="PASS",control_group="not_yet",method="DRIPW",n_observations=len(data),n_clusters=data.cluster.nunique(),max_estimate_gap=float(np.ptp(out.estimate)),max_native_se_gap=float(np.ptp(out.native_se)),engines=out.engine.tolist())
(root/"validation.json").write_text(json.dumps(summary,indent=2),encoding="utf-8")
print(out.to_string(index=False))
print(json.dumps(summary))
