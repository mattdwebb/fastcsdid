from pathlib import Path
import pandas as pd
from fastcsdidpy import fit

data = pd.read_csv(Path(__file__).resolve().parents[1]/"validation/frozen/fixture.csv")
result = fit(data,"y","time",gvar="g",cluster="cluster",covariates=["x","z"],
             weights="w",inference="standard",data_mode="rcs")
print(result["aggregation"].to_string(index=False))
print("Native SE:",result["native_se"],"CV3J SE:",result["jackknife"]["se"])
