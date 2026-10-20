"""Development-only comparison of actual Stata batch output to frozen R."""
from pathlib import Path
import json
import numpy as np
import pandas as pd

here = Path(__file__).resolve().parent
root = here/"frozen"
rows = pd.read_csv(here/"actual.csv")
assert len(rows) == 16
gaps = []
for row in rows.itertuples():
    name = f"{row.method}_{row.control}_{row.weighted}"
    expected = pd.read_csv(root/"expected.csv")
    ref = expected[(expected.method==row.method)&(expected.control==row.control)&(expected.weighted==row.weighted)].iloc[0]
    gaps.append({"case":name,"att_gap":abs(row.estimate-ref.estimate),"se_gap":abs(row.se-ref.se),
                 "jk_gap":abs(row.jk-ref.jk)})
maximum = {k:max(r[k] for r in gaps) for k in ("att_gap","se_gap","jk_gap")}
print(maximum)
assert max(maximum.values()) < 2e-8
actual_components = pd.read_csv(here/"components.csv")
expected_components = pd.read_csv(root/"components.csv")
pd.testing.assert_frame_equal(actual_components.drop(columns="estimate"),expected_components.drop(columns="estimate"))
maximum["component_gap"] = float(np.max(np.abs(actual_components.estimate-expected_components.estimate)))
assert maximum["component_gap"] < 2e-8
(here/"results.json").write_text(json.dumps({"status":"pass","cases":gaps,"maxima":maximum},indent=2)+"\n")
