# Python: fastcsdidpy 1.0.0

Install from the repository root with `python -m pip install "./python[test]"`.
Run the suite with `python -m pytest python/tests -q`.

```python
import pandas as pd
from fastcsdidpy import fit

d = pd.read_csv("python/validation/frozen/fixture.csv")
result = fit(d, "y", "time", gvar="g", cluster="cluster",
             covariates=["x", "z"], weights="w", method="dr",
             data_mode="rcs", control_group="never", inference="standard")
print(result["estimate"], result["native_se"], result["native_p"])
print(result["jackknife"]["se"], result["jackknife"]["p"])
```

Methods: `none`, `reg`/`ra`, `ipw`/`stdipw`, `dr`/`dripw`. Encode categorical
covariates before fitting. `method="none"` rejects supplied covariates.

Inference choices: `native` (default), `jackknife`, `standard` (native plus
jackknife), `point`, or `none`. Analytical covariance is computed even under
`point`/`none`; those settings omit additional inference. `null` changes the
coefficient hypothesis. Unsupported inference choices and options raise errors.

The jackknife is a complete refit with mean-centered CV3J variance and
Student t(G−1) inference. `jackknife_support="changing"` recomputes cohort
weights and allows a singleton cohort to disappear when the aggregate remains
identified. `"fixed"` requires the full-sample group-time support to survive.
Failed deletions remain failures; they are not dropped.

`treatment="D"` can replace `gvar` for an absorbing binary history.
`event_window=(0, 2)` restricts post-treatment comparisons.
`data_mode="auto"` rejects repeated unit IDs as unsupported panel data;
explicit `"rcs"` selects the RCS equations.

Results include group-time and aggregation DataFrames, native covariance,
observation/cluster influence functions, diagnostics, and optional jackknife
results. `influence_function(result)` and `rif(result)` expose influence arrays.
Invalid identification retains explicit failure status and NaN uncertainty.

The frozen R/DRDID fixtures retain their original generating version metadata;
they were not regenerated to make the v1 implementation pass.
