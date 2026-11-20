# Fast CSDID in Python

Suppose several states introduce merit scholarships in different years. You
want to estimate how college enrollment changes after adoption, using states
that have not adopted as a comparison group. Fast CSDID computes an effect for
each adoption cohort and year, then combines those effects into an overall
average treatment effect on the treated (ATT).

The package implements stationary repeated-cross-section equations. It does not
require following the same individuals over time. The validity of the comparison
still depends on the study's identifying assumptions, including appropriate
parallel trends; running the command does not establish those assumptions.

## Install

Python 3.10 or later, NumPy and pandas are required.
From the downloaded `fastcsdid` repository root, run
`python -m pip install ./python`. The public GitHub installation is
`python -m pip install "git+https://github.com/mattdwebb/fastcsdid.git#subdirectory=python"`.
Run the example below from the downloaded repository root.

## A merit-scholarship example

The bundled MERIT example has 612 state-year rows, covering 51 state identifiers
from 1989 through 2000. Each row contains:

| Column | Meaning |
|---|---|
| `state` | Jurisdiction identifier; inference is clustered at this level |
| `year` | Observation year |
| `cohort` | First scholarship year, repeated on every row; 0 means never treated |
| `outcome` | College-enrollment share, between 0 and 1 |
| `weight` | Number of underlying observations in that state-year |

These are collapsed data from the merit-scholarship example used in the UN-DID
paper, originally distributed by Conley and Taber (2011). The [original data
source](https://economics.uwo.ca/people/conley_docs/code_to_download.html) and the
included fixture retain the original state codes. This compact example is for
learning the commands; it is not a new empirical analysis.

```python
import pandas as pd
from fastcsdidpy import fit

data = pd.read_csv("python/validation/frozen/merit.csv")
result = fit(data, "outcome", "year", gvar="cohort", cluster="state",
             weights="weight", method="none", inference="jackknife")
print(f"ATT: {result['estimate']:.6f}")
print(f"Jackknife SE: {result['jackknife']['se']:.6f}")
print(f"Jackknife p-value: {result['jackknife']['p']:.6f}")
result["aggregation"].to_csv("merit_effects.csv", index=False)
```

Actual output:

```text
ATT: 0.047782
Jackknife SE: 0.014617
Jackknife p-value: 0.001957
```

The estimate of **0.047782 is 4.78 percentage points of enrollment**, not a
4.78 percent proportional increase. The ordinary delete-state jackknife standard
error is **1.46 percentage points**. It measures uncertainty around the same
estimate by leaving each of the 51 jurisdictions out in turn and refitting.
Four deletions remove some cohort-time effects; their surviving aggregates remain
estimable and all 51 deletions are included. A required deletion that cannot
produce an aggregate is a reported failure, not a draw silently discarded.

## Choose an estimator and inference

Use the unadjusted estimator when no additional covariates are required. `reg`
uses outcome-regression adjustment, `ipw` uses propensity-score weighting, and
`dr` (`dripw`) combines the two. Supply pre-specified, appropriately coded
covariates for adjusted methods; the program does not decide which controls are
causally appropriate. The default comparison is never-treated units; `not_yet`
also permits later-treated cohorts while they remain eligible.

Analytical clustered standard errors and jackknife standard errors are separate
outputs. The analytical calculation uses cluster-summed influence functions and
a normal reference. The ordinary jackknife uses mean-centered CV3J and a
Student t reference with G−1 degrees of freedom. It recomputes nuisance models
and aggregation weights after deletion. A larger or smaller standard error is
not, by itself, evidence that one procedure is always better.

Supported results include post-treatment cohort-time effects and simple, cohort,
calendar, and event-time aggregations. The package does not implement the
panel-specific equations, anticipation, pre-treatment event-study coefficients,
or every option of the official CSDID packages. Keep the estimand and observation
weights explicit when comparing results with other software.

The public Python API supports point, native, jackknife and standard inference.
`standard` combines native and jackknife results. Unsupported requests fail
explicitly.

For a new dataset with numeric covariates, use
`covariates=["age", "education"]` and `method="dr"`; replace those names with
actual columns in your data. MERIT's collapsed teaching fixture has no individual
covariates, so do not treat its state-year averages as individual-level controls.

Use `result["group_time"]` for component effects, `result["aggregation"]` for
summaries, and `result["jackknife"]["deleted_estimates"]` for the complete deletion
vector. The returned `metadata` records the selected target and options.
Read [development history](DEVELOPMENT_REFERENCE.md) for assignment layouts and
additional options. Run `python -m pytest` from the Python package directory to
execute the regression suite after installing pytest.

## Methods, help, and license

- Callaway and Sant'Anna (2021), [Difference-in-Differences with Multiple Time
  Periods](https://doi.org/10.1016/j.jeconom.2020.12.001).
- Sant'Anna and Zhao (2020), [Doubly Robust Difference-in-Differences
  Estimators](https://doi.org/10.1016/j.jeconom.2020.06.003).
- Karim, Nielsen, MacKinnon, and Webb (2026), [Improved Inference for CSDID
  Using the Cluster Jackknife](https://arxiv.org/abs/2602.12043).

Developed by Matthew D. Webb. MIT licensed; see the repository's license and
provenance files. Numerical checks and known scope are documented in
`VALIDATION.md`; development history is separate from this introductory guide.
