# Fast CSDID in Stata

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

Stata 16 or later is required. No other ado packages are needed.

```stata
net install fastcsdid, from("https://raw.githubusercontent.com/mattdwebb/fastcsdid/main/stata") replace
```

For the bundled example, also download the repository using **Code → Download ZIP**,
extract it, and change Stata's working directory to its `stata` folder.

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

```stata
set type double
import delimited using "validation/frozen/merit.csv", clear
fastcsdid outcome, time(year) gvar(cohort) cluster(state) weights(weight) ///
    method(none) jackknife
display "Jackknife SE = " %9.6f e(jk_se)
display "Jackknife p-value = " %9.6f e(jk_p)
```

Actual numerical output:

```text
RCS aggregate ATT =   0.04778200
Native SE (normal reference) =   0.01305415
Jackknife SE =  0.014617
Jackknife p-value =  0.001957
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

With covariates in your own data, add `covariates(age education)` and
`method(dripw)`; use the actual variable names. Encode categorical variables as
appropriate indicators first. The teaching MERIT fixture has no individual-level
covariates.

`matrix list e(attgt)` displays cohort-time effects; `matrix list e(aggregation)`
shows the four aggregation families. `e(b)` and `e(V)` contain the simple ATT
and its analytical variance. Jackknife results are separate: `e(jk_se)`,
`e(jk_p)`, `e(jk_df)`, and `e(deleted)`. Requesting `jackknife` does not replace
`e(V)` with the jackknife variance.

Read `help fastcsdid` and [development history](DEVELOPMENT_REFERENCE.md).
From `validation`, `do oracle.do` followed by `python check.py` checks frozen
numerical comparisons; `do dr_derivative.do` checks analytical covariance against
independent weight perturbations. Tests require a licensed Stata installation.

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
