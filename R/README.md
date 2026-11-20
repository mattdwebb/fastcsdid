# Fast CSDID in R

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

R 4.1 or later is required. Estimation uses base R; test and vignette packages
are optional.

```r
install.packages("remotes")
remotes::install_github("mattdwebb/fastcsdid", subdir = "R")
```

From a downloaded repository, `R CMD INSTALL R` is the local alternative.

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

```r
library(fastrcsdid)
data <- read.csv(system.file("extdata", "merit.csv", package = "fastrcsdid"))
fit <- fast_csdid(data, "outcome", "year", "cohort", "state",
                  weights = "weight", est_method = "none",
                  inference = "jackknife")
att <- unname(fit$aggregation$estimates["simple"])
se <- unname(fit$jackknife$aggregate_covariance$se_cv3j["simple"])
cat(sprintf("ATT: %.6f\nJackknife SE: %.6f\n", att, se))
write.csv(fit$aggregation$metadata, "merit_effects.csv", row.names = FALSE)
```

Actual output:

```text
ATT: 0.047782
Jackknife SE: 0.014617
```

The estimate of **0.047782 is 4.78 percentage points of enrollment**, not a
4.78 percent proportional increase. The ordinary delete-state jackknife standard
error is **1.46 percentage points**. It measures uncertainty around the same
estimate by leaving each of the 51 jurisdictions out in turn and refitting.
Four deletions remove some cohort-time effects; their surviving aggregates remain
estimable and all 51 deletions are included. A required deletion that cannot
produce an aggregate is a reported failure, not a draw silently discarded.

## Add regression adjustment

For a fully runnable adjusted example, the package also includes small synthetic
individual-level data with covariates:

```r
data(fastrcsdid_example)
adjusted <- fast_csdid(fastrcsdid_example, "y", "time", "cohort", "region",
                       weights = "weight", xformla = ~x1 + x2,
                       est_method = "dr", inference = "jackknife",
                       jackknife_method = "refit")
inference_comparison(adjusted)
```

`xformla` names the adjustment variables. For nonlinear IPW/DR models,
`jackknife_method="refit"` refits every deleted propensity model to convergence.
The default `"linearized"` option uses a faster approximation; `"both"` retains
both so you can compare them. Linear/unadjusted deletion updates are exact.

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
Student t reference with Gâˆ’1 degrees of freedom. It recomputes nuisance models
and aggregation weights after deletion. A larger or smaller standard error is
not, by itself, evidence that one procedure is always better.

Supported results include post-treatment cohort-time effects and simple, cohort,
calendar, and event-time aggregations. The package does not implement the
panel-specific equations, anticipation, pre-treatment event-study coefficients,
or every option of the official CSDID packages. Keep the estimand and observation
weights explicit when comparing results with other software.

`summary(fit)` reports effects, while `jackknife_influence_summary(fit)` and
`jackknife_cluster_table(fit, estimand="simple")` help identify influential states.
Use `?fast_csdid` for arguments. [Development history](DEVELOPMENT_REFERENCE.md)
contains advanced input and inference workflows. The regression suite is under
`tests/testthat`; development tools are not runtime dependencies.

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
