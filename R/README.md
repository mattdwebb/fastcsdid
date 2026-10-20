# R: fastrcsdid 1.0.0

Install from the repository root with `R CMD INSTALL R`.

```r
library(fastrcsdid)
data(fastrcsdid_example)
fit <- fast_csdid(
  fastrcsdid_example, "y", "time", "cohort", "region",
  weights = "weight", xformla = ~x1 + x2, est_method = "dr",
  control_group = "never", inference = "jackknife",
  jackknife_method = "refit"
)
fit$aggregation$estimates
inference_comparison(fit)
fit$jackknife$aggregate_covariance$se_cv3j
```

`est_method` is `"none"`, `"reg"`, `"ipw"`, or `"dr"`. Covariate methods require
`xformla` (use `~1` for an intercept-only model). `inference="none"` computes
point estimates and analytical cluster influence covariance;
`inference="jackknife"` additionally computes cluster deletions.

For regression and no-covariate models, deletion is exact. Nonlinear models
default to `jackknife_method="linearized"`; `"refit"` computes exact refits and
`"both"` retains both **jackknife algorithms**. The separate `inference="both"`
option from development versions is no longer supported.

`jackknife_support="changing"` recomputes weights and permits zero-mass cohorts
to disappear if the aggregate remains identified. `"fixed"` additionally
requires the original contributing group-time cells to survive. Consult support
diagnostics before interpreting uncertainty. Mean-centered CV3J and
full-estimate-centered CV3 outputs have distinct names; use the intended
variance convention explicitly.

`collapse_rcs()` and `fit_csdid_moments()` provide a pre-collapsed input route.
Public timings include raw-data preparation unless explicitly labeled otherwise.

To test from the repository root:

```r
pkgload::load_all("R")
testthat::test_dir("R/tests/testthat", stop_on_failure = TRUE)
```

Optional equivalence tests use the R `did` package. Base estimation has no
non-base runtime dependency.
