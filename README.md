# Fast CSDID

**Estimate policy effects when different states adopt in different years.**

Suppose several states introduce merit scholarships in different years. You
want to estimate how college enrollment changes after adoption, using states
that have not adopted as a comparison group. Fast CSDID computes an effect for
each adoption cohort and year, then combines those effects into an overall
average treatment effect on the treated (ATT).

The package implements stationary repeated-cross-section equations. It does not
require following the same individuals over time. The validity of the comparison
still depends on the study's identifying assumptions, including appropriate
parallel trends; running the command does not establish those assumptions.

## Choose your language

Each implementation runs independently. The guides use the same merit-scholarship
data and explain the resulting estimates and standard errors.

| Language | Start here | Requirements |
|---|---|---|
| Stata | [Commands, example and output](stata/README.md) | Stata 16+ |
| R | [Installation and worked example](R/README.md) | R 4.1+; base R estimation |
| Python | [Installation and worked example](python/README.md) | Python 3.10+, NumPy, pandas |

## What do the results look like?

The bundled example uses college-enrollment shares and merit-scholarship adoption
in 51 jurisdictions during 1989–2000. The unadjusted simple ATT is **0.047782**,
or **4.78 percentage points**, with a delete-state jackknife SE of **0.014617**,
or **1.46 percentage points**. All three language guides reproduce these numbers.
The jackknife refits after each state is removed; it does not treat 612 state-year
rows as 612 independent clusters.

For example, after installing the R package:

```r
library(fastrcsdid)
d <- read.csv(system.file("extdata", "merit.csv", package="fastrcsdid"))
fit <- fast_csdid(d, "outcome", "year", "cohort", "state",
                  weights="weight", est_method="none", inference="jackknife")
fit$aggregation$estimates["simple"]
fit$jackknife$aggregate_covariance$se_cv3j["simple"]
```

The [language guides](#choose-your-language) give installation commands, explain
the data columns, show adjusted examples, and identify the returned results.

## Scope

Fast CSDID provides unadjusted, outcome-regression, IPW and doubly robust
estimation for stationary repeated cross sections. It reports cohort-time effects,
four aggregation families, analytical clustered inference, and ordinary
delete-cluster jackknife inference. It does not implement the panel equations or
every option of the official CSDID packages.

R uses exact moment updates for linear estimators. Its nonlinear jackknife offers
both linearized and exact-refit calculations; request `jackknife_method="refit"`
for an exact comparison with Stata or Python. Performance depends on the estimator,
sample and inference requested; [benchmark reports](benchmarks/README.md) state
those choices explicitly.

## Correction and validation

Release 1.0.1 corrects the pre-period outcome-regression contribution to
DR analytical inference in all three languages. Point estimates and exact
jackknife estimates are unchanged. Recompute DR analytical standard errors made
with 1.0.0. The correction is checked against independent weight derivatives and
official Stata CSDID 2.0.3; see [validation](VALIDATION.md) and
[change history](CHANGELOG.md). The published [v1.0.1 release](https://github.com/mattdwebb/fastcsdid/releases/tag/v1.0.1) is tagged at `bf4a7ad`.

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

[Source provenance](PROVENANCE.md) records the implementation lineage.
