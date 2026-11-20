# Release validation: 19 November 2026

The repeated-cross-section DR influence correction is checked independently by
perturbing each cluster's empirical observation weights and refitting the entire
point estimator. The numerical derivatives match cluster influence contributions
and the full covariance of cohort-time and simple/group/calendar/event aggregates
under both control rules, with and without sampling weights. These regressions
would fail with the previous pre-period sign.

The corrected pre-period control-outcome coefficient loading is the treated-pre
covariate mean minus the propensity-weighted control-pre covariate mean. The old
formula reversed this loading. The outcome-regression coefficient itself and the
point estimate are unchanged.

Historical `DRDID` 1.3.0 / `did` 2.5.1 DR point estimates remain comparators, but
their historical DR influence formula shares the old sign. Current correctness
uses the empirical-weight derivatives and official Stata CSDID 2.0.3 instead.
The non-DR influence comparisons with the installed R references remain active.
Old point and exact-deletion fixtures are retained and explicitly checked during
DR fixture refresh; assignment matrices and generation seeds are unchanged.

## Independent external comparison

Official Stata CSDID 2.0.3 commit
`3bfcfbb53250b534f5c7d279969d14dba3bd834b` was compared with corrected Fast CSDID
on the same synthetic RCS samples, with 40 clusters, eight periods, two covariates,
never-treated controls and analytical inference. Sample sizes were 3,200, 32,000
and 320,000. The maximum ATT gap was 2.72e-11 and maximum analytical SE gap was
1.35e-12. This check concerns agreement, not a newly measured speedup.

| Observations | Official analytical SE | Corrected Fast analytical SE |
|---:|---:|---:|
| 3,200 | 0.0949831096136118 | 0.0949831096149549 |
| 32,000 | 0.0286030444628190 | 0.0286030444628184 |
| 320,000 | 0.0190933812088823 | 0.0190933812088804 |

The previous public comparator was Fast CSDID commit
`e34ad3f0e3a821119863c5503ef78f564eafd1b0`. The official correction is documented
in [csdid's NEWS](https://github.com/pedrohcgs/csdid-stata/blob/main/NEWS.md).

## Reproduce the checks

- Python: install pytest, then run `python -m pytest -q` in the Python package.
  The derivative regression is `tests/test_dr_influence_derivative.py`.
- R: run the `tests/testthat` suite from the package, or `R CMD check` on a source
  build. The independent derivative regression is `test-dr-influence-derivative.R`.
- Stata: from `validation`, run `do oracle.do`, `python check.py`,
  `do dr_derivative.do` and `do merit.do`.

The tests cover the documented stationary RCS scope. They do not establish causal
identification for a particular empirical application or feature parity with every
official CSDID option. Global packages, source datasets and running simulations were
not modified. Canonical dirty working trees were preserved; release work was isolated.

## This repository

Public candidate: **69 Python tests** and the complete R suite (**1779 expectations**) passed. R source build and `R CMD check --no-manual --no-vignettes` reported **Status: OK**. Stata passed all 16 frozen cases, independent derivative checks, MERIT51-deletion and public-API tests. Python and R README examples were executed. The public API source audit confirms that private inference extensions were not added.

The Python wheel built successfully and was installed into its own clean target.
A fresh subprocess imported that installed wheel and reproduced the MERIT point
estimate and exact jackknife standard error.
