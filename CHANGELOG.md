# Changelog

## 1.0.1 candidate — 2026-11-19

- Correct the pre-period outcome-regression contribution to repeated-cross-section
  doubly robust analytical influence functions in Stata, R and Python. Point
  estimates and exact delete-cluster jackknife estimates are unchanged. Recompute
  DR analytical standard errors, covariance and associated tests from 1.0.0.
- Validate the corrected influence and all aggregation covariance entries using
  empirical cluster-weight derivatives, independently of copied influence formulas.
  Corrected Stata analytical inference matches official csdid 2.0.3.
- Refresh obsolete DR influence fixtures and explain why DRDID 1.3.0's historical
  DR influence is no longer the correctness oracle. Retain point/deletion fixtures.
- Replace developer-first landing pages with installation, runnable MERIT examples,
  actual results, plain-language interpretation, and explicit supported scope.
- Bump package metadata and CI build filename consistently; no release tag created.

## 1.0.0 — 2026-10-20

- Initial release of Stata, R, and Python implementations.
- Unadjusted, regression, IPW, and doubly robust estimation for stationary
  repeated cross sections.
- Analytical cluster inference and delete-one-cluster jackknife inference.
- Sampling weights, never-treated and not-yet-treated controls, and simple,
  cohort, calendar, and event-time aggregations.
- Reproducible timings against official Stata csdid 2.0.0.
- Scale-aware logit line searches in Python and Stata for large samples.

DRDID 1.4.0 and 1.4.1 include the corresponding upstream sign correction. Historical references to 1.3.0 do not apply to those versions. Both exact-refit and R linearized deletion estimates are unchanged.
