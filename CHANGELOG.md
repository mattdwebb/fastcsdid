# Changelog

## 1.0.0 — 2026-10-20

- Consolidate independent Stata, R, and Python implementations.
- Retain RCS estimation, analytical cluster inference, and cluster jackknife.
- Remove randomization inference, RI-jackknife, assignment generators, their
  API options, and implementation source. These features are deferred to v2.
- Reject removed options rather than silently substituting another procedure.
- Add reproducible timings against the pinned official Stata csdid 2.0.0.
- Scale the logit line-search objective comparison to double-precision
  roundoff in Python and Stata, preventing large-sample numerical stalls
  while retaining the original convergence and identification thresholds.
