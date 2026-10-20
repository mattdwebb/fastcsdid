# Version 1 validation

Validated locally on 2026-10-20. Numerical fixtures were retained from the
development repositories rather than regenerated against the edited code.

| Check | Result |
|---|---|
| Python retained-estimator and release-API suite | 65 tests passed |
| R built-package testthat suite | 1,747 expectations passed; no failures, warnings, or skips |
| `R CMD check --no-manual` | Status: OK |
| Stata frozen-oracle comparisons | All 16 method/control/weight combinations passed |
| Stata post-treatment component estimates | Maximum absolute gap 6.67e-15 |
| Stata aggregate ATT / native SE / CV3J SE | Maximum gaps 5.33e-15 / 6.25e-16 / 2.46e-13 |
| Stata removed-options and support-policy tests | Passed |
| Stata MERIT changing-support example | Passed; 51 deletions |
| Python built-wheel import and native/jackknife fit | Passed |
| Runtime/API source audit for removed v2 features | Passed |

The R suite includes independent equivalence checks against installed
`did`/`DRDID`, ordinary exact/linearized jackknife tests, influence functions,
sampling weights, and support failures. Python and Stata retain frozen
cross-language numerical fixtures. The source audit checks the executable
implementations and exported API rather than just help-file visibility.

The large-N benchmark revealed a numerical logit line-search stall in the
Stata port: an absolute objective tolerance was smaller than roundoff in a
summed likelihood. Stata and Python now use a scale-aware floating-point
comparison, retaining the original Newton-step convergence threshold,
separation checks, and iteration cap. The formerly failing 320,000-observation
IPW fit succeeds, and the frozen-oracle comparisons remain within tolerance.

These are validations of the documented RCS scope, not claims of full feature
parity with the official `csdid`/`did` packages. Native analytical inference and
cluster-jackknife inference are distinct outputs. In particular, a runtime
comparison with native analytical csdid does not establish jackknife speed.

GitHub Actions runs Python on 3.10, 3.12, and 3.14 and runs R validation.
Stata tests require a licensed installation and are supplied as do-files.
See [benchmarks/results/TIMINGS.md](benchmarks/results/TIMINGS.md) for
separately measured timing and agreement results.
