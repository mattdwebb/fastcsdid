# Version 1 validation

Version 1.0.0 was validated on 2026-10-20 against frozen numerical references
and independent implementations.

| Check | Result |
|---|---|
| Python estimator and API suite | 65 tests passed |
| R built-package testthat suite | 1,747 expectations passed; no failures, warnings, or skips |
| `R CMD check --no-manual` | Status: OK |
| Stata frozen-oracle comparisons | All 16 method/control/weight combinations passed |
| Stata post-treatment component estimates | Maximum absolute gap 6.67e-15 |
| Stata aggregate ATT / native SE / CV3J SE | Maximum gaps 5.33e-15 / 6.25e-16 / 2.46e-13 |
| Stata argument validation and support-policy tests | Passed |
| Stata MERIT changing-support example | Passed; 51 deletions |
| Python built-wheel import and native/jackknife fit | Passed |
| Stata installation from the public GitHub URL into an isolated ado directory | Passed |
| Source and exported-API checks | Passed |
| Matched native-inference timing run | 192 calls; maximum ATT gap 8.19e-11; maximum SE gap 4.52e-12 |

The R suite includes independent equivalence checks against installed
`did`/`DRDID`, exact/linearized jackknife tests, influence functions,
sampling weights, and support failures. Python and Stata use frozen
cross-language numerical references. API checks verify supported inference
methods and rejection of unsupported arguments.

Large-sample validation includes IPW fits with 320,000 observations. Stata and
Python use scale-aware logit line-search comparisons to accommodate roundoff
in the summed likelihood, with separate Newton-step convergence and
identification checks.

These are validations of the documented RCS scope, not claims of full feature
parity with the official `csdid`/`did` packages. Native analytical inference and
cluster-jackknife inference are distinct outputs. In particular, a runtime
comparison with native analytical csdid does not establish jackknife speed.

GitHub Actions runs Python on 3.10, 3.12, and 3.14 and runs R validation.
Stata tests require a licensed installation and are supplied as do-files.
See [benchmarks/results/TIMINGS.md](benchmarks/results/TIMINGS.md) for
separately measured timing and agreement results.
