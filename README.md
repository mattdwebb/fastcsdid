# fastcsdid

Fast group-time difference-in-differences for **stationary repeated cross
sections**, with Stata, R, and Python implementations in one repository.
Version 1.0.0 provides unadjusted, outcome-regression, inverse-probability-weighted,
and doubly robust estimation, analytical cluster inference, and ordinary
delete-one-cluster jackknife inference.

| Implementation | Directory | Package/command | Requirement |
|---|---|---|---|
| Stata | [stata](stata) | `fastcsdid` | Stata 16+ |
| R | [R](R) | `fastrcsdid` | R 4.1+ |
| Python | [python](python) | `fastcsdidpy` | Python 3.10+, NumPy, pandas |

The existing language package names are retained for compatibility. There is one
Git repository and one release version. Each implementation runs independently.

## Install

```stata
net install fastcsdid, from("https://raw.githubusercontent.com/mattdwebb/fastcsdid/main/stata") replace
```

```r
remotes::install_github("mattdwebb/fastcsdid", subdir = "R")
```

```sh
python -m pip install "git+https://github.com/mattdwebb/fastcsdid.git#subdirectory=python"
```

For reproducible analyses, replace `main` in the Stata URL, or append `@COMMIT`
to the GitHub repository reference in R/Python, with a recorded commit SHA.
Local installation: `R CMD INSTALL R` or `python -m pip install ./python`.

## Scope and inference

Treatment is absorbing and assigned at the cluster/jurisdiction level; each
cluster has one adoption cohort. Zero denotes never treated. Covariates are
contemporaneous repeated-cross-section covariates, not frozen panel baselines.
Never-treated and not-yet-treated control groups are supported. The common
outputs are post-treatment ATT(g,t) and simple, cohort, calendar, and event-time
aggregations.

Analytical standard errors use cluster-summed influence functions. The ordinary
cluster jackknife remains available. Stata and Python completely refit each
deletion. R uses exact moment updates for linear estimators and offers both
linearized and exact-refit nonlinear jackknife paths; choose `jackknife_method =
"refit"` for an exact cross-language comparison.

**Randomization inference and RI-jackknife are absent from v1 source and APIs.**
Their options, generators, exports, and implementation code have been removed;
old calls fail rather than silently selecting another inference method. They
are planned for version 2. This does not remove ordinary cluster jackknife.

These implementations are RCS estimators. They do not implement panel-specific
equations, pre-treatment event-study estimation, anticipation, or all features
of the official `csdid`/`did` packages. The Python and Stata implementations
remain straightforward reference ports; the name does not guarantee a speed
advantage on every problem.

## Examples and validation

See the language READMEs for executable examples and return values.
[VALIDATION.md](VALIDATION.md) describes the checks and limitations.
[benchmarks/README.md](benchmarks/README.md) documents the timing design and
the pinned official Stata `csdid` comparator. Benchmark data are synthetic.

The MERIT example fixture is a state-year collapse from Conley and Taber
(2011), also used in the UNDID paper. Its public source is
<https://economics.uwo.ca/people/conley_docs/code_to_download.html>.
It is included for regression testing; see the language tests for its expected
estimate and jackknife uncertainty.

## References and license

Callaway, B. and P. H. C. Sant'Anna (2021), “Difference-in-Differences with
Multiple Time Periods,” *Journal of Econometrics*, 225(2), 200–230.
<https://doi.org/10.1016/j.jeconom.2020.12.001>.

Sant'Anna, P. H. C. and J. Zhao (2020), “Doubly Robust Difference-in-Differences
Estimators,” *Journal of Econometrics*, 219(1), 101–122.
<https://doi.org/10.1016/j.jeconom.2020.06.003>.

MIT; see [LICENSE](LICENSE). Developed by Matthew D. Webb.
Source lineage is recorded in [PROVENANCE.md](PROVENANCE.md).
