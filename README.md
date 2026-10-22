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

Each implementation runs independently; install the package for your preferred
language.

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

Analytical standard errors use cluster-summed influence functions. The cluster
jackknife deletes one cluster at a time. Stata and Python completely refit each
deletion. R uses exact moment updates for linear estimators and offers both
linearized and exact-refit nonlinear jackknife paths; choose
`jackknife_method = "refit"` for an exact cross-language comparison.

These implementations are RCS estimators. They do not implement panel-specific
equations, pre-treatment event-study estimation, anticipation, or all features
of the official `csdid`/`did` packages. See the
[timing results](benchmarks/results/TIMINGS.md) for performance across languages,
methods, and sample sizes.

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

Karim, S. R., M. Ø. Nielsen, J. G. MacKinnon, and M. D. Webb (2026),
“Improved Inference for CSDID Using the Cluster Jackknife,” arXiv:2602.12043.
<https://arxiv.org/abs/2602.12043>.

Sant'Anna, P. H. C. and J. Zhao (2020), “Doubly Robust Difference-in-Differences
Estimators,” *Journal of Econometrics*, 219(1), 101–122.
<https://doi.org/10.1016/j.jeconom.2020.06.003>.

MIT; see [LICENSE](LICENSE). Developed by Matthew D. Webb.
Source lineage is recorded in [PROVENANCE.md](PROVENANCE.md).
