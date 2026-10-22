# Timing fastcsdid against Stata csdid 2.0.0

The comparator is **csdid 2.0.0**, distributed 2026-09-27, from
<https://github.com/pedrohcgs/csdid-stata>, pinned at
`31a577382d7c09757ed3baa095659f45cc63c107`.
The benchmark loads this version from a separate checkout.

Read [results/TIMINGS.md](results/TIMINGS.md) for measured results.
These results were measured with fastcsdid v1.0.0; `results/config.json`
records the SHA-256 hashes of the source files used.

## Matched workload

- Synthetic stationary RCS with 40 clusters, 8 periods, and 10 clusters in each
  of cohorts 3, 5, 7, and never treated.
- 3,200, 32,000, and 320,000 observations; no sampling weights.
- Unadjusted (`none`, mapped to intercept-only `reg` in official csdid),
  regression, IPW, and doubly robust methods. The latter three use two
  covariates, `x` and `z`.
- Never-treated controls, last pre-adoption base, and all post-treatment
  comparisons. Simple, group, calendar, and event aggregations are requested.
- Analytical cluster influence-function uncertainty; no bootstrap or
  jackknife. R `inference="none"` still computes analytical covariance.
- csdid specifies `rcs nevertreated base_period(universal) analytical pointwise`
  explicitly. This avoids its different comparison-group and bootstrap defaults.

csdid additionally computes pre-treatment cells and a pre-test, which the fast
implementations do not provide. Consequently this is a comparison of complete
supported command workflows with matching post-treatment estimands, not an
isolated identical kernel microbenchmark. Event outputs differ in pre-treatment
coverage. The reported numerical agreement checks the simple ATT and its native
standard error; separate regression tests cover post-treatment component effects.

## Timing protocol

Each engine/workload runs serially, with one warm-up followed by three measured
calls. Input reading, process startup, package/source loading, and output writing
are excluded. Raw-data preparation inside each estimator is included. BLAS and
OpenMP thread counts are set to one; Stata is SE, with `set processors 1` used
where supported. R garbage collection is performed before each call. Python GC
is likewise collected before timing.

Warm-ups are retained as repetition zero in `timings.csv` and excluded from
medians. `summary.csv` includes minima, maxima, and
`csdid_over_runtime`: values above one mean the fast implementation was faster.
The summary refuses to publish ratios if any matched ATT or SE differs by
1e-8 or more. CPU, memory, OS, runtime versions, source hashes, and data hashes
are recorded alongside results. Other machine activity can affect wall time;
these results are specific to the recorded machine and workloads.

## Reproduce

Use Python with NumPy/pandas, R, a licensed Stata installation, and the pinned
comparator checkout. Run from the repository root:

```sh
git clone https://github.com/pedrohcgs/csdid-stata.git ../csdid-stata
git -C ../csdid-stata checkout 31a577382d7c09757ed3baa095659f45cc63c107
python benchmarks/run.py --stata "/path/to/stata" --rscript "/path/to/Rscript" --csdid ../csdid-stata
python tools/summarize_benchmarks.py
```

The supplied Stata executable must accept batch arguments `/e do` (the measured
run uses Windows Stata SE). Adapt the launcher for a non-Windows Stata binary.
No global Stata installation changes are needed. Completed workload checkpoints
are reused only with the same source hashes and configuration. Change `--work`
for a fresh timing run. Generated data and raw Stata logs are excluded from
version control; published results include timings, source/data hashes, and
machine specifications.
