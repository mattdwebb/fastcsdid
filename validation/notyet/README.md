# Not-yet-treated DRIPW validation

This synthetic repeated-cross-section fixture has 1,152 observations, twelve
clusters, two treated cohorts, and two covariates. All cohort-time comparisons use
never-treated or not-yet-treated controls as eligible at the comparison date.

`compare.R` compares each cell estimate and cluster-summed influence function
with CRAN DRDID 1.4.1, then independently aggregates the cells and estimated
cohort weights. `compare.do` checks the aggregate against official Stata csdid
2.0.3. `compare.py` reconciles all five outputs. The maximum estimate gap is
6.77e-12 and maximum native-SE gap is 3.06e-9 (tolerance 1e-8).

Run the R script with fastrcsdid and DRDID >=1.4.0 installed in a test library.
Set `CSDID_OFFICIAL_ADO` to the frozen official comparator directory, run the
Stata do-file from this directory, then run the Python script with fastcsdidpy
installed. `provenance.json` records the data seed, source hashes and comparator
commit. This checks the documented stationary RCS method, not all CSDID options.

`R_cells.csv` also records the maximum absolute and root-mean-square differences
between the two cluster-summed influence vectors for each cell. The largest
cell-level absolute difference is 2.59e-8, below the 1e-7 cell-influence tolerance
in `compare.R`. This differs from the 1e-8 aggregate estimate/SE tolerance above.
The fixture SHA-256 is
`76791831e0510d3492f94e1d4a62c60b25edbc16bfdf9376dae66675401923c8`.

To regenerate the additional discrepancy columns in an R session in this folder:

```r
source("compare.R")
gaps <- cell - reference
write.csv(data.frame(map, fast_att = unname(fit$att_gt), upstream_att = values,
  max_abs_cluster_if_gap = apply(abs(gaps), 2, max),
  rms_cluster_if_gap = sqrt(colMeans(gaps^2))),
  "R_cells.csv", row.names = FALSE)
```
