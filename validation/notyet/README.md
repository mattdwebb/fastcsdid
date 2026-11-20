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
