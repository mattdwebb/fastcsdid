# Stata: fastcsdid 1.0.0

When upgrading a package-managed development installation, first run
`ado uninstall fastcsdid` so obsolete helper files are removed, then install
the new package using the URL in the repository README. If you previously used
`adopath` with a development checkout, switch that path to this directory.

From this directory:

```stata
adopath ++ "."
import delimited using "validation/frozen/fixture.csv", clear
fastcsdid y, time(time) gvar(g) cluster(cluster) covariates(x z) ///
    weights(w) method(dr) control(never) jackknife
matrix list e(attgt)
matrix list e(aggregation)
display e(native_se)
display e(jk_se)
display e(jk_p)
```

`e(V)` always contains the analytical cluster variance, even when `jackknife`
is requested. Jackknife results are stored separately as `e(jk_se)`, `e(jk_p)`,
`e(jk_df)`, and `e(deleted)`. The jackknife uses exact deletion refits,
mean-centered CV3J, and Student t(G−1) inference.

From `validation/`, run `do oracle.do`, then
`python check.py`; also run `do release_api.do` and `do merit.do`.
The command depends only on Stata 16+ and its Mata runtime.
