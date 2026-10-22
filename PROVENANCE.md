# Source lineage

The implementations originate from the following repositories and revisions:

| Component | Source repository | Commit |
|---|---|---|
| Stata | mattdwebb/fastcsdidstata | `b1eb2e5f4081390c0c3b6f3d8fca66718913a529` |
| R | mattdwebb/fastrcsdid | `21dd70ae2dbbfa8c52daf73f0673cc052bc47b43` |
| Python | mattdwebb/fastcsdidpy | `810672235c820fcf1125653a406150eeed63a82b` |

The frozen numerical reference outputs were generated using R, `did`, and
`DRDID`. Their software versions and the R source revision are recorded in
[versions.txt](python/validation/frozen/versions.txt). These describe the
reference calculations, rather than the current package dependencies.
