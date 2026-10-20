# Source lineage

The v1 monorepository imports distribution source from these development commits:

| Component | Source repository | Commit |
|---|---|---|
| Stata | mattdwebb/fastcsdidstata | `b1eb2e5f4081390c0c3b6f3d8fca66718913a529` |
| R | mattdwebb/fastrcsdid | `21dd70ae2dbbfa8c52daf73f0673cc052bc47b43` |
| Python | mattdwebb/fastcsdidpy | `810672235c820fcf1125653a406150eeed63a82b` |

The release keeps the established estimator implementations and numerical
fixtures, removes randomization-inference code and APIs, and reorganizes
installation, documentation, tests, and benchmarks. Historical research scripts
and source-repository Git histories are not included. Future work belongs in
this repository; the old development repositories remain archival inputs.

Existing frozen numerical oracle metadata live in
`python/validation/frozen/versions.txt`. Those oracle generation versions precede
the imported source commits; their expected numeric values have been retained.
