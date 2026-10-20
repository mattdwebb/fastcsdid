# Native inference timings

Measured on 2026-10-20. Medians of three repetitions after one recorded warm-up.
All times are seconds. CSV files retain individual runs and unrounded values.

| N | Method | csdid 2.0.0 Stata | Fast Stata | Fast R | Fast Python |
|---:|---|---:|---:|---:|---:|
| 3,200 | dr | 0.0770 | 0.0230 | 0.0600 | 0.0348 |
| 3,200 | ipw | 0.0690 | 0.0200 | 0.0500 | 0.0239 |
| 3,200 | none | 0.0710 | 0.0160 | 0.0200 | 0.0183 |
| 3,200 | reg | 0.0830 | 0.0170 | 0.0200 | 0.0246 |
| 32,000 | dr | 0.4510 | 0.2020 | 0.2000 | 0.0793 |
| 32,000 | ipw | 0.3970 | 0.1610 | 0.1800 | 0.0562 |
| 32,000 | none | 0.4830 | 0.1400 | 0.0200 | 0.0522 |
| 32,000 | reg | 0.3820 | 0.1500 | 0.0700 | 0.0537 |
| 320,000 | dr | 4.8920 | 2.0720 | 1.8100 | 0.8591 |
| 320,000 | ipw | 4.3320 | 1.8340 | 1.5400 | 0.6714 |
| 320,000 | none | 3.9790 | 1.4400 | 0.0900 | 0.5674 |
| 320,000 | reg | 4.0730 | 1.6510 | 0.6000 | 0.6799 |

Across all 192 calls, maximum absolute ATT gap: 8.19e-11; maximum analytical SE gap: 4.52e-12.

See ../README.md for scope, comparator revision, reproduction, and timing limitations.
These timings compare analytical-inference workflows. They do not compare cluster-jackknife runtimes.
