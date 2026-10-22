"""Validate matched ATT/SE values before publishing timing ratios."""
from pathlib import Path
import argparse
import json
import shutil
import numpy as np
import pandas as pd

root = Path(__file__).resolve().parents[1]
p = argparse.ArgumentParser()
p.add_argument("--work", type=Path, default=root / "benchmarks/work")
p.add_argument("--output", type=Path, default=root / "benchmarks/results")
args = p.parse_args()
work = args.work
out = args.output
config = json.loads((work / "config.json").read_text())
data = pd.read_csv(work / "timings.csv")
expected = len(config["sizes"]) * 4 * 4 * (config["reps"] + 1)
assert len(data) == expected, (len(data), expected)
assert not data.duplicated(["engine", "n", "method", "rep"]).any()
assert set(data.engine) == {"csdid_stata", "fast_stata", "fast_R", "fast_python"}
assert np.isfinite(data[["seconds", "estimate", "native_se"]]).all().all()
assert (data.seconds >= 0).all()
baseline = data.query("engine == 'csdid_stata' and rep == 1")[["n", "method", "estimate", "native_se"]]
joined = data.merge(baseline, on=["n", "method"], suffixes=("", "_reference"), validate="many_to_one")
joined["estimate_gap"] = abs(joined.estimate - joined.estimate_reference)
joined["native_se_gap"] = abs(joined.native_se - joined.native_se_reference)
assert joined.estimate_gap.max() < 1e-8, joined.estimate_gap.max()
assert joined.native_se_gap.max() < 1e-8, joined.native_se_gap.max()
measured = joined.query("rep > 0")
summary = measured.groupby(["n", "method", "engine"], as_index=False).agg(
    median_seconds=("seconds", "median"), min_seconds=("seconds", "min"),
    max_seconds=("seconds", "max"), repetitions=("rep", "count"),
    estimate=("estimate", "first"), native_se=("native_se", "first"),
    max_estimate_gap=("estimate_gap", "max"), max_native_se_gap=("native_se_gap", "max"))
reference = summary.query("engine == 'csdid_stata'")[["n", "method", "median_seconds"]].rename(columns={"median_seconds": "csdid_seconds"})
summary = summary.merge(reference, on=["n", "method"], validate="many_to_one")
summary["csdid_over_runtime"] = summary.csdid_seconds / summary.median_seconds
out.mkdir(parents=True, exist_ok=True)
data.to_csv(out / "timings.csv", index=False, float_format="%.17g")
summary.to_csv(out / "summary.csv", index=False, float_format="%.17g")
for name in ("config.json", "data_hashes.json"):
    shutil.copyfile(work / name, out / name)
report = ["# Native inference timings", "",
          "Measured on 2026-10-20. Medians of three repetitions after one recorded warm-up.",
          "All times are seconds. CSV files retain individual runs and unrounded values.", "",
          "| N | Method | csdid 2.0.0 Stata | Fast Stata | Fast R | Fast Python |",
          "|---:|---|---:|---:|---:|---:|"]
for (n, method), group in summary.groupby(["n", "method"], sort=False):
    values = group.set_index("engine").median_seconds
    report.append(f"| {n:,} | {method} | " + " | ".join(f"{values[e]:.4f}" for e in ("csdid_stata", "fast_stata", "fast_R", "fast_python")) + " |")
report += ["", f"Across all {len(joined)} calls, maximum absolute ATT gap: {joined.estimate_gap.max():.3g}; "
           f"maximum analytical SE gap: {joined.native_se_gap.max():.3g}.",
           "", "See the [benchmark guide](../README.md) for scope, comparator revision, reproduction, and timing limitations.",
           "These timings compare analytical-inference workflows. They do not compare cluster-jackknife runtimes."]
(out / "TIMINGS.md").write_text("\n".join(report) + "\n")
print("\n".join(report))
