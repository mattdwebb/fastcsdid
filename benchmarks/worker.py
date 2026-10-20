from pathlib import Path
import csv
import gc
import sys
from time import perf_counter
import pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "python"))
from fastcsdidpy import fit

source, method, repetitions, output = sys.argv[1:]
d = pd.read_csv(source, float_precision="round_trip")
with open(output, "w", newline="") as stream:
    writer = csv.writer(stream)
    writer.writerow(["engine", "n", "method", "rep", "seconds", "estimate", "native_se"])
    for rep in range(int(repetitions) + 1):
        gc.collect()
        start = perf_counter()
        result = fit(d, "y", "time", gvar="g", cluster="cluster",
                     covariates=None if method == "none" else ["x", "z"],
                     method=method, control_group="never", data_mode="rcs", inference="native")
        elapsed = perf_counter() - start
        if result["status"] != "pass":
            raise RuntimeError(result["status"])
        writer.writerow(["fast_python", len(d), method, rep, elapsed,
                         result["estimate"], result["native_se"]])
        stream.flush()
