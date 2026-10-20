"""Serial, checkpointed comparisons against pinned csdid-stata 2.0.0.

Run from any directory. Supply a licensed Stata executable, Rscript, and an
unmodified checkout of the official comparator. Never changes installed ado files.
"""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import os
import platform
import subprocess
import sys
import time

for variable in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS",
                 "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[variable] = "1"

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
PIN = "31a577382d7c09757ed3baa095659f45cc63c107"
FIELDS = ["engine", "n", "method", "rep", "seconds", "estimate", "native_se"]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run_process(command, cwd, log):
    kwargs = {}
    if os.name == "nt":
        si = subprocess.STARTUPINFO()
        si.dwFlags |= subprocess.STARTF_USESHOWWINDOW
        si.wShowWindow = 0
        kwargs["startupinfo"] = si
        kwargs["creationflags"] = subprocess.CREATE_NO_WINDOW
    with log.open("w") as output:
        subprocess.run(command, cwd=cwd, stdout=output, stderr=subprocess.STDOUT,
                       check=True, **kwargs)


def generate(n, work):
    if n % 320:
        raise ValueError("N must be a multiple of 40 clusters x 8 periods")
    rng = np.random.default_rng(20261020 + n)
    c = np.repeat(np.arange(1, 41), n // 40)
    t = np.tile(np.repeat(np.arange(1, 9), n // 320), 40)
    g = np.array([3] * 10 + [5] * 10 + [7] * 10 + [0] * 10)[c - 1]
    x = rng.normal(size=n) + .15 * np.sin(c)
    z = rng.binomial(1, .4, size=n).astype(float)
    treated = (g > 0) & (t >= g)
    y = .2 * c + .3 * t + .4 * x - .2 * z + treated * (1 + .1 * (t - g)) + rng.normal(size=n)
    d = pd.DataFrame(dict(y=y, x=x, z=z, cluster=c, time=t, g=g))
    d.to_csv(work / f"data_{n}.csv", index=False, float_format="%.17g")
    d.to_stata(work / f"data_{n}.dta", write_index=False, version=118)


def stata_job(engine, n, method, args, work):
    key = f"{engine}_{n}_{method}"
    output = work / f"{key}.csv"
    cov = "" if method == "none" else "x z"
    lines = [
        "version 16", "clear all", "set more off", "set type double",
        "capture set processors 1",
        f'adopath ++ "{(args.csdid / "pkg").as_posix()}"',
        f'adopath ++ "{(ROOT / "stata").as_posix()}"',
        f'use "{(work / f"data_{n}.dta").as_posix()}", clear',
        "tempname result", f'postfile `result\' str20 engine long n str8 method int rep double seconds double estimate double native_se using "{key}.dta", replace',
    ]
    if engine == "csdid_stata":
        lines += ["which csdid", "csdid version"]
    lines += [f"forvalues rep=0/{args.reps} {{", "timer clear 1", "timer on 1"]
    if engine == "csdid_stata":
        mapped = "reg" if method == "none" else method
        lines += [
            f"quietly csdid y {cov}, time(time) gvar(g) cluster(cluster) rcs nevertreated base_period(universal) method({mapped}) analytical pointwise",
            "quietly estat group", "quietly estat calendar", "quietly estat event",
            "quietly estat simple, post",
            "scalar bench_att = _b[ATT]", "scalar bench_se = _se[ATT]",
        ]
    else:
        covopt = f"covariates({cov})" if cov else ""
        lines += [
            f"quietly fastcsdid y, time(time) gvar(g) cluster(cluster) method({method}) control(never) {covopt}",
            "scalar bench_att = _b[ATT]", "scalar bench_se = e(native_se)",
        ]
    lines += ["timer off 1", "quietly timer list 1",
              f'post `result\' ("{engine}") ({n}) ("{method}") (`rep\') (r(t1)) (bench_att) (bench_se)',
              "}", "postclose `result'", f'use "{key}.dta", clear',
              "format seconds estimate native_se %24.17g",
              f'export delimited using "{output.as_posix()}", replace datafmt',
              'display "FASTCSDID_BENCHMARK_COMPLETE"', "exit, clear"]
    dofile = work / f"{key}.do"
    dofile.write_text("\n".join(lines) + "\n")
    run_process([str(args.stata), "/e", "do", dofile.name], work, work / f"{key}.stdout")
    log = work / f"{key}.log"
    # Stata batch can exit zero even when a do-file failed.
    if not output.exists() or "FASTCSDID_BENCHMARK_COMPLETE\n" not in log.read_text(errors="replace"):
        raise RuntimeError(f"Stata benchmark failed; inspect {log}")
    return output


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--stata", type=Path, required=True)
    p.add_argument("--rscript", type=Path, required=True)
    p.add_argument("--csdid", type=Path, required=True)
    p.add_argument("--sizes", type=int, nargs="+", default=[3200, 32000, 320000])
    p.add_argument("--reps", type=int, default=3)
    p.add_argument("--work", type=Path, default=ROOT / "benchmarks/work")
    args = p.parse_args()
    args.csdid = args.csdid.resolve()
    work = args.work.resolve()
    work.mkdir(parents=True, exist_ok=True)
    head = subprocess.check_output(["git", "-C", str(args.csdid), "rev-parse", "HEAD"], text=True).strip()
    dirty = subprocess.check_output(["git", "-C", str(args.csdid), "status", "--porcelain"], text=True)
    if head != PIN or dirty:
        raise RuntimeError(f"Comparator must be an unmodified checkout of {PIN}")
    source_hashes = {p.relative_to(ROOT).as_posix(): digest(p)
                     for directory, suffix in (("python/fastcsdidpy", ".py"), ("R/R", ".R"), ("stata", ".ado"), ("benchmarks", ".py"), ("benchmarks", ".R"))
                     for p in (ROOT / directory).glob("*" + suffix)}
    config = {"csdid_commit": head, "sizes": args.sizes, "reps": args.reps,
              "source_hashes": source_hashes, "python": sys.version,
              "numpy": np.__version__, "pandas": pd.__version__, "os": platform.platform(),
              "processor": platform.processor(), "threads": 1,
              "design": "stationary unweighted RCS; 40 clusters; 8 periods; cohorts 3/5/7/never; 10 clusters each",
              "timing": "loaded data; one untimed warm-up recorded as rep 0; all four aggregates; native cluster IF inference",
              "extra_comparator_work": "csdid also estimates pre-treatment cells and computes a pre-test"}
    config_path = work / "config.json"
    if config_path.exists() and json.loads(config_path.read_text()) != config:
        raise RuntimeError("Benchmark inputs changed. Use a new --work directory.")
    config_path.write_text(json.dumps(config, indent=2))
    rows = []
    for n in args.sizes:
        if not (work / f"data_{n}.dta").exists():
            generate(n, work)
        for method in ("none", "reg", "ipw", "dr"):
            for engine in ("csdid_stata", "fast_stata", "fast_R", "fast_python"):
                key = f"{engine}_{n}_{method}"
                checkpoint = work / f"{key}.complete.json"
                output = work / f"{key}.csv"
                if checkpoint.exists() and output.exists() and json.loads(checkpoint.read_text())["sha256"] == digest(output):
                    print("RESUME", key, flush=True)
                else:
                    print("RUN", key, flush=True)
                    if engine.endswith("stata"):
                        output = stata_job(engine, n, method, args, work)
                    else:
                        command = ([str(args.rscript), str(ROOT / "benchmarks/worker.R")] if engine == "fast_R"
                                   else [sys.executable, str(ROOT / "benchmarks/worker.py")])
                        command += [str(work / f"data_{n}.csv"), method, str(args.reps), str(output)]
                        run_process(command, ROOT, work / f"{key}.stdout")
                    result = pd.read_csv(output)
                    if list(result.columns) != FIELDS or len(result) != args.reps + 1 or not np.isfinite(result[["seconds", "estimate", "native_se"]]).all().all():
                        raise RuntimeError(f"Invalid output {output}")
                    checkpoint.write_text(json.dumps({"sha256": digest(output)}))
                rows.extend(pd.read_csv(output).to_dict("records"))
                pd.DataFrame(rows).to_csv(work / "timings.csv", index=False, float_format="%.17g")
                print("DONE", key, flush=True)
    inputs = {p.name: digest(p) for p in work.glob("data_*") if p.suffix in (".csv", ".dta")}
    (work / "data_hashes.json").write_text(json.dumps(inputs, indent=2))
    print("BENCHMARK_RUN_COMPLETE", flush=True)


if __name__ == "__main__":
    main()
