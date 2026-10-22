"""Check that distributed source exposes only the supported inference API."""
from pathlib import Path
import re

root = Path(__file__).resolve().parents[1]
patterns = re.compile(
    r"\b(?:rijack|randjack|ri_beta|ri_t_native|ri_t_jk|ri_reps|ri_seed|"
    r"ri_assignments|ri_strata|randomization_inference|prepare_ri_assignments|"
    r"did_assignments|did_placebo|fcs_ri|fcs_check_strata)\b", re.I)
files = (list((root / "R/R").glob("*.R")) +
         list((root / "python/fastcsdidpy").glob("*.py")) +
         list((root / "stata").glob("*.ado")) +
         [root / "R/NAMESPACE", root / "stata/fastcsdid.pkg"])
failures = []
for path in files:
    for number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        if patterns.search(line):
            failures.append(f"{path.relative_to(root)}:{number}: {line}")
for forbidden in ("R/R/ri.R", "R/R/ri-evaluate.R", "R/R/ri-subtract.R",
                  "R/R/ri-covariate.R", "stata/fastcsdid_assignments.ado",
                  "stata/fastcsdid_assignments.sthlp"):
    if (root / forbidden).exists():
        failures.append("Unsupported inference source: " + forbidden)
if failures:
    raise SystemExit("\n".join(failures))
print(f"PASS: {len(files)} runtime/API files conform to the supported inference API.")
