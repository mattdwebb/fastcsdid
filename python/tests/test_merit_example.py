from pathlib import Path

import pandas as pd
import pytest

from fastcsdidpy import fit


def test_merit_paper_example_keeps_all_defined_scalar_deletions():
    data = pd.read_csv(Path(__file__).resolve().parents[1] / "validation/frozen/merit.csv")
    result = fit(
        data, "outcome", "year", gvar="cohort", cluster="state",
        weights="weight", method="none", inference="jackknife"
    )
    jackknife = result["jackknife"]

    assert result["estimate"] == pytest.approx(0.0477819967004955, abs=1e-12)
    assert jackknife["se"] == pytest.approx(0.0146168254236831, abs=1e-12)
    assert len(jackknife["deleted_estimates"]) == 51
    assert sum(jackknife["support_changed"]) == 4
    assert jackknife["status"] == "pass"
