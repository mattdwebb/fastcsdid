"""Validate DR inference by perturbing data weights, not a copied IF formula."""
from pathlib import Path
import numpy as np
import pandas as pd
import pytest
from fastcsdidpy import fit

DATA = Path(__file__).resolve().parents[1] / "validation/frozen/fixture.csv"


@pytest.mark.parametrize("control", ["never", "not_yet"])
@pytest.mark.parametrize("weighted", [False, True])
def test_dr_cluster_influence_is_empirical_weight_derivative(control, weighted):
    data = pd.read_csv(DATA)
    if not weighted:
        data["w"] = 1.0
    options = dict(gvar="g", cluster="cluster", weights="w", covariates=["x", "z"],
                   method="dr", control_group=control, inference="native")
    result = fit(data, "y", "time", **options)
    derivative = []
    step = 1e-5
    for cluster in result["cluster_ids"]:
        mask = (data.cluster == cluster).to_numpy()
        values = []
        for sign in [1, -1]:
            perturbed = data.copy()
            perturbed["w"] = data.w * (1 + sign * step * mask)
            refit = fit(perturbed, "y", "time", **options)
            values.append(np.r_[refit["group_time"].estimate,
                                refit["aggregation"].estimate])
        derivative.append((values[0] - values[1]) / (2 * step))
    derivative = np.asarray(derivative)
    np.testing.assert_allclose(result["cluster_influence"] / len(data), derivative,
                               atol=2e-7, rtol=0)
    np.testing.assert_allclose(result["covariance"], derivative.T @ derivative,
                               atol=2e-8, rtol=0)
    simple = len(result["group_time"])
    assert result["native_se"] == pytest.approx(
        np.linalg.norm(derivative[:, simple]), abs=2e-8)
