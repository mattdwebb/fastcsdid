"""Changing-support CV3J: physical deletion refits, never dropped failures."""
from pathlib import Path
import numpy as np
import pandas as pd
import pytest
from fastcsdidpy import fit
from fastcsdidpy.probability import student_two_sided


@pytest.fixture
def data():
    return pd.read_csv(Path(__file__).resolve().parents[1] / "validation/frozen/fixture.csv")


def fit_data(data, method="none", **options):
    return fit(data, "y", "time", gvar="g", cluster="cluster", method=method,
               covariates=None if method == "none" else ["x", "z"], weights="w", **options)


@pytest.mark.parametrize("method", ["none", "reg", "ipw", "dr"])
@pytest.mark.parametrize("control", ["never", "not_yet"])
def test_singleton_default_is_physical_refit(data, method, control):
    data = data[~data.cluster.isin([2, 3, 4])]
    common = dict(control_group=control, inference="standard")
    result = fit_data(data, method, **common)
    explicit = fit_data(data, method, jackknife_support="changing", **common)
    fixed = fit_data(data, method, jackknife_support="fixed", **common)
    jk = result["jackknife"]
    assert result["metadata"]["jackknife_support"] == "changing"
    assert jk["jackknife_support"] == "changing"
    assert jk["status"] == "pass"
    assert sum(jk["support_changed"]) == 1
    assert fixed["jackknife"]["status"] == "IDENTIFICATION_FAILURE"
    assert len(fixed["jackknife"]["failures"]) == 1
    assert np.isnan(fixed["jackknife"]["p"])
    for field in ("estimate", "native_se", "native_p"):
        assert result[field] == fixed[field] == explicit[field]
    expected = [fit_data(data[data.cluster != h], method, control_group=control,
                         inference="point")["estimate"] for h in result["cluster_ids"]]
    np.testing.assert_allclose(jk["deleted_estimates"], expected, atol=1e-12)
    np.testing.assert_array_equal(jk["deleted_estimates"], explicit["jackknife"]["deleted_estimates"])
    count = len(expected)
    variance = (count - 1) / count * np.sum((expected - np.mean(expected)) ** 2)
    assert jk["variance"] == pytest.approx(variance)
    assert jk["df"] == count - 1
    assert jk["p"] == pytest.approx(student_two_sided(result["estimate"] / np.sqrt(variance), count - 1))
    # Deleting a member of a surviving cohort must change its aggregation mass.
    full = fit_data(data, method, inference="point")
    reduced = fit_data(data[data.cluster != 5], method, inference="point")
    assert not np.allclose(full["aggregation_weights"][0], reduced["aggregation_weights"][0])


@pytest.mark.parametrize("method", ["none", "reg", "ipw", "dr"])
def test_unchanged_support_policies_agree(data, method):
    a = fit_data(data, method, inference="jackknife")
    b = fit_data(data, method, inference="jackknife", jackknife_support="fixed")
    assert not any(a["jackknife"]["support_changed"])
    np.testing.assert_array_equal(a["jackknife"]["deleted_estimates"], b["jackknife"]["deleted_estimates"])
    assert a["jackknife"]["p"] == b["jackknife"]["p"]


@pytest.mark.parametrize("policy", ["changing", "fixed"])
def test_undefined_deletion_not_dropped(data, policy):
    control = data.loc[data.g == 0, "cluster"].iloc[0]
    data = data[(data.g != 0) | (data.cluster == control)]
    result = fit_data(data, inference="jackknife", jackknife_support=policy)
    jk = result["jackknife"]
    assert np.isfinite(result["estimate"])
    assert jk["status"] == "IDENTIFICATION_FAILURE"
    assert len(jk["deleted_estimates"]) == data.cluster.nunique()
    assert jk["df"] == data.cluster.nunique() - 1
    assert any(f["cluster"] == control for f in jk["failures"])
    assert np.isnan(jk["deleted_estimates"][list(result["cluster_ids"]).index(control)])
    assert np.isnan(jk["p"])


@pytest.mark.parametrize("invalid", [None, True, "adaptive", ["changing"]])
def test_invalid_policy_rejected(data, invalid):
    with pytest.raises(ValueError, match="jackknife_support"):
        fit_data(data, inference="point", jackknife_support=invalid)




def test_fixed_support_ignores_already_omitted_cohort(data):
    data.loc[data.cluster == 1, "g"] = 1  # Always treated: omitted by _point.
    a = fit_data(data, inference="jackknife", jackknife_support="fixed")
    assert a["jackknife"]["status"] == "pass"
    assert not any(a["jackknife"]["support_changed"])
