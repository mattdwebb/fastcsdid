"""Removed v2 features must fail instead of silently choosing another test."""
import numpy as np
import pandas as pd
import pytest
import fastcsdidpy


def data():
    d = pd.DataFrame([(c, t, j) for c in range(8) for t in range(1, 5)
                      for j in range(3)], columns=["c", "t", "j"])
    d["g"] = np.where(d.c < 4, 3, 0)
    d["y"] = d.c / 10 + d.t + np.sin(d.j + d.c) + ((d.g > 0) & (d.t >= d.g))
    return d


@pytest.mark.parametrize("method", ["ri", "ri_beta", "ri_t_native", "rijack", "randjack", "ri_t_jk", "full"])
def test_removed_inference_rejected(method):
    with pytest.raises(ValueError, match="inference"):
        fastcsdidpy.fit(data(), "y", "t", gvar="g", cluster="c", method="none", inference=method)


@pytest.mark.parametrize("option", ["assignments", "ri_assignments", "reps", "ri_reps",
                                  "seed", "ri_seed", "strata", "ri_strata"])
def test_removed_options_rejected(option):
    with pytest.raises(TypeError, match="Unimplemented"):
        fastcsdidpy.fit(data(), "y", "t", gvar="g", cluster="c", **{option: 1})


def test_removed_export_absent():
    assert not hasattr(fastcsdidpy, "prepare_ri_assignments")
