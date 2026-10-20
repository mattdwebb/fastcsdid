"""Complete-refit cluster jackknife (mean-centered CV3J)."""
import numpy as np
from .probability import student_two_sided
from .core import _point, IdentificationError


def _jackknife(d, method, control, window, full, null, support="changing"):
    ids = full["cluster_ids"]
    if len(ids) < 2:
        return dict(status="IDENTIFICATION_FAILURE", reason="at least two clusters required",
                    jackknife_support=support)
    estimates, failures, support_changed = [], [], []
    columns = ["group", "base", "post", "event"]
    for h in ids:
        changed = None
        try:
            deleted = _point(d.loc[d.c != h], method, control, window)
            changed = not np.array_equal(deleted["group_time"][columns].to_numpy(), full["group_time"][columns].to_numpy())
            if support == "fixed" and changed:
                raise IdentificationError("deletion changed required ATT support under fixed policy")
            if deleted["status"] != "pass" or not np.isfinite(deleted["estimate"]):
                raise IdentificationError("deletion aggregate is unidentified")
            estimates.append(deleted["estimate"])
        except (ValueError, np.linalg.LinAlgError) as exc:
            estimates.append(np.nan)
            failures.append(dict(cluster=h, reason=str(exc)))
        support_changed.append(changed)
    a = np.array(estimates)
    diagnostics = dict(jackknife_support=support, support_changed=support_changed,
                       delete_clusters=ids.copy())
    if failures:
        return dict(status="IDENTIFICATION_FAILURE", reason="every deletion is required",
                    failures=failures, deleted_estimates=a, se=np.nan, p=np.nan, df=len(ids)-1,
                    **diagnostics)
    v = (len(ids)-1)/len(ids)*np.sum((a-a.mean())**2)
    se = np.sqrt(v)
    stat = (full["estimate"]-null)/se if se > 0 else np.nan
    status = "pass" if np.isfinite(stat) else ("ZERO_VARIANCE" if se == 0 else "IDENTIFICATION_FAILURE")
    return dict(status=status, deleted_estimates=a,
                variance=float(v), se=float(se), statistic=float(stat), df=len(ids)-1,
                p=student_two_sided(stat, len(ids)-1) if np.isfinite(stat) else np.nan,
                method="complete exact refit CV3J", **diagnostics)


def add_inference(result, d, method, control, window, options):
    requested = options.get("inference", "native")
    bundles = {"none": [], "point": [], "standard": ["native", "jackknife"]}
    methods = bundles.get(requested, [requested]) if isinstance(requested, str) else list(requested)
    if set(methods) - {"native", "jackknife"}:
        raise ValueError("inference must be native, jackknife, standard, point, or none")
    result["inference"] = {}
    if "jackknife" in methods:
        result["jackknife"] = _jackknife(
            d, method, control, window, result, float(options.get("null", 0)),
            options.get("jackknife_support", "changing"))
        result["inference"]["jackknife"] = result["jackknife"]
