"""RCS equations ported from fastrcsdid (MIT, Matthew Webb, 2026)."""
from time import perf_counter
import numpy as np
import pandas as pd
from .probability import normal_two_sided


def expit(x):
    return np.exp(-np.logaddexp(0, -x))


class IdentificationError(ValueError):
    """Required comparison or nuisance model is not identified."""


def _solve(A, b):
    if np.linalg.matrix_rank(A) < A.shape[0]:
        raise IdentificationError("singular nuisance design")
    return np.linalg.solve(A, b)


def _ratio(w, z):
    s = w.sum()
    if s <= 0:
        raise IdentificationError("zero weighted cell normalizer")
    a = w @ z / s
    return a, w * (z - a) / w.mean()


def _wls(X, y, w):
    b = _solve(X.T @ (w[:, None] * X), X.T @ (w * y))
    IF = _solve(X.T @ (w[:, None] * X) / len(y),
                (X * (w * (y - X @ b))[:, None]).T).T
    return X @ b, IF


def _logit(X, D, w):
    b = np.zeros(X.shape[1])
    for _ in range(100):
        p = expit(X @ b)
        step = _solve(X.T @ ((w * p * (1-p))[:, None] * X), X.T @ (w * (D-p)))
        loss = np.sum(w * (np.logaddexp(0, X @ b) - D * (X @ b)))
        scale = 1.
        while scale > 2**-25:
            new = b + scale * step
            eta = X @ new
            new_loss = np.sum(w * (np.logaddexp(0, eta) - D * eta))
            # Compare the summed objective at its own floating-point scale.
            # This does not relax the Newton-step convergence threshold.
            roundoff = 64 * np.finfo(float).eps * max(1., abs(loss), abs(new_loss))
            if new_loss <= loss + roundoff:
                break
            scale /= 2
        b = new
        if np.max(np.abs(step)) < 1e-10:
            break
        if np.max(np.abs(b)) > 40:
            raise IdentificationError("propensity separation or extreme coefficients")
    else:
        raise IdentificationError("propensity did not converge")
    p = np.minimum(expit(X @ b), 1-1e-6)
    IF = _solve(X.T @ ((w*p*(1-p))[:, None] * X) / len(D),
                (X * (w*(D-p))[:, None]).T).T
    return p, IF


def _cell(y, P, D, X, w, method):
    w = w / w.mean()
    w10, w11 = w*D*(1-P), w*D*P
    w00, w01 = w*(1-D)*(1-P), w*(1-D)*P
    if method == "none":
        rr = [_ratio(z, y) for z in (w11, w10, w01, w00)]
        return rr[0][0]-rr[1][0]-rr[2][0]+rr[3][0], rr[0][1]-rr[1][1]-rr[2][1]+rr[3][1]
    if method in ("reg", "dr"):
        m00, o00 = _wls(X, y, w00)
        m01, o01 = _wls(X, y, w01)
    if method == "reg":
        a, ai = _ratio(w11, y)
        b, bi = _ratio(w10, y)
        c, ci = _ratio(w*D, m01-m00)
        ex = X.T @ (w*D) / (w*D).sum()
        return a-b-c, ai-bi-ci-(o01-o00) @ ex
    p, pi = _logit(X, D, w)
    trim = (D == 1) | (p < .995)
    odds = p/(1-p)
    w10, w11 = trim*w10, trim*w11
    w00, w01 = trim*odds*w00, trim*odds*w01
    z = y if method == "ipw" else y-P*m01-(1-P)*m00
    r10, i10 = _ratio(w10, z)
    r11, i11 = _ratio(w11, z)
    r00, i00 = _ratio(w00, z)
    r01, i01 = _ratio(w01, z)
    att = r11-r10-r01+r00
    grad = X.T @ (w01*(z-r01))/w01.sum() - X.T @ (w00*(z-r00))/w00.sum()
    inf = i11-i10-i01+i00-pi @ grad
    if method == "dr":
        m10, o10 = _wls(X, y, w*D*(1-P))
        m11, o11 = _wls(X, y, w*D*P)
        wd = trim*w*D
        ep, eip = _ratio(wd, m11-m01)
        e1, ei1 = _ratio(w11, m11-m01)
        eb, eib = _ratio(wd, m10-m00)
        e0, ei0 = _ratio(w10, m10-m00)
        att += ep-e1-eb+e0
        mp = X.T @ (wd/wd.sum()-w11/w11.sum())
        mb = X.T @ (wd/wd.sum()-w10/w10.sum())
        # Include the outcome-model influence terms for both periods.
        inf += o01 @ (X.T @ (w01/w01.sum()-w11/w11.sum()))
        inf += o00 @ (X.T @ (w00/w00.sum()-w10/w10.sum()))
        inf += eip-ei1-eib+ei0+(o11-o01) @ mp-(o10-o00) @ mb
    return att, inf


def _point(d, method, control, window):
    N = len(d)
    times = np.sort(d.t.unique())
    groups = np.sort(d.loc[np.isfinite(d.g), "g"].unique())
    y, t, g, w = [d[k].to_numpy() for k in ("y", "t", "g", "w")]
    X = np.column_stack([np.ones(N)] + [d[k].to_numpy() for k in d if k.startswith("x:")])
    records, cols = [], []
    omitted = []
    for group in groups:
        pre = times[times < group]
        post = times[times >= group]
        if not len(pre) or not len(post):
            omitted.append(float(group))
            continue
        for period in post:
            event = period-group
            if window is not None and not window[0] <= event <= window[1]:
                continue
            base = pre[-1]
            controls = np.isinf(g) if control == "never" else g > period
            rows = ((g == group) | controls) & ((t == base) | (t == period))
            inf = np.full(N, np.nan)
            try:
                a, local = _cell(y[rows], (t[rows] == period).astype(float),
                                 (g[rows] == group).astype(float), X[rows], w[rows], method)
                inf[:] = 0
                inf[rows] = N/rows.sum()*local
                status, reason = "pass", ""
            except (IdentificationError, np.linalg.LinAlgError) as exc:
                a, status, reason = np.nan, "IDENTIFICATION_FAILURE", str(exc)
            records.append(dict(group=group, base=base, post=period, event=event,
                                estimate=a, status=status, reason=reason))
            cols.append(inf)
    if not records:
        raise IdentificationError("No valid post-treatment comparisons in requested window")
    cells = pd.DataFrame(records)
    I = np.column_stack(cols)
    mass = np.array([w[g == z].sum() for z in cells.group])
    aggregates, agif, weights = [], [], []
    selectors = [("simple", np.nan, np.ones(len(cells), dtype=bool))]
    for kind, col in (("group", "group"), ("calendar", "post"), ("dynamic", "event")):
        selectors.extend((kind, z, (cells[col] == z).to_numpy()) for z in sorted(cells[col].unique()))
    for kind, index, sel in selectors:
        raw = np.where(sel, 1 if kind == "group" else mass, 0.)
        active = raw > 0
        valid = active.any() and np.isfinite(cells.estimate[active]).all()
        W = raw/raw.sum() if valid else np.full(len(cells), np.nan)
        estimate = W[active] @ cells.estimate[active] if valid else np.nan
        influence = I[:, active] @ W[active] if valid else np.full(N, np.nan)
        if valid and kind != "group":
            for cohort in np.unique(cells.group[active]):
                use = active & (cells.group.to_numpy() == cohort)
                influence += N*w*(g == cohort)*(cells.estimate[use]-estimate).sum()/raw.sum()
        label = "simple" if kind == "simple" else f"{'event' if kind == 'dynamic' else kind}:{index:g}"
        aggregates.append(dict(label=label, type=kind, index=index, estimate=estimate))
        agif.append(influence)
        weights.append(W)
    A = np.column_stack(agif)
    joint = np.column_stack([I, A])
    ids = pd.unique(d.c)
    cluster_if = np.stack([joint[d.c.to_numpy() == h].sum(axis=0) for h in ids])
    cov = cluster_if.T @ cluster_if/N**2
    agg = pd.DataFrame(aggregates)
    j = len(cells)
    se = np.sqrt(cov[j, j])
    theta = agg.estimate.iloc[0]
    return dict(estimate=float(theta), att=float(theta), target_label="cohort-mass weighted group-time ATT",
                status="pass" if np.isfinite(theta) else "IDENTIFICATION_FAILURE",
                group_time=cells, aggregation=agg, aggregation_weights=np.stack(weights),
                influence_att=I, influence_aggregate=A, cluster_influence=cluster_if,
                covariance=cov, cluster_ids=ids, native_se=float(se),
                native_p=normal_two_sided(theta/se) if se > 0 else np.nan,
                native_distribution="normal", nobs=N, n_clusters=len(ids),
                diagnostics={"omitted_groups": omitted})


def fit(data, y, time, unit=None, gvar=None, treatment=None, cluster=None,
        covariates=None, event_window=None, data_mode="auto", **options):
    """Fit the stationary RCS estimator; panel is explicitly unsupported.

    Options: method/est_method, control_group, weights, never_treated,
    inference (native, jackknife, standard, point, or none), null,
    jackknife_support ("changing" by default; "fixed" rejects changed ATT support).
    Covariates are numeric column names, observed in each RCS period.
    """
    start = perf_counter()
    allowed = {"method", "est_method", "control_group", "weights", "never_treated", "inference",
               "null",
               "anticipation", "jackknife_method", "jackknife_support"}
    unknown = set(options)-allowed
    if unknown:
        raise TypeError(f"Unimplemented options: {sorted(unknown)}")
    if options.get("anticipation", 0) != 0:
        raise NotImplementedError("fastcsdidpy supports anticipation=0 only")
    if options.get("jackknife_method", "refit") not in ("refit", "exact"):
        raise NotImplementedError("only complete exact/refit jackknife is implemented")
    support = options.get("jackknife_support", "changing")
    if not isinstance(support, str) or support not in ("changing", "fixed"):
        raise ValueError("jackknife_support must be changing or fixed")
    if data_mode not in ("auto", "rcs", "panel"):
        raise ValueError("data_mode must be auto, rcs, or panel")
    frame = pd.DataFrame(data)
    mode = ("panel" if unit is not None and frame[unit].duplicated().any() else "rcs") if data_mode == "auto" else data_mode
    if mode == "panel":
        raise NotImplementedError("fastcsdidpy supports repeated cross sections; panel data are unsupported")
    if cluster is None:
        raise ValueError("cluster is required and must have a constant treatment cohort")
    if gvar is None and treatment is None:
        raise ValueError("provide gvar or absorbing treatment")
    method = options.get("method", options.get("est_method", "dripw"))
    method = {"dripw": "dr", "ra": "reg", "stdipw": "ipw"}.get(method, method)
    if method not in ("none", "reg", "ipw", "dr"):
        raise ValueError("method must be none, reg/ra, ipw/stdipw, or dr/dripw")
    covariates = [covariates] if isinstance(covariates, str) else list(covariates or [])
    if covariates and method == "none":
        raise ValueError("method=none cannot silently ignore covariates")
    control = options.get("control_group", "never")
    control = {"nevertreated": "never", "notyettreated": "not_yet"}.get(control, control)
    if control not in ("never", "not_yet"):
        raise ValueError("control_group must be never or not_yet")
    d = pd.DataFrame({"y": frame[y].to_numpy(dtype=float), "t": frame[time].to_numpy(dtype=float),
                      "c": frame[cluster].to_numpy()})
    wc = options.get("weights")
    d["w"] = frame[wc].to_numpy(dtype=float) if wc is not None else 1.
    for col in covariates:
        d[f"x:{col}"] = frame[col].to_numpy(dtype=float)
    if d.c.isna().any() or not np.isfinite(d.drop(columns="c").to_numpy()).all() or (d.w < 0).any() or not len(d):
        raise ValueError("outcome, time, covariates, cluster and nonnegative weights must be observed")
    if gvar is not None:
        d["g"] = frame[gvar].to_numpy(dtype=float)
        d.loc[d.g == options.get("never_treated", 0), "g"] = np.inf
        if d.g.isna().any() or np.isneginf(d.g).any() or (d.groupby("c").g.nunique() != 1).any():
            raise ValueError("cohort must be observed and constant within cluster")
    if treatment is not None:
        D = frame[treatment].to_numpy()
        if not np.isin(D, [0, 1]).all():
            raise ValueError("treatment must be binary")
        treated_time = d.t.where(D == 1)
        inferred = treated_time.groupby(d.c, sort=False).min().fillna(np.inf)
        gi = d.c.map(inferred).to_numpy()
        expected = (d.t.to_numpy() >= gi).astype(int)
        if not np.array_equal(expected, D):
            raise ValueError("treatment must be absorbing and constant within cluster-period")
        if gvar is None:
            d["g"] = gi
        elif not np.array_equal((d.t >= d.g).astype(int).to_numpy(), D):
            raise ValueError("treatment and cohort are inconsistent")
    if event_window is not None:
        if len(event_window) != 2 or not np.isfinite(event_window).all() or event_window[0] > event_window[1]:
            raise ValueError("event_window must be finite ordered (lower, upper)")
        if event_window[0] < 0:
            raise NotImplementedError("pre-treatment event effects are unsupported; event_window must start at zero or later")
    result = _point(d, method, control, event_window)
    result["metadata"] = dict(estimator="CSDID", method=method, common_method="dripw" if method == "dr" else method,
                              jackknife_support=support,
                              data_mode=mode, control_group=control, event_window=event_window,
                              covariates_requested=covariates, covariates_used=covariates,
                              covariate_timing="contemporaneous stationary RCS; no baseline freezing",
                              aggregation_weighting="total sampling weight per cohort across all periods",
                              target_compatible=False, native_vce="cluster IF crossproduct / N^2; no correction")
    null = float(options.get("null", 0))
    if result["native_se"] > 0:
        result["native_p"] = normal_two_sided((result["estimate"]-null)/result["native_se"])
    from .inference import add_inference
    add_inference(result, d, method, control, event_window, options)
    result["runtime"] = perf_counter()-start
    return result


def influence_function(result, type="att", cluster_sum=False):
    if type not in ("att", "aggregate"):
        raise ValueError("type must be att or aggregate")
    k = len(result["group_time"])
    if cluster_sum:
        return result["cluster_influence"][:, :k if type == "att" else None] if type == "att" else result["cluster_influence"][:, k:]
    return result["influence_att" if type == "att" else "influence_aggregate"].copy()


def rif(result, type="att"):
    values = result["group_time" if type == "att" else "aggregation"].estimate.to_numpy()
    return influence_function(result, type)+values
