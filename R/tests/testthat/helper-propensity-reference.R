reference_nonlinear_prepare <- function(data, xformla, weights = NULL,
                                        never_treated = 0) {
  model_frame <- stats::model.frame(
    xformla, data = data, na.action = stats::na.pass,
    drop.unused.levels = FALSE
  )
  core <- data[c("y", "time", "cohort", "region")]
  if (!is.null(weights)) core[[weights]] <- data[[weights]]
  model_ok <- stats::complete.cases(model_frame)
  for (column in model_frame) {
    if (is.matrix(column) && is.numeric(column)) {
      model_ok <- model_ok & apply(is.finite(column), 1L, all)
    } else if (is.numeric(column) || is.integer(column)) {
      model_ok <- model_ok & is.finite(column)
    }
  }
  keep <- stats::complete.cases(core) & model_ok
  for (column in core) {
    if (is.numeric(column) || is.integer(column)) keep <- keep & is.finite(column)
  }
  data <- droplevels(data[keep, , drop = FALSE])
  model_frame <- droplevels(model_frame[keep, , drop = FALSE])
  X <- stats::model.matrix(attr(model_frame, "terms"), model_frame)
  cohort <- numeric_labels(data$cohort)
  cohort[cohort == numeric_labels(never_treated)] <- Inf
  list(
    data = data, X = X, y = numeric_labels(data$y),
    time = numeric_labels(data$time), cohort = cohort,
    region = data$region,
    w = if (is.null(weights)) rep(1, nrow(data)) else
      numeric_labels(data[[weights]])
  )
}

reference_att_sample <- function(prepared, cell,
                                 control = c("never", "not_yet")) {
  control <- match.arg(control)
  treated <- prepared$cohort == cell$group
  controls <- if (control == "never") is.infinite(prepared$cohort) else
    is.infinite(prepared$cohort) | prepared$cohort > cell$post
  rows <- which((treated | controls) &
                  prepared$time %in% c(cell$base, cell$post))
  list(
    rows = rows,
    y = prepared$y[rows],
    X = prepared$X[rows, , drop = FALSE],
    w = prepared$w[rows],
    D = as.numeric(treated[rows]),
    post = as.numeric(prepared$time[rows] == cell$post)
  )
}

reference_logit <- function(X, D, w, start = NULL) {
  fit <- suppressWarnings(stats::glm.fit(
    x = X, y = D, weights = w, start = start,
    family = stats::binomial(), intercept = FALSE,
    control = stats::glm.control(epsilon = 1e-11, maxit = 200L)
  ))
  if (!isTRUE(fit$converged) || anyNA(fit$coefficients)) return(NULL)
  list(coefficients = unname(fit$coefficients),
       probabilities = unname(fit$fitted.values), iterations = fit$iter)
}

reference_ratio <- function(a, z) {
  denominator <- sum(a)
  if (!is.finite(denominator) || denominator <= 0) return(NA_real_)
  sum(a * z) / denominator
}

reference_nonlinear_att <- function(sample, method = c("ipw", "dr")) {
  method <- match.arg(method)
  overlap <- reference_logit(sample$X, sample$D, rep(1, length(sample$D)))
  if (is.null(overlap) || max(overlap$probabilities) >= 0.999) return(NA_real_)
  ps <- reference_logit(sample$X, sample$D, sample$w)
  if (is.null(ps)) return(NA_real_)
  p <- pmin(ps$probabilities, 1 - 1e-6)
  trim <- sample$D == 1 | p < 0.995
  odds <- p / (1 - p)
  a10 <- trim * sample$w * sample$D * (1 - sample$post)
  a11 <- trim * sample$w * sample$D * sample$post
  a00 <- trim * sample$w * odds * (1 - sample$D) * (1 - sample$post)
  a01 <- trim * sample$w * odds * (1 - sample$D) * sample$post
  if (method == "ipw") {
    return((reference_ratio(a11, sample$y) - reference_ratio(a10, sample$y)) -
      (reference_ratio(a01, sample$y) - reference_ratio(a00, sample$y)))
  }

  fit_wls <- function(d, s) {
    rows <- sample$D == d & sample$post == s & sample$w > 0
    if (sum(rows) < ncol(sample$X)) return(NULL)
    fit <- stats::lm.wfit(
      sample$X[rows, , drop = FALSE], sample$y[rows], sample$w[rows]
    )
    if (fit$rank != ncol(sample$X) || anyNA(fit$coefficients)) NULL else
      unname(fit$coefficients)
  }
  beta00 <- fit_wls(0, 0); beta01 <- fit_wls(0, 1)
  beta10 <- fit_wls(1, 0); beta11 <- fit_wls(1, 1)
  if (any(vapply(list(beta00, beta01, beta10, beta11), is.null, logical(1)))) {
    return(NA_real_)
  }
  m00 <- as.numeric(sample$X %*% beta00)
  m01 <- as.numeric(sample$X %*% beta01)
  m10 <- as.numeric(sample$X %*% beta10)
  m11 <- as.numeric(sample$X %*% beta11)
  m0 <- sample$post * m01 + (1 - sample$post) * m00
  residual_part <-
    (reference_ratio(a11, sample$y - m0) -
       reference_ratio(a10, sample$y - m0)) -
    (reference_ratio(a01, sample$y - m0) -
       reference_ratio(a00, sample$y - m0))
  ad <- trim * sample$w * sample$D
  correction <-
    reference_ratio(ad, m11 - m01) - reference_ratio(a11, m11 - m01) -
    reference_ratio(ad, m10 - m00) + reference_ratio(a10, m10 - m00)
  residual_part + correction
}

reference_nonlinear_fit <- function(data, xformla, method = c("ipw", "dr"),
                                    control = c("never", "not_yet"),
                                    weights = NULL,
                                    att_map = reference_reg_map(data)) {
  method <- match.arg(method)
  control <- match.arg(control)
  prepared <- reference_nonlinear_prepare(data, xformla, weights)
  att <- vapply(seq_len(nrow(att_map)), function(k) {
    reference_nonlinear_att(
      reference_att_sample(prepared, att_map[k, ], control), method
    )
  }, numeric(1))
  names(att) <- paste0("ATT(", att_map$group, ",", att_map$post, ")")
  sizes <- reference_group_size(prepared$data, sort(unique(att_map$group)),
                                weights)
  list(
    att = att,
    aggregate = reference_aggregate(att, att_map, sizes),
    att_map = att_map
  )
}

did_nonlinear_reference <- function(data, xformla, method, control,
                                    weighted = FALSE,
                                    compute_inffunc = FALSE,
                                    clustered = FALSE) {
  did::att_gt(
    yname = "y", tname = "time", idname = "id", gname = "cohort",
    xformla = xformla, data = data, panel = FALSE,
    control_group = if (control == "never") "nevertreated" else
      "notyettreated",
    anticipation = 0, weightsname = if (weighted) "w" else NULL,
    clustervars = if (clustered) "region" else NULL,
    bstrap = FALSE, cband = FALSE, est_method = method,
    base_period = "universal", faster_mode = TRUE,
    compute_inffunc = compute_inffunc
  )
}

reference_did_aggregates <- function(reference) {
  simple <- did::aggte(reference, type = "simple", bstrap = FALSE, cband = FALSE)
  out <- c(simple = simple$overall.att)
  for (type in c("group", "calendar", "dynamic")) {
    result <- did::aggte(reference, type = type, bstrap = FALSE, cband = FALSE)
    keep <- if (type == "dynamic") result$egt >= 0 else
      rep(TRUE, length(result$egt))
    prefix <- switch(type, group = "group:", calendar = "calendar:",
                     dynamic = "event:")
    values <- result$att.egt[keep]
    names(values) <- paste0(prefix, result$egt[keep])
    out <- c(out, values)
  }
  out
}

make_nonlinear_fixture <- function(H = 20L, periods = 1:6,
                                   observations_per_cell = 15L,
                                   weighted = TRUE) {
  cells <- expand.grid(region = seq_len(H), time = periods,
                       observation = seq_len(observations_per_cell))
  cohort_by_region <- rep(0, H)
  treated_count <- floor(0.6 * H)
  cohort_by_region[seq_len(treated_count)] <- rep(
    c(3, 4, 5), length.out = treated_count
  )
  cells$cohort <- cohort_by_region[cells$region]
  cells$id <- seq_len(nrow(cells))
  cells$x1 <- sin(cells$id / 11) + cos(cells$region / 7) / 5 +
    cells$time / 30
  cells$x2 <- cos(cells$id / 17) - sin(cells$region / 9) / 6
  cells$female <- as.integer((cells$observation + cells$region) %% 2 == 0)
  cells$sector <- factor(c("a", "b", "c")[
    1 + (cells$observation + cells$region) %% 3
  ])
  treated <- cells$cohort > 0 & cells$time >= cells$cohort
  cells$y <- 0.2 * cells$region + 0.15 * cells$time +
    0.4 * cells$x1 - 0.25 * cells$x2 + 0.1 * cells$female +
    0.8 * treated + sin(cells$id / 5) / 4
  cells$w <- if (weighted) 0.2 + (cells$id %% 23) / 13 else 1
  cells
}
