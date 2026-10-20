make_regression_fixture <- function(unbalanced = FALSE, weighted = FALSE,
                                    missing_cells = FALSE) {
  data <- did_fixture(unbalanced = unbalanced, weighted = weighted,
                      missing_cells = missing_cells)
  time_index <- match(data$time, sort(unique(data$time)))
  data$x1 <- sin(data$id / 7) + data$region / 9 + time_index / 13
  data$x2 <- cos(data$id / 11) - data$region / 17 + data$id %% 5L / 8
  data$female <- as.integer((data$id + data$region) %% 2L == 0L)
  data$sector <- factor(c("agriculture", "manufacturing", "services")[
    1L + (data$id + 2L * data$region) %% 3L
  ])
  data$y <- data$y + 0.45 * data$x1 - 0.2 * data$x2 + 0.3 * data$female +
    c(agriculture = -0.1, manufacturing = 0.15, services = 0.4)[data$sector]
  data
}

reference_reg_map <- function(data, never_treated = 0) {
  time <- sort(unique(numeric_labels(data$time)))
  cohort <- numeric_labels(data$cohort)
  never <- numeric_labels(never_treated)
  groups <- sort(unique(cohort[cohort != never]))
  rows <- lapply(groups, function(group) {
    pre <- time[time < group]
    post <- time[time >= group]
    if (!length(pre) || !length(post)) return(NULL)
    data.frame(
      group = group, base = max(pre), post = post, event = post - group,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out$att_id <- seq_len(nrow(out))
  out[, c("att_id", "group", "base", "post", "event")]
}

reference_reg_prepare <- function(data, xformla, weights = NULL,
                                  never_treated = 0) {
  model_frame <- stats::model.frame(
    xformla, data = data, na.action = stats::na.pass, drop.unused.levels = TRUE
  )
  X <- stats::model.matrix(xformla, data = model_frame)
  list(
    data = data,
    X = X,
    y = numeric_labels(data$y),
    time = numeric_labels(data$time),
    cohort = replace(
      numeric_labels(data$cohort),
      numeric_labels(data$cohort) == numeric_labels(never_treated), Inf
    ),
    w = if (is.null(weights)) rep(1, nrow(data)) else numeric_labels(data[[weights]])
  )
}

reference_wls <- function(X, y, w) {
  positive <- is.finite(w) & w > 0
  X <- X[positive, , drop = FALSE]
  y <- y[positive]
  w <- w[positive]
  if (!nrow(X) || sum(w) <= 0 || nrow(X) < ncol(X)) return(NULL)
  Q <- crossprod(X, w * X)
  if (!all(is.finite(Q)) || rcond(Q) < .Machine$double.eps) return(NULL)
  fit <- stats::lm.wfit(X, y, w)
  if (fit$rank != ncol(X) || anyNA(fit$coefficients)) return(NULL)
  unname(fit$coefficients)
}

reference_reg_att <- function(prepared, att_map,
                              control = c("never", "not_yet")) {
  control <- match.arg(control)
  with(prepared, vapply(seq_len(nrow(att_map)), function(k) {
    group <- att_map$group[k]
    base <- att_map$base[k]
    post <- att_map$post[k]
    treated <- cohort == group
    controls <- if (control == "never") is.infinite(cohort) else
      is.infinite(cohort) | cohort > post
    control_two_period <- controls & time %in% c(base, post)
    pooled_X <- X[control_two_period, , drop = FALSE]
    if (nrow(pooled_X) < ncol(X) ||
        rcond(crossprod(pooled_X)) < .Machine$double.eps) return(NA_real_)

    base_control <- controls & time == base
    post_control <- controls & time == post
    beta_base <- reference_wls(X[base_control, , drop = FALSE],
                               y[base_control], w[base_control])
    beta_post <- reference_wls(X[post_control, , drop = FALSE],
                               y[post_control], w[post_control])
    treated_base <- treated & time == base & w > 0
    treated_post <- treated & time == post & w > 0
    if (is.null(beta_base) || is.null(beta_post) ||
        !any(treated_base) || !any(treated_post) ||
        sum(w[treated_base]) <= 0 || sum(w[treated_post]) <= 0) {
      return(NA_real_)
    }
    y_base <- sum(w[treated_base] * y[treated_base]) / sum(w[treated_base])
    y_post <- sum(w[treated_post] * y[treated_post]) / sum(w[treated_post])
    treated_pool <- treated_base | treated_post
    x_pool <- colSums(X[treated_pool, , drop = FALSE] * w[treated_pool]) /
      sum(w[treated_pool])
    as.numeric(y_post - y_base - crossprod(x_pool, beta_post - beta_base))
  }, numeric(1)))
}

reference_reg_fit <- function(data, xformla, control = "never", weights = NULL,
                              att_map = reference_reg_map(data)) {
  prepared <- reference_reg_prepare(data, xformla, weights)
  att <- reference_reg_att(prepared, att_map, control)
  names(att) <- paste0("ATT(", att_map$group, ",", att_map$post, ")")
  sizes <- reference_group_size(data, sort(unique(att_map$group)), weights)
  list(
    att = att,
    aggregate = reference_aggregate(att, att_map, sizes),
    att_map = att_map,
    group_size = sizes
  )
}

did_reg_reference <- function(data, xformla, control, weighted = FALSE,
                              compute_inffunc = FALSE) {
  did::att_gt(
    yname = "y", tname = "time", idname = "id", gname = "cohort",
    xformla = xformla, data = data, panel = FALSE,
    control_group = if (control == "never") "nevertreated" else "notyettreated",
    anticipation = 0, weightsname = if (weighted) "w" else NULL,
    bstrap = FALSE, cband = FALSE, est_method = "reg",
    base_period = "universal", faster_mode = TRUE,
    compute_inffunc = compute_inffunc
  )
}
