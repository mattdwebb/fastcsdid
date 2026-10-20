.if_solve_rows <- function(system, row_scores) {
  if (!nrow(row_scores)) return(row_scores)
  solved <- tryCatch(
    t(solve(system, t(row_scores))),
    error = function(e) matrix(NA_real_, nrow(row_scores), ncol(row_scores))
  )
  unname(solved)
}

.retain_influence_data <- function(moments, data, outcome, time, cohort,
                                   cluster, weights) {
  if (!is.null(moments$regression_data)) {
    regression <- moments$regression_data
    return(list(
      y = regression$y, weight = regression$weight,
      region_index = regression$region_index,
      time_index = regression$time_index,
      assignment = regression$assignment,
      original_rows = regression$original_rows,
      original_n = nrow(data)
    ))
  }
  region_index <- match(data[[cluster]], moments$regions)
  time_index <- match(.as_numeric_design(data[[time]], "time"), moments$times)
  list(
    y = .as_numeric_design(data[[outcome]], "outcome"),
    weight = if (is.null(weights)) rep.int(1, nrow(data)) else
      .as_numeric_design(data[[weights]], "weights"),
    region_index = region_index, time_index = time_index,
    assignment = moments$assignment[region_index],
    original_rows = seq_len(nrow(data)), original_n = nrow(data)
  )
}

.if_ratio <- function(weight, value) {
  mean_weight <- mean(weight)
  if (!is.finite(mean_weight) || mean_weight <= 0) {
    return(list(estimate = NA_real_, influence = rep(NA_real_, length(weight))))
  }
  estimate <- mean(weight * value) / mean_weight
  list(
    estimate = estimate,
    influence = weight * (value - estimate) / mean_weight
  )
}

.if_wls_representation <- function(X, y, weight, coefficients) {
  n <- nrow(X)
  residual <- y - as.numeric(X %*% coefficients)
  scores <- X * (weight * residual)
  information <- crossprod(X, weight * X) / n
  .if_solve_rows(information, scores)
}

.if_propensity_representation <- function(X, D, weight, probability) {
  n <- nrow(X)
  scores <- X * (weight * (D - probability))
  information <- crossprod(X, (weight * probability * (1 - probability)) * X) / n
  .if_solve_rows(information, scores)
}

.local_none_if <- function(y, post, D, weight) {
  treated_pre <- .if_ratio(weight * D * (1 - post), y)
  treated_post <- .if_ratio(weight * D * post, y)
  control_pre <- .if_ratio(weight * (1 - D) * (1 - post), y)
  control_post <- .if_ratio(weight * (1 - D) * post, y)
  treated_post$influence - treated_pre$influence -
    control_post$influence + control_pre$influence
}

.local_reg_if <- function(y, post, D, X, weight, coefficients) {
  n <- length(y)
  weight <- weight / mean(weight)
  control_pre <- as.numeric(X %*% coefficients$control_base)
  control_post <- as.numeric(X %*% coefficients$control_post)
  w_treated_pre <- weight * D * (1 - post)
  w_treated_post <- weight * D * post
  w_treated <- weight * D
  eta_pre <- mean(w_treated_pre * y) / mean(w_treated_pre)
  eta_post <- mean(w_treated_post * y) / mean(w_treated_post)
  eta_control <- mean(w_treated * (control_post - control_pre)) /
    mean(w_treated)
  treated_pre <- w_treated_pre * (y - eta_pre) / mean(w_treated_pre)
  treated_post <- w_treated_post * (y - eta_post) / mean(w_treated_post)

  ols_pre <- .if_wls_representation(
    X, y, weight * (1 - D) * (1 - post), coefficients$control_base
  )
  ols_post <- .if_wls_representation(
    X, y, weight * (1 - D) * post, coefficients$control_post
  )
  M1 <- as.numeric(crossprod(w_treated, X)) / n
  control_direct <- w_treated *
    (control_post - control_pre - eta_control)
  control <- (control_direct + as.numeric(ols_post %*% M1) -
                as.numeric(ols_pre %*% M1)) / mean(w_treated)
  treated_post - treated_pre - control
}

.local_ipw_if <- function(y, post, D, X, weight, probability) {
  n <- length(y)
  weight <- weight / mean(weight)
  probability <- pmin(probability, 1 - 1e-6)
  trim <- D == 1 | probability < 0.995
  odds <- probability / (1 - probability)
  w10 <- trim * weight * D * (1 - post)
  w11 <- trim * weight * D * post
  w00 <- trim * weight * odds * (1 - D) * (1 - post)
  w01 <- trim * weight * odds * (1 - D) * post
  r10 <- .if_ratio(w10, y); r11 <- .if_ratio(w11, y)
  r00 <- .if_ratio(w00, y); r01 <- .if_ratio(w01, y)
  propensity <- .if_propensity_representation(
    X, D, weight, probability
  )
  M2_pre <- as.numeric(crossprod(w00 * (y - r00$estimate), X)) /
    n / mean(w00)
  M2_post <- as.numeric(crossprod(w01 * (y - r01$estimate), X)) /
    n / mean(w01)
  propensity_effect <- as.numeric(propensity %*% (M2_post - M2_pre))
  (r11$influence - r10$influence) -
    (r01$influence - r00$influence + propensity_effect)
}

.local_dr_if <- function(y, post, D, X, weight, probability, coefficients) {
  n <- length(y)
  weight <- weight / mean(weight)
  probability <- pmin(probability, 1 - 1e-6)
  trim <- D == 1 | probability < 0.995
  odds <- probability / (1 - probability)
  m00 <- as.numeric(X %*% coefficients$control_base)
  m01 <- as.numeric(X %*% coefficients$control_post)
  m10 <- as.numeric(X %*% coefficients$treated_base)
  m11 <- as.numeric(X %*% coefficients$treated_post)
  m0 <- post * m01 + (1 - post) * m00

  w10 <- trim * weight * D * (1 - post)
  w11 <- trim * weight * D * post
  w00 <- trim * weight * odds * (1 - D) * (1 - post)
  w01 <- trim * weight * odds * (1 - D) * post
  wd <- trim * weight * D
  wd1 <- trim * weight * D * post
  wd0 <- trim * weight * D * (1 - post)
  r10 <- .if_ratio(w10, y - m0); r11 <- .if_ratio(w11, y - m0)
  r00 <- .if_ratio(w00, y - m0); r01 <- .if_ratio(w01, y - m0)
  ed_post <- .if_ratio(wd, m11 - m01)
  ed1_post <- .if_ratio(wd1, m11 - m01)
  ed_pre <- .if_ratio(wd, m10 - m00)
  ed0_pre <- .if_ratio(wd0, m10 - m00)

  ols00 <- .if_wls_representation(
    X, y, weight * (1 - D) * (1 - post), coefficients$control_base
  )
  ols01 <- .if_wls_representation(
    X, y, weight * (1 - D) * post, coefficients$control_post
  )
  ols10 <- .if_wls_representation(
    X, y, weight * D * (1 - post), coefficients$treated_base
  )
  ols11 <- .if_wls_representation(
    X, y, weight * D * post, coefficients$treated_post
  )
  propensity <- .if_propensity_representation(X, D, weight, probability)

  M1_post <- -as.numeric(crossprod(w11, X)) / n / mean(w11)
  M1_pre <- -as.numeric(crossprod(w10, X)) / n / mean(w10)
  treated_or <- as.numeric(ols01 %*% M1_post) +
    as.numeric(ols00 %*% M1_pre)

  M2_pre <- as.numeric(crossprod(
    w00 * (y - m0 - r00$estimate), X
  )) / n / mean(w00)
  M2_post <- as.numeric(crossprod(
    w01 * (y - m0 - r01$estimate), X
  )) / n / mean(w01)
  control_ps <- as.numeric(propensity %*% (M2_post - M2_pre))
  M3_post <- -as.numeric(crossprod(w01, X)) / n / mean(w01)
  M3_pre <- -as.numeric(crossprod(w00, X)) / n / mean(w00)
  control_or <- as.numeric(ols01 %*% M3_post) +
    as.numeric(ols00 %*% M3_pre)

  efficient_direct <- (ed_post$influence - ed1_post$influence) -
    (ed_pre$influence - ed0_pre$influence)
  moment_post <- as.numeric(crossprod(
    wd / mean(wd) - wd1 / mean(wd1), X
  )) / n
  moment_pre <- as.numeric(crossprod(
    wd / mean(wd) - wd0 / mean(wd0), X
  )) / n
  efficient_or <- as.numeric((ols11 - ols01) %*% moment_post) -
    as.numeric((ols10 - ols00) %*% moment_pre)

  treated <- r11$influence - r10$influence + treated_or
  control <- r01$influence - r00$influence + control_ps + control_or
  treated - control + efficient_direct + efficient_or
}

.influence_observation_data <- function(moments) {
  data <- moments$influence_data
  if (is.null(data)) {
    stop(
      "Observation-level influence functions are unavailable: this fit was created from moments without retained microdata. Cluster-summed influence remains available.",
      call. = FALSE
    )
  }
  data
}

.att_if_column <- function(factorized, k) {
  moments <- factorized$moments
  map <- factorized$att_map
  N <- moments$N
  if (!factorized$support[k]) return(rep(NA_real_, N))
  method <- factorized$est_method
  if (method == "none") {
    data <- .influence_observation_data(moments)
    treated <- data$assignment == map$group[k]
    controls <- if (factorized$control_group == "never") {
      is.infinite(data$assignment)
    } else {
      is.infinite(data$assignment) | data$assignment > map$post[k]
    }
    rows <- which((treated | controls) &
                    data$time_index %in% match(c(map$base[k], map$post[k]),
                                                moments$times))
    D <- as.numeric(treated[rows])
    post <- as.numeric(data$time_index[rows] == match(map$post[k], moments$times))
    local <- .local_none_if(data$y[rows], post, D, data$weight[rows])
  } else if (method == "reg") {
    data <- moments$regression_data
    definition <- .propensity_sample_rows(
      data, moments, map[k, ], factorized$control_group
    )
    rows <- definition$rows
    D <- as.numeric(data$assignment[rows] == map$group[k])
    post <- as.numeric(data$time_index[rows] == match(map$post[k], moments$times))
    nuisance <- factorized$evaluation$nuisance
    coefficients <- list(
      control_base = nuisance$blocks[[nuisance$cell_blocks[k, "base_block"]]]$coefficients,
      control_post = nuisance$blocks[[nuisance$cell_blocks[k, "post_block"]]]$coefficients
    )
    local <- .local_reg_if(
      data$y[rows], post, D, data$X[rows, , drop = FALSE],
      data$weight[rows], coefficients
    )
  } else {
    data <- moments$regression_data
    block <- factorized$evaluation$propensity$blocks[[k]]
    rows <- block$rows
    X <- data$X[rows, , drop = FALSE]
    probability <- pmin(stats::plogis(as.numeric(X %*% block$coefficients)),
                        1 - 1e-6)
    if (method == "ipw") {
      local <- .local_ipw_if(
        data$y[rows], block$post, block$D, X, data$weight[rows], probability
      )
    } else {
      coefficients <- .nonlinear_outcome_coefficients(
        k, factorized$evaluation$control_nuisance,
        factorized$evaluation$treated_nuisance
      )
      local <- .local_dr_if(
        data$y[rows], block$post, block$D, X, data$weight[rows],
        probability, coefficients
      )
    }
  }
  out <- numeric(N)
  out[rows] <- N / length(rows) * local
  out
}

.cluster_none_influence <- function(moments, att_map, control_group, support) {
  H <- length(moments$regions)
  K <- nrow(att_map)
  N <- moments$N
  out <- matrix(NA_real_, H, K,
                dimnames = list(as.character(moments$regions), .att_names(att_map)))
  for (k in seq_len(K)) {
    if (!support[k]) next
    b <- match(att_map$base[k], moments$times)
    t <- match(att_map$post[k], moments$times)
    treated <- moments$assignment == att_map$group[k]
    control <- if (control_group == "never") is.infinite(moments$assignment) else
      is.infinite(moments$assignment) | moments$assignment > att_map$post[k]
    cell <- function(regions, period, sign) {
      total_weight <- sum(moments$weight_sum[regions, period])
      mean_y <- sum(moments$outcome_sum[regions, period]) / total_weight
      sign * N * as.numeric(regions) *
        (moments$outcome_sum[, period] - mean_y * moments$weight_sum[, period]) /
        total_weight
    }
    out[, k] <- cell(treated, t, 1) + cell(treated, b, -1) +
      cell(control, t, -1) + cell(control, b, 1)
  }
  out
}

.rowsum_if <- function(values, region_index, H) {
  table <- rowsum(values, region_index, reorder = FALSE)
  out <- numeric(H)
  out[as.integer(rownames(table))] <- table[, 1L]
  out
}

.aggregate_selector <- function(att_map, metadata_row) {
  switch(
    metadata_row$type,
    simple = rep(TRUE, nrow(att_map)),
    group = att_map$group == metadata_row$index,
    calendar = att_map$post == metadata_row$index,
    dynamic = att_map$event == metadata_row$index
  )
}

.aggregate_weight_cluster_if <- function(factorized, a) {
  aggregate <- factorized$aggregate
  metadata <- aggregate$metadata[a, ]
  H <- length(factorized$moments$regions)
  if (!is.finite(metadata$estimate) || metadata$type == "group") return(numeric(H))
  selector <- .aggregate_selector(factorized$att_map, metadata)
  groups <- factorized$att_map$group[selector]
  effects <- factorized$att[selector] - metadata$estimate
  masses <- factorized$evaluation$group_size[as.character(groups)]
  denominator <- sum(masses)
  contribution_by_group <- tapply(effects, groups, sum)
  assignment <- factorized$moments$assignment
  contribution <- numeric(H)
  finite <- is.finite(assignment) & as.character(assignment) %in%
    names(contribution_by_group)
  contribution[finite] <- factorized$moments$N *
    factorized$moments$region_total_weight[finite] *
    contribution_by_group[as.character(assignment[finite])] / denominator
  contribution
}

.aggregate_weight_observation_if <- function(factorized, a) {
  data <- .influence_observation_data(factorized$moments)
  aggregate <- factorized$aggregate
  metadata <- aggregate$metadata[a, ]
  if (!is.finite(metadata$estimate) || metadata$type == "group") {
    return(numeric(factorized$moments$N))
  }
  selector <- .aggregate_selector(factorized$att_map, metadata)
  groups <- factorized$att_map$group[selector]
  effects <- factorized$att[selector] - metadata$estimate
  masses <- factorized$evaluation$group_size[as.character(groups)]
  denominator <- sum(masses)
  contribution_by_group <- tapply(effects, groups, sum)
  out <- numeric(factorized$moments$N)
  finite <- is.finite(data$assignment) & as.character(data$assignment) %in%
    names(contribution_by_group)
  out[finite] <- factorized$moments$N * data$weight[finite] *
    contribution_by_group[as.character(data$assignment[finite])] / denominator
  out
}

.analytical_covariance <- function(cluster, N) {
  labels <- colnames(cluster)
  blank <- matrix(NA_real_, ncol(cluster), ncol(cluster),
                  dimnames = list(labels, labels))
  keep <- apply(is.finite(cluster), 2L, all)
  if (any(keep)) blank[keep, keep] <- crossprod(cluster[, keep, drop = FALSE]) / N^2
  blank
}

.build_influence_system <- function(moments, att_map, evaluation, aggregate,
                                    est_method, control_group) {
  cluster_start <- proc.time()[["elapsed"]]
  factorized <- list(
    moments = moments, att_map = att_map, evaluation = evaluation,
    aggregate = aggregate, att = evaluation$att,
    support = evaluation$support$supported, est_method = est_method,
    control_group = control_group
  )
  H <- length(moments$regions)
  K <- nrow(att_map)
  cluster_att <- if (est_method == "none") {
    .cluster_none_influence(
      moments, att_map, control_group, evaluation$support$supported
    )
  } else {
    out <- matrix(NA_real_, H, K,
                  dimnames = list(as.character(moments$regions), .att_names(att_map)))
    region_index <- moments$regression_data$region_index
    for (k in seq_len(K)) {
      column <- .att_if_column(factorized, k)
      if (all(is.finite(column))) out[, k] <- .rowsum_if(column, region_index, H)
    }
    out
  }
  A <- length(aggregate$estimates)
  cluster_aggregate <- matrix(
    NA_real_, H, A,
    dimnames = list(as.character(moments$regions), names(aggregate$estimates))
  )
  for (a in seq_len(A)) {
    weights <- aggregate$weights[a, ]
    if (all(is.finite(weights))) {
      needed <- which(weights != 0)
      # Unsupported ATT columns are NA. Subset before multiplication because
      # R propagates NA * 0, which would otherwise contaminate a supported
      # calendar/event aggregate that does not use an unsupported cell.
      if (length(needed) && all(is.finite(cluster_att[, needed, drop = FALSE]))) {
        cluster_aggregate[, a] <- as.numeric(
          cluster_att[, needed, drop = FALSE] %*% weights[needed]
        ) + .aggregate_weight_cluster_if(factorized, a)
      }
    }
  }
  cluster_elapsed <- proc.time()[["elapsed"]] - cluster_start
  covariance_start <- proc.time()[["elapsed"]]
  joint <- cbind(cluster_att, cluster_aggregate)
  raw <- .analytical_covariance(joint, moments$N)
  df_adjusted <- length(moments$regions) / (length(moments$regions) - 1) * raw
  covariance_elapsed <- proc.time()[["elapsed"]] - covariance_start
  list(
    cluster_att = cluster_att,
    cluster_aggregate = cluster_aggregate,
    metadata = list(
      scale = "global empirical influence: estimator perturbation is N^-1 times the influence sum",
      normalization = "DRDID local influence multiplied by N / n_att",
      observations = moments$N, clusters = length(moments$regions),
      att_order = colnames(cluster_att),
      aggregate_order = colnames(cluster_aggregate),
      observation_available = !is.null(moments$influence_data),
      aggregation_weight_influence = TRUE,
      timing = c(cluster_influence = cluster_elapsed,
                 analytical_covariance = covariance_elapsed)
    ),
    factorized = factorized,
    covariance = list(
      analytical_raw = raw,
      analytical = raw,
      analytical_df = df_adjusted,
      metadata = list(
        reference = "did 2.5.1 getSE cluster convention",
        raw_scale = "crossprod(cluster influence) / N^2",
        analytical_adjustment = "none",
        df_adjustment = "H / (H - 1)",
        parameter_count_adjustment = "none"
      )
    )
  )
}

.resolve_influence_estimands <- function(fit, estimand, type) {
  labels <- if (type == "att") colnames(fit$influence$cluster_att) else
    colnames(fit$influence$cluster_aggregate)
  if (is.null(estimand)) return(seq_along(labels))
  position <- match(estimand, labels)
  if (anyNA(position)) stop("Unknown ", type, " estimand: ",
                            paste(estimand[is.na(position)], collapse = ", "),
                            call. = FALSE)
  position
}

#' Extract influence functions
#'
#' Observation output uses the global empirical convention
#' `estimate - target = mean(influence) + o_p(N^-1/2)`. Cluster output is the
#' stored within-cluster sum of those observation values; it is not divided by
#' `N`. Thus the raw analytical covariance is the cross-product of cluster sums
#' divided by `N^2`.
#'
#' @param fit A fitted `fastcsdid` object.
#' @param estimand Optional ATT or aggregate labels. `NULL` selects all labels.
#' @param type Either `"att"` or `"aggregate"`.
#' @param cluster_sum Return stored cluster sums instead of observation values.
#' @param batch_size Maximum ATT columns generated in one observation batch.
#' @return A matrix. Observation output carries an `observation_metadata`
#'   attribute recording original row, cluster, period, and cohort. A
#'   moments-only fit without retained microdata can return cluster sums but
#'   cannot return observation values.
#' @export
influence_function <- function(fit, estimand = NULL,
                               type = c("att", "aggregate"),
                               cluster_sum = FALSE, batch_size = 32L) {
  if (!inherits(fit, "fastcsdid") || is.null(fit$influence)) {
    stop("fit does not contain influence-function ingredients.", call. = FALSE)
  }
  type <- match.arg(type)
  position <- .resolve_influence_estimands(fit, estimand, type)
  stored <- if (type == "att") fit$influence$cluster_att else
    fit$influence$cluster_aggregate
  if (cluster_sum) return(stored[, position, drop = FALSE])
  data <- .influence_observation_data(fit$influence$factorized$moments)
  batch_size <- as.integer(batch_size)
  if (length(batch_size) != 1L || is.na(batch_size) || batch_size < 1L) {
    stop("batch_size must be a positive integer.", call. = FALSE)
  }
  factorized <- fit$influence$factorized
  if (type == "att") {
    out <- matrix(NA_real_, fit$moments$N, length(position),
                  dimnames = list(NULL, colnames(stored)[position]))
    starts <- seq.int(1L, length(position), by = batch_size)
    for (start in starts) {
      local <- start:min(length(position), start + batch_size - 1L)
      out[, local] <- vapply(position[local], function(k) {
        .att_if_column(factorized, k)
      }, numeric(fit$moments$N))
    }
  } else {
    out <- matrix(NA_real_, fit$moments$N, length(position),
                  dimnames = list(NULL, colnames(stored)[position]))
    for (j in seq_along(position)) {
      a <- position[j]
      weights <- fit$aggregation$weights[a, ]
      if (!all(is.finite(weights))) next
      needed <- which(weights != 0)
      att_if <- vapply(needed, function(k) .att_if_column(factorized, k),
                       numeric(fit$moments$N))
      out[, j] <- as.numeric(att_if %*% weights[needed]) +
        .aggregate_weight_observation_if(factorized, a)
    }
  }
  attr(out, "observation_metadata") <- data.frame(
    original_row = data$original_rows,
    cluster = fit$moments$regions[data$region_index],
    period = fit$moments$times[data$time_index],
    cohort = fit$moments$assignment_raw[data$region_index],
    stringsAsFactors = FALSE
  )
  out
}

#' Extract recentered influence functions
#'
#' A recentered influence function is the fitted estimate plus its observation
#' influence value. Its sample column mean is therefore the fitted estimate up
#' to numerical tolerance. Aggregate RIFs include estimated cohort-mass weight
#' effects where applicable.
#'
#' @inheritParams influence_function
#' @param estimand One or more ATT or aggregate labels.
#' @param type `"auto"` infers the estimand family, or specify `"att"` or
#'   `"aggregate"`.
#' @param include_estimate Add the fitted estimand to each influence value.
#' @return A matrix in original retained-observation order.
#' @export
rif <- function(fit, estimand, type = c("auto", "att", "aggregate"),
                include_estimate = TRUE, batch_size = 32L) {
  type <- match.arg(type)
  if (type == "auto") {
    in_att <- estimand %in% names(fit$att_gt)
    in_aggregate <- estimand %in% names(fit$aggregation$estimates)
    if (all(in_att) && !any(in_aggregate)) type <- "att" else
      if (all(in_aggregate) && !any(in_att)) type <- "aggregate" else
        stop("estimand labels do not identify one unique estimand type.", call. = FALSE)
  }
  out <- influence_function(
    fit, estimand = estimand, type = type, cluster_sum = FALSE,
    batch_size = batch_size
  )
  if (isTRUE(include_estimate)) {
    metadata <- attr(out, "observation_metadata")
    estimates <- if (type == "att") fit$att_gt[colnames(out)] else
      fit$aggregation$estimates[colnames(out)]
    out <- sweep(out, 2L, estimates, `+`)
    attr(out, "observation_metadata") <- metadata
  }
  out
}

#' Compare analytical and jackknife inference
#'
#' Analytical intervals use a normal reference distribution and the
#' `did`-compatible unadjusted cluster influence covariance. Jackknife intervals
#' use a `t(H - 1)` reference distribution and explicitly labelled pairwise
#' partial-support CV3 results. Missing methods remain `NA`.
#'
#' @param fit A fitted `fastcsdid` object.
#' @param level Confidence level.
#' @return A data frame with analytical and available jackknife results.
#' @export
inference_comparison <- function(fit, level = 0.95) {
  if (!inherits(fit, "fastcsdid") || is.null(fit$vcov)) {
    stop("fit does not contain analytical covariance results.", call. = FALSE)
  }
  estimate <- c(fit$att_gt, fit$aggregation$estimates)
  labels <- names(estimate)
  analytical_se <- sqrt(diag(fit$vcov$analytical))[labels]
  z <- stats::qnorm(1 - (1 - level) / 2)
  out <- data.frame(
    estimand = labels, estimate = unname(estimate),
    analytical_se = unname(analytical_se),
    analytical_conf_low = unname(estimate - z * analytical_se),
    analytical_conf_high = unname(estimate + z * analytical_se),
    analytical_statistic = unname(estimate / analytical_se),
    analytical_p_value = 2 * stats::pnorm(-abs(estimate / analytical_se)),
    linearized_cv3_se = NA_real_, linearized_cv3j_se = NA_real_,
    refit_cv3_se = NA_real_, refit_cv3j_se = NA_real_,
    linearized_conf_low = NA_real_, linearized_conf_high = NA_real_,
    linearized_statistic = NA_real_, linearized_p_value = NA_real_,
    refit_conf_low = NA_real_, refit_conf_high = NA_real_,
    refit_statistic = NA_real_, refit_p_value = NA_real_,
    linearized_refit_ratio = NA_real_, analytical_refit_ratio = NA_real_,
    n_valid_deletions = NA_integer_,
    support_status = c(fit$att_support$status,
                       ifelse(is.finite(fit$aggregation$estimates), "valid",
                              "undefined_aggregation")),
    analytical_distribution = "normal", analytical_df = Inf,
    jackknife_distribution = paste0("t(", length(fit$moments$regions) - 1L, ")"),
    jackknife_covariance_label = "pairwise partial-support CV3 (not complete-support CV3 when deletions are missing)",
    stringsAsFactors = FALSE
  )
  if (!is.null(fit$jackknife)) {
    add_method <- function(result, prefix) {
      se3 <- c(result$att_covariance$pairwise_se_cv3,
               result$aggregate_covariance$pairwise_se_cv3)
      se3j <- c(result$att_covariance$pairwise_se_cv3j,
                result$aggregate_covariance$pairwise_se_cv3j)
      out[[paste0(prefix, "_cv3_se")]] <<- unname(se3[labels])
      out[[paste0(prefix, "_cv3j_se")]] <<- unname(se3j[labels])
    }
    if (!is.null(fit$jackknife$methods$linearized)) {
      add_method(fit$jackknife$methods$linearized, "linearized")
    } else if (identical(fit$jackknife$method, "exact")) {
      # Linear estimators have an exact sparse deletion, so there is no
      # approximation/refit distinction. Put that result in the primary
      # comparison columns while retaining method = "exact" on the fit.
      add_method(fit$jackknife, "linearized")
    }
    if (!is.null(fit$jackknife$methods$refit)) {
      add_method(fit$jackknife$methods$refit, "refit")
    }
    valid <- c(fit$jackknife$att_covariance$n_valid_by_stat,
               fit$jackknife$aggregate_covariance$n_valid_by_stat)
    out$n_valid_deletions <- unname(valid[labels])
    out$linearized_refit_ratio <- out$linearized_cv3_se / out$refit_cv3_se
    out$analytical_refit_ratio <- out$analytical_se / out$refit_cv3_se
    degrees <- length(fit$moments$regions) - 1L
    critical <- stats::qt(1 - (1 - level) / 2, degrees)
    for (prefix in c("linearized", "refit")) {
      standard_error <- out[[paste0(prefix, "_cv3_se")]]
      statistic <- out$estimate / standard_error
      out[[paste0(prefix, "_conf_low")]] <- out$estimate -
        critical * standard_error
      out[[paste0(prefix, "_conf_high")]] <- out$estimate +
        critical * standard_error
      out[[paste0(prefix, "_statistic")]] <- statistic
      out[[paste0(prefix, "_p_value")]] <-
        2 * stats::pt(-abs(statistic), degrees)
    }
  }
  out
}
