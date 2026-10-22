#' Fit the estimator from pre-collapsed moments
#'
#' Uses pre-computed region-period moments to avoid revisiting individual
#' observations.
#'
#' @param moments An object returned by [collapse_rcs()].
#' @param control_group `"never"` or `"not_yet"`.
#' @param aggregation Requested aggregation types.
#' @param inference `"none"` or `"jackknife"`.
#' @param est_method `"none"`, `"reg"`, `"ipw"`, or `"dr"`.
#' @param ps_fit Propensity starting strategy: a master-model one-step warm
#'   start or a conventional zero start. Exact child logits converge fully.
#' @param jackknife_method Nonlinear propensity deletion method. `"linearized"`
#'   uses a deleted-information Newton step and first-order propensity effect;
#'   `"refit"` converges every affected deleted logit; `"both"` retains both
#'   paths for direct comparison. Linear estimators use exact
#'   deletion.
#' @param jackknife_support `"changing"` (default) recomputes cohort masses and
#'   aggregation after each deletion, allowing zero-mass cohorts to disappear.
#'   `"fixed"` additionally requires every positive-weight full-sample ATT cell
#'   needed by an aggregate to survive. Both policies record support changes;
#'   undefined aggregates still fail and are never silently dropped.
#' @return An object of class `fastcsdid`.
#' @export
fit_csdid_moments <- function(moments,
                              control_group = c("never", "not_yet"),
                              aggregation = c("simple", "group", "calendar", "dynamic"),
                              inference = c("jackknife", "none"),
                              est_method = c("none", "reg", "ipw", "dr"),
                              ps_fit = c("warm_start", "standard"),
                              jackknife_method = c("linearized", "refit", "both"),
                              jackknife_support = c("changing", "fixed")) {
  if (!inherits(moments, "fastrcsdid_moments")) {
    stop("moments must be returned by collapse_rcs().", call. = FALSE)
  }
  control_group <- match.arg(control_group)
  inference <- match.arg(inference)
  est_method <- match.arg(est_method)
  ps_fit <- match.arg(ps_fit)
  jackknife_method <- match.arg(jackknife_method)
  jackknife_support <- match.arg(jackknife_support)
  aggregation <- .match_aggregation_types(aggregation)
  if (est_method %in% c("reg", "ipw", "dr") &&
      is.null(moments$regression_data)) {
    stop("est_method = '", est_method,
         "' requires moments created with xformla.", call. = FALSE)
  }
  att_map <- compile_att_map(moments)
  regression <- NULL
  propensity_public <- NULL
  if (est_method == "none") {
    evaluation <- evaluate_assignment(
      moments, moments$assignment, att_map, control_group
    )
  } else if (est_method == "reg") {
    nuisance <- .fit_nuisance_blocks(moments, att_map, control_group)
    evaluation <- .evaluate_att_from_blocks(moments, att_map, nuisance)
    regression <- list(
      formula = moments$regression_data$formula,
      terms = moments$regression_data$terms,
      column_names = moments$regression_data$column_names,
      contrasts = moments$regression_data$contrasts,
      blocks = nuisance$metadata,
      cell_blocks = nuisance$cell_blocks,
      coefficients = do.call(rbind, lapply(nuisance$blocks, `[[`, "coefficients")),
      timing = c(moments$regression_data$timing, nuisance$timing,
                 evaluation$timing)
    )
    rownames(regression$coefficients) <- paste0("block:", nuisance$metadata$block_id)
    colnames(regression$coefficients) <- regression$column_names
  } else {
    propensity <- .fit_propensity_blocks(
      moments, att_map, control_group, ps_fit
    )
    control_nuisance <- treated_nuisance <- NULL
    if (est_method == "dr") {
      control_nuisance <- .fit_nuisance_blocks(
        moments, att_map, control_group
      )
      treated_nuisance <- .fit_treated_regression_blocks(moments, att_map)
      regression <- list(
        formula = moments$regression_data$formula,
        terms = moments$regression_data$terms,
        column_names = moments$regression_data$column_names,
        contrasts = moments$regression_data$contrasts,
        blocks = control_nuisance$metadata,
        control_blocks = control_nuisance$metadata,
        control_cell_blocks = control_nuisance$cell_blocks,
        treated_blocks = treated_nuisance$metadata,
        treated_cell_blocks = treated_nuisance$cell_blocks,
        control_coefficients = do.call(
          rbind, lapply(control_nuisance$blocks, `[[`, "coefficients")
        ),
        treated_coefficients = do.call(
          rbind, lapply(treated_nuisance$blocks, `[[`, "coefficients")
        ),
        timing = c(
          moments$regression_data$timing,
          stats::setNames(control_nuisance$timing,
                          paste0("control_", names(control_nuisance$timing))),
          stats::setNames(treated_nuisance$timing,
                          paste0("treated_", names(treated_nuisance$timing)))
        )
      )
      colnames(regression$control_coefficients) <-
        colnames(regression$treated_coefficients) <- regression$column_names
    }
    evaluation <- .evaluate_nonlinear_blocks(
      moments, att_map, propensity, est_method,
      control_nuisance, treated_nuisance
    )
    coefficients <- do.call(
      rbind, lapply(propensity$blocks, `[[`, "coefficients")
    )
    rownames(coefficients) <- paste0("block:", seq_len(nrow(coefficients)))
    colnames(coefficients) <- moments$regression_data$column_names
    propensity_public <- list(
      formula = moments$regression_data$formula,
      column_names = moments$regression_data$column_names,
      ps_fit = ps_fit,
      blocks = propensity$metadata,
      masters = propensity$master_metadata,
      coefficients = coefficients,
      timing = c(
        master_fits = sum(propensity$master_metadata$elapsed_seconds),
        att_specific_fits = sum(propensity$metadata$elapsed_seconds)
      )
    )
  }
  aggregate <- aggregate_att(evaluation$att, att_map, evaluation$group_size, aggregation)
  influence <- .build_influence_system(
    moments, att_map, evaluation, aggregate, est_method, control_group
  )
  analytical_vcov <- influence$covariance
  analytical_vcov$analytical_metadata <- analytical_vcov$metadata
  influence$covariance <- NULL

  jackknife <- NULL
  if (inference == "jackknife") {
    build_one <- function(method) {
      delete_att <- if (est_method == "none") {
        .delete_att_matrix(moments, att_map, evaluation$att, control_group)
      } else if (est_method == "reg") {
        .delete_regression_att_matrix(moments, att_map, evaluation)
      } else {
        .delete_nonlinear_att_matrix(
          moments, att_map, evaluation, est_method, method
        )
      }
      delete_aggregate <- .delete_aggregate_matrix(
        moments, att_map, delete_att$estimates, aggregation, jackknife_support
      )
      result <- list(
        att = delete_att$estimates,
        att_changes = delete_att$changes,
        att_status = delete_att$status,
        att_status_label = delete_att$status_label,
        att_reason = delete_att$reason,
        aggregate = delete_aggregate$estimates,
        aggregate_status = delete_aggregate$status,
        aggregate_status_label = delete_aggregate$status_label,
        aggregate_reason = delete_aggregate$reason,
        support = jackknife_support,
        support_changed = delete_aggregate$support_changed,
        status_legend = list(
          att = delete_att$status_legend,
          aggregate = delete_aggregate$status_legend
        ),
        att_covariance = jackknife_covariance(
          delete_att$estimates, evaluation$att, strict = TRUE
        ),
        aggregate_covariance = jackknife_covariance(
          delete_aggregate$estimates, aggregate$estimates, strict = TRUE
        ),
        method = if (est_method %in% c("ipw", "dr")) method else "exact"
      )
      if (!is.null(delete_att$influence)) result$att_influence <- delete_att$influence
      if (!is.null(delete_att$nuisance_diagnostics)) {
        result$nuisance_diagnostics <- delete_att$nuisance_diagnostics
      }
      result
    }
    if (est_method %in% c("ipw", "dr") && jackknife_method == "both") {
      linearized_jackknife <- build_one("linearized")
      refit_jackknife <- build_one("refit")
      jackknife <- linearized_jackknife
      jackknife$methods <- list(
        linearized = linearized_jackknife,
        refit = refit_jackknife
      )
    } else {
      method <- if (est_method %in% c("ipw", "dr")) jackknife_method else "exact"
      jackknife <- build_one(method)
      jackknife$methods <- stats::setNames(list(jackknife), method)
    }
  }

  result <- list(
    call = match.call(),
    moments = moments,
    att_map = att_map,
    att_gt = evaluation$att,
    att_support = evaluation$support,
    aggregation = aggregate,
    jackknife = jackknife,
    regression = regression,
    propensity = propensity_public,
    influence = influence,
    vcov = analytical_vcov,
    diagnostics = list(
      omitted_groups = attr(att_map, "omitted_groups"),
      unsupported_att = data.frame(
        att_map[!evaluation$support$supported, , drop = FALSE],
        evaluation$support[!evaluation$support$supported, , drop = FALSE],
        row.names = NULL
      ),
      jackknife_att_failures = if (is.null(jackknife)) NULL else
        .matrix_failures(jackknife$att_status, jackknife$att_reason,
                          status_label = jackknife$att_status_label),
      jackknife_aggregate_failures = if (is.null(jackknife)) NULL else
        .matrix_failures(jackknife$aggregate_status, jackknife$aggregate_reason,
                          status_label = jackknife$aggregate_status_label),
      support_status_vocabulary = .support_status_vocabulary()
    ),
    settings = list(
      control_group = control_group,
      aggregation = aggregation,
      inference = inference,
      est_method = est_method,
      ps_fit = ps_fit,
      jackknife_method = if (est_method %in% c("ipw", "dr"))
        jackknife_method else "exact",
      weighting = moments$weighting,
      jackknife_support = jackknife_support
    )
  )
  class(result) <- "fastcsdid"
  if (!is.null(jackknife)) {
    result$diagnostics$jackknife_influence <- jackknife_influence_summary(result)
  }
  result
}

#' Fast CSDID for repeated cross sections
#'
#' @param data A data frame containing repeated-cross-section microdata.
#' @param outcome,time,cohort,cluster Column names supplied as strings.
#' @param weights Optional sampling-weight column.
#' @param never_treated Value identifying never-treated regions.
#' @param control_group `"never"` or `"not_yet"`.
#' @param aggregation Requested aggregation types.
#' @param inference `"none"` or `"jackknife"`.
#' @param xformla Optional one-sided covariate formula used by covariate-based
#'   estimators.
#' @param est_method `"none"` estimates effects without covariate adjustment;
#'   `"reg"`, `"ipw"`, and `"dr"` use the corresponding repeated-CS estimator.
#' @param ps_fit,jackknife_method,jackknife_support See [fit_csdid_moments()].
#' @return An object of class `fastcsdid`.
#' @examples
#' d <- expand.grid(region = 1:6, time = 1:4, observation = 1:2)
#' cohort_by_region <- c(2, 2, 3, 3, 0, 0)
#' d$cohort <- cohort_by_region[d$region]
#' d$y <- d$region / 10 + d$time / 5 +
#'   (d$cohort > 0 & d$time >= d$cohort)
#' fit <- fast_csdid(
#'   d, outcome = "y", time = "time", cohort = "cohort", cluster = "region",
#'   inference = "none"
#' )
#' fit$aggregation$estimates
#' d$x <- sin(d$observation + d$region / 3)
#' adjusted <- fast_csdid(
#'   d, outcome = "y", time = "time", cohort = "cohort", cluster = "region",
#'   xformla = ~x, est_method = "reg", inference = "jackknife"
#' )
#' adjusted$aggregation$estimates
#' \dontrun{
#' fit_ipw <- fast_csdid(
#'   earnings_data, "earnings", "year", "first_treat", "state",
#'   weights = "survey_weight",
#'   xformla = ~age + I(age^2) + education + female,
#'   est_method = "ipw", control_group = "never",
#'   inference = "jackknife", jackknife_method = "linearized"
#' )
#' fit_dr <- fast_csdid(
#'   earnings_data, "earnings", "year", "first_treat", "state",
#'   weights = "survey_weight",
#'   xformla = ~age + I(age^2) + education + female,
#'   est_method = "dr", control_group = "never",
#'   inference = "jackknife", jackknife_method = "linearized"
#' )
#' # Use exact refits for nonlinear jackknife inference:
#' fit_dr_exact <- update(fit_dr, jackknife_method = "refit")
#' }
#' @export
fast_csdid <- function(data, outcome, time, cohort, cluster,
                       weights = NULL, never_treated = 0,
                       control_group = c("never", "not_yet"),
                       aggregation = c("simple", "group", "calendar", "dynamic"),
                       inference = c("jackknife", "none"),
                       xformla = NULL,
                       est_method = c("none", "reg", "ipw", "dr"),
                       ps_fit = c("warm_start", "standard"),
                       jackknife_method = c("linearized", "refit", "both"),
                       jackknife_support = c("changing", "fixed")) {
  est_method <- match.arg(est_method)
  ps_fit <- match.arg(ps_fit)
  jackknife_method <- match.arg(jackknife_method)
  jackknife_support <- match.arg(jackknife_support)
  if (est_method != "none" && is.null(xformla)) {
    stop("est_method = '", est_method,
         "' requires xformla; use xformla = ~1 for an intercept-only model.",
         call. = FALSE)
  }
  if (est_method == "none" && !is.null(xformla)) {
    stop("xformla is only used with est_method = 'reg', 'ipw', or 'dr'.",
         call. = FALSE)
  }
  moments <- collapse_rcs(
    data, outcome, time, cohort, cluster, weights, never_treated,
    xformla = xformla
  )
  moments$influence_data <- .retain_influence_data(
    moments, data, outcome, time, cohort, cluster, weights
  )

  fit <- fit_csdid_moments(
    moments = moments,
    control_group = control_group,
    aggregation = aggregation,
    inference = inference,
    est_method = est_method,
    ps_fit = ps_fit,
    jackknife_method = jackknife_method,
    jackknife_support = jackknife_support
  )
  fit$call <- match.call()
  fit
}

#' @export
print.fastcsdid <- function(x, ...) {
  cat("Fast repeated-cross-section CSDID\n")
  cat("Control group:", x$settings$control_group, "\n")
  cat("Estimator:", x$settings$est_method, "\n")
  cat("Regions:", length(x$moments$regions),
      " Periods:", length(x$moments$times),
      " ATT(g,t) cells:", length(x$att_gt), "\n\n")
  if (!is.null(x$regression)) {
    cat("Regression columns:", length(x$regression$column_names),
        " Nuisance blocks:", nrow(x$regression$blocks), "\n\n")
  }
  if (!is.null(x$propensity)) {
    cat("Propensity columns:", length(x$propensity$column_names),
        " Master blocks:", nrow(x$propensity$masters),
        " ATT-specific blocks:", nrow(x$propensity$blocks), "\n\n")
  }
  print(data.frame(
    group = x$att_map$group,
    base = x$att_map$base,
    post = x$att_map$post,
    event = x$att_map$event,
    estimate = unname(x$att_gt),
    supported = x$att_support$supported
  ), row.names = FALSE)
  cat("\nAggregates\n")
  print(x$aggregation$metadata, row.names = FALSE)
  if (!is.null(x$jackknife)) {
    cat("\nJackknife aggregate CV3J standard errors (mean-centered)\n")
    cat("Support policy:", x$settings$jackknife_support, "\n")
    print(x$jackknife$aggregate_covariance$se_cv3j)
  }
  if (!is.null(x$vcov)) {
    aggregate_labels <- names(x$aggregation$estimates)
    cat("\nAnalytical cluster-IF aggregate standard errors\n")
    print(sqrt(diag(x$vcov$analytical))[aggregate_labels])
  }
  invisible(x)
}

#' @export
summary.fastcsdid <- function(object, ...) {
  list(
    att_gt = data.frame(object$att_map, estimate = unname(object$att_gt),
                        supported = object$att_support$supported),
    aggregates = object$aggregation$metadata,
    regression = object$regression,
    influence = object$influence[c("cluster_att", "cluster_aggregate", "metadata")],
    vcov = object$vcov,
    inference_comparison = inference_comparison(object),
    jackknife = object$jackknife
  )
}
