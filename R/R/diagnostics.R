.validate_jackknife_fit <- function(fit) {
  if (!inherits(fit, "fastcsdid")) {
    stop("fit must be returned by fast_csdid() or fit_csdid_moments().",
         call. = FALSE)
  }
  if (is.null(fit$jackknife)) {
    stop("fit does not contain jackknife results; use inference = 'jackknife'.",
         call. = FALSE)
  }
}

.support_status_vocabulary <- function() {
  data.frame(
    label = c(
      "valid", "unaffected", "updated", "missing_base_period",
      "missing_post_period", "cohort_removed", "control_group_removed",
      "zero_weight_treated_cell", "zero_weight_control_cell",
      "zero_weight_regression_sample", "singular_outcome_regression",
      "rank_changed", "missing_regression_sample",
      "undefined_prediction_average", "undefined_aggregation",
      "missing_propensity_sample", "invalid_propensity_sample",
      "missing_propensity_class", "singular_propensity_information",
      "singular_deleted_information", "propensity_line_search_failed",
      "propensity_nonconvergence", "propensity_separation",
      "overlap_violation", "zero_propensity_normalizer"
    ),
    meaning = c(
      "The full-sample estimand has all required support.",
      "Deleting the cluster does not enter the estimand.",
      "The delete-one estimate is supported and recomputed.",
      "A treated cohort has no observed pre-treatment base period.",
      "A treated cohort has no observed post-treatment period.",
      "The deletion removes all sampling-weight mass for a treated cohort.",
      "The deletion removes all sampling-weight mass for the required control group.",
      "A required treated base or post cell has zero sampling-weight mass.",
      "A required control base or post cell has zero sampling-weight mass.",
      "An outcome-regression block has zero positive sampling-weight mass.",
      "An outcome-regression block is singular or numerically ill-conditioned.",
      "Cluster deletion changes the identified outcome-regression rank.",
      "The required outcome-regression estimation sample is absent.",
      "The treated covariate prediction average is undefined.",
      "A required ATT is unsupported or the aggregate weighting mass is zero.",
      "The ATT-specific propensity sample is absent.",
      "The propensity sample contains invalid values or zero total weight.",
      "The propensity sample lacks a positive-weight treated or control class.",
      "The full-sample propensity information is singular or ill-conditioned.",
      "The deleted propensity information is singular or ill-conditioned.",
      "Step halving could not improve the weighted log likelihood.",
      "The ATT-specific propensity fit did not converge.",
      "The fitted propensity model has a separation diagnostic.",
      "The did-compatible overlap prefit reaches the 0.999 threshold.",
      "A normalized treated or propensity-weighted moment has zero mass."
    ),
    stringsAsFactors = FALSE
  )
}

.validate_influence_options <- function(relative_threshold, level) {
  if (length(relative_threshold) != 1L || !is.finite(relative_threshold) ||
      relative_threshold < 0) {
    stop("relative_threshold must be one nonnegative finite number.", call. = FALSE)
  }
  if (length(level) != 1L || !is.finite(level) || level <= 0 || level >= 1) {
    stop("level must be strictly between zero and one.", call. = FALSE)
  }
}

.jackknife_estimand <- function(fit, estimand) {
  if (length(estimand) != 1L || is.na(estimand)) {
    stop("estimand must be one ATT or aggregate label.", call. = FALSE)
  }
  estimand <- as.character(estimand)
  att_position <- match(estimand, names(fit$att_gt))
  aggregate_position <- match(estimand, names(fit$aggregation$estimates))
  if (!is.na(att_position)) {
    return(list(
      type = "att", position = att_position,
      full = unname(fit$att_gt[att_position]),
      deleted = fit$jackknife$att[, att_position],
      status = fit$jackknife$att_status[, att_position],
      status_label = fit$jackknife$att_status_label[, att_position],
      reason = fit$jackknife$att_reason[, att_position],
      se = unname(fit$jackknife$att_covariance$pairwise_se_cv3[att_position]),
      component_positions = att_position
    ))
  }
  if (!is.na(aggregate_position)) {
    return(list(
      type = "aggregate", position = aggregate_position,
      full = unname(fit$aggregation$estimates[aggregate_position]),
      deleted = fit$jackknife$aggregate[, aggregate_position],
      status = fit$jackknife$aggregate_status[, aggregate_position],
      status_label = fit$jackknife$aggregate_status_label[, aggregate_position],
      reason = fit$jackknife$aggregate_reason[, aggregate_position],
      se = unname(fit$jackknife$aggregate_covariance$pairwise_se_cv3[aggregate_position]),
      component_positions = which(
        is.finite(fit$aggregation$weights[aggregate_position, ]) &
          fit$aggregation$weights[aggregate_position, ] != 0
      )
    ))
  }
  stop("unknown estimand label. Use names(fit$att_gt) or names(fit$aggregation$estimates).",
       call. = FALSE)
}

#' Inspect delete-one-cluster estimates for one estimand
#'
#' @param fit A fitted `fastcsdid` object containing jackknife results.
#' @param estimand One label from `names(fit$att_gt)` or
#'   `names(fit$aggregation$estimates)`.
#' @param relative_threshold A full estimate with absolute value no larger than
#'   this fraction of the observed estimate scale is treated as near zero, and
#'   relative changes are reported as `NA`.
#' @param level Confidence level for the fixed-standard-error significance-flip
#'   diagnostic.
#' @return A data frame with one row per deleted cluster. Significance flips use
#'   the full-sample pairwise-CV3 standard error for both the full and deleted
#'   point estimates; they are sensitivity diagnostics, not deletion-specific
#'   hypothesis tests.
#' @export
jackknife_cluster_table <- function(fit, estimand, relative_threshold = 1e-8,
                                    level = 0.95) {
  .validate_jackknife_fit(fit)
  .validate_influence_options(relative_threshold, level)
  resolved <- .jackknife_estimand(fit, estimand)
  deleted <- as.numeric(resolved$deleted)
  change <- deleted - resolved$full
  finite <- is.finite(deleted) & is.finite(resolved$full)
  estimate_scale <- max(c(1, abs(resolved$full), abs(deleted[is.finite(deleted)])))
  denominator_ok <- is.finite(resolved$full) &&
    abs(resolved$full) > relative_threshold * estimate_scale
  relative_change <- if (denominator_ok) change / abs(resolved$full) else
    rep(NA_real_, length(change))
  sign_flip <- rep(NA, length(change))
  sign_comparable <- finite & resolved$full != 0 & deleted != 0
  sign_flip[sign_comparable] <- sign(deleted[sign_comparable]) != sign(resolved$full)

  significance_flip <- rep(NA, length(change))
  if (is.finite(resolved$se) && resolved$se > 0 && is.finite(resolved$full)) {
    critical <- stats::qnorm(1 - (1 - level) / 2) * resolved$se
    significance_flip[finite] <-
      (abs(deleted[finite]) > critical) != (abs(resolved$full) > critical)
  }
  out <- data.frame(
    estimand_type = resolved$type,
    estimand = as.character(estimand),
    cluster = rownames(as.matrix(resolved$deleted)),
    full_estimate = resolved$full,
    delete_estimate = deleted,
    change = change,
    relative_change = relative_change,
    sign_flip = sign_flip,
    significance_flip = significance_flip,
    valid = finite,
    status_code = as.integer(resolved$status),
    status = as.character(resolved$status_label),
    reason = as.character(resolved$reason),
    jackknife_se = resolved$se,
    stringsAsFactors = FALSE
  )
  nuisance <- fit$jackknife$nuisance_diagnostics
  out$rank_changed <- NA
  out$failed_nuisance_block <- NA_character_
  out$max_coefficient_change <- NA_real_
  out$influence_source <- NA_character_
  out$jackknife_method <- if (is.null(fit$settings$jackknife_method))
    "exact" else fit$settings$jackknife_method
  out$max_propensity_coefficient_change <- NA_real_
  out$deleted_score_norm <- NA_real_
  out$deleted_information_condition <- NA_real_
  out$failed_propensity_block <- NA_character_
  out$propensity_status <- NA_character_
  out$direct_moment_change <- NA_real_
  out$outcome_regression_change <- NA_real_
  out$propensity_change <- NA_real_
  out$cluster_score_norm <- NA_real_
  out$propensity_score_min <- NA_real_
  out$propensity_score_max <- NA_real_
  out$effective_treated_weight <- NA_real_
  out$effective_comparison_weight <- NA_real_
  out$refit_delete_estimate <- NA_real_
  out$refit_valid <- NA
  out$refit_status <- NA_character_
  out$refit_reason <- NA_character_
  out$linearized_refit_difference <- NA_real_
  out$linearized_refit_sign_change <- NA
  out$linearized_refit_significance_change <- NA
  if (!is.null(nuisance$att_rank_changed) &&
      length(resolved$component_positions)) {
    components <- resolved$component_positions
    out$rank_changed <- apply(
      nuisance$att_rank_changed[, components, drop = FALSE], 1L, any
    )
    out$failed_nuisance_block <- apply(
      nuisance$att_failed_block[, components, drop = FALSE], 1L,
      function(values) paste(sort(unique(values[nzchar(values)])), collapse = ",")
    )
    out$max_coefficient_change <- apply(
      nuisance$att_max_coefficient_change[, components, drop = FALSE], 1L,
      function(values) {
        values <- values[is.finite(values)]
        if (length(values)) max(values) else NA_real_
      }
    )
    out$influence_source <- apply(
      nuisance$att_influence_source[, components, drop = FALSE], 1L,
      function(values) {
        values <- sort(unique(values[values != "unaffected"]))
        if (length(values)) paste(values, collapse = ",") else "unaffected"
      }
    )
  }
  if (!is.null(nuisance$propensity_coefficient_change) &&
      length(resolved$component_positions)) {
    components <- resolved$component_positions
    row_max <- function(matrix_value) apply(
      matrix_value[, components, drop = FALSE], 1L,
      function(values) {
        values <- values[is.finite(values)]
        if (length(values)) max(abs(values)) else NA_real_
      }
    )
    row_labels <- function(matrix_value, exclude = character()) apply(
      matrix_value[, components, drop = FALSE], 1L,
      function(values) {
        values <- sort(unique(values[
          !is.na(values) & nzchar(values) & !(values %in% exclude)
        ]))
        if (length(values)) paste(values, collapse = ",") else ""
      }
    )
    out$max_propensity_coefficient_change <- row_max(
      nuisance$propensity_coefficient_change
    )
    out$deleted_score_norm <- row_max(nuisance$deleted_score_norm)
    out$deleted_information_condition <- row_max(
      nuisance$deleted_information_condition
    )
    out$failed_propensity_block <- row_labels(
      nuisance$failed_propensity_block
    )
    out$propensity_status <- row_labels(
      nuisance$propensity_status, exclude = "unaffected"
    )
    out$direct_moment_change <- row_max(nuisance$direct_moment_change)
    out$outcome_regression_change <- row_max(
      nuisance$outcome_regression_change
    )
    out$propensity_change <- row_max(nuisance$propensity_change)
    factorized <- fit$influence$factorized$evaluation$propensity
    probabilities <- unlist(lapply(components, function(k) {
      factorized$blocks[[k]]$probability
    }), use.names = FALSE)
    out$propensity_score_min <- min(probabilities, na.rm = TRUE)
    out$propensity_score_max <- max(probabilities, na.rm = TRUE)
    out$effective_treated_weight <- min(
      factorized$metadata$treated_weight[components], na.rm = TRUE
    )
    out$effective_comparison_weight <- min(
      factorized$metadata$control_weight[components], na.rm = TRUE
    )
    out$cluster_score_norm <- apply(
      do.call(cbind, lapply(components, function(k) {
        sqrt(rowSums(factorized$blocks[[k]]$cluster_score^2))
      })), 1L, max
    )
  }
  if (!is.null(fit$jackknife$methods$refit) &&
      !is.null(fit$jackknife$methods$linearized)) {
    refit <- fit$jackknife$methods$refit
    refit_deleted <- if (resolved$type == "att") {
      refit$att[, resolved$position]
    } else {
      refit$aggregate[, resolved$position]
    }
    refit_status <- if (resolved$type == "att") {
      refit$att_status_label[, resolved$position]
    } else {
      refit$aggregate_status_label[, resolved$position]
    }
    refit_reason <- if (resolved$type == "att") {
      refit$att_reason[, resolved$position]
    } else {
      refit$aggregate_reason[, resolved$position]
    }
    out$refit_delete_estimate <- as.numeric(refit_deleted)
    out$refit_valid <- is.finite(refit_deleted)
    out$refit_status <- as.character(refit_status)
    out$refit_reason <- as.character(refit_reason)
    out$linearized_refit_difference <- deleted - as.numeric(refit_deleted)
    comparable <- is.finite(deleted) & is.finite(refit_deleted) &
      deleted != 0 & refit_deleted != 0
    out$linearized_refit_sign_change[comparable] <-
      sign(deleted[comparable]) != sign(refit_deleted[comparable])
    refit_se <- if (resolved$type == "att") {
      refit$att_covariance$pairwise_se_cv3[resolved$position]
    } else {
      refit$aggregate_covariance$pairwise_se_cv3[resolved$position]
    }
    if (is.finite(resolved$se) && resolved$se > 0 &&
        is.finite(refit_se) && refit_se > 0) {
      critical <- stats::qnorm(1 - (1 - level) / 2)
      out$linearized_refit_significance_change[comparable] <-
        (abs(deleted[comparable]) / resolved$se > critical) !=
        (abs(refit_deleted[comparable]) / refit_se > critical)
    }
  }
  out
}

#' Summarize cluster influence across ATT and aggregate estimands
#'
#' @inheritParams jackknife_cluster_table
#' @return A compact data frame reporting extrema of the deletion estimates,
#'   the largest absolute and relative changes, the most influential cluster,
#'   sign and fixed-standard-error significance flips, supported-deletion
#'   counts, pairwise-CV3 standard errors, and failure reasons.
#' @export
jackknife_influence_summary <- function(fit, relative_threshold = 1e-8,
                                        level = 0.95) {
  .validate_jackknife_fit(fit)
  .validate_influence_options(relative_threshold, level)
  labels <- c(names(fit$att_gt), names(fit$aggregation$estimates))
  rows <- lapply(labels, function(label) {
    table <- jackknife_cluster_table(fit, label, relative_threshold, level)
    valid <- table$valid
    valid_change <- which(valid & is.finite(table$change))
    influential <- if (length(valid_change)) {
      table$cluster[valid_change[which.max(abs(table$change[valid_change]))]]
    } else {
      NA_character_
    }
    failures <- sort(unique(table$reason[!valid & nzchar(table$reason)]))
    data.frame(
      estimand_type = table$estimand_type[1L],
      estimand = label,
      estimate = table$full_estimate[1L],
      min_delete = if (any(valid)) min(table$delete_estimate[valid]) else NA_real_,
      max_delete = if (any(valid)) max(table$delete_estimate[valid]) else NA_real_,
      max_abs_change = if (length(valid_change))
        max(abs(table$change[valid_change])) else NA_real_,
      most_influential_cluster = influential,
      largest_relative_change = if (any(is.finite(table$relative_change)))
        max(abs(table$relative_change), na.rm = TRUE) else NA_real_,
      sign_flip = if (any(!is.na(table$sign_flip)))
        any(table$sign_flip, na.rm = TRUE) else NA,
      significance_flip = if (any(!is.na(table$significance_flip)))
        any(table$significance_flip, na.rm = TRUE) else NA,
      n_valid = sum(valid),
      n_failed = sum(!valid),
      jackknife_se = table$jackknife_se[1L],
      rank_changed = if (any(!is.na(table$rank_changed)))
        any(table$rank_changed, na.rm = TRUE) else NA,
      failed_nuisance_blocks = paste(
        sort(unique(table$failed_nuisance_block[
          !is.na(table$failed_nuisance_block) & nzchar(table$failed_nuisance_block)
        ])), collapse = ","
      ),
      max_coefficient_change = if (any(is.finite(table$max_coefficient_change)))
        max(table$max_coefficient_change, na.rm = TRUE) else NA_real_,
      influence_sources = paste(
        sort(unique(table$influence_source[
          !is.na(table$influence_source) & table$influence_source != "unaffected"
        ])), collapse = ","
      ),
      jackknife_method = table$jackknife_method[1L],
      max_propensity_coefficient_change =
        if (any(is.finite(table$max_propensity_coefficient_change)))
          max(table$max_propensity_coefficient_change, na.rm = TRUE) else
            NA_real_,
      max_deleted_score_norm = if (any(is.finite(table$deleted_score_norm)))
        max(table$deleted_score_norm, na.rm = TRUE) else NA_real_,
      min_deleted_information_condition =
        if (any(is.finite(table$deleted_information_condition)))
          min(table$deleted_information_condition, na.rm = TRUE) else NA_real_,
      failed_propensity_blocks = paste(sort(unique(
        table$failed_propensity_block[
          !is.na(table$failed_propensity_block) &
            nzchar(table$failed_propensity_block)
        ]
      )), collapse = ","),
      propensity_statuses = paste(sort(unique(
        table$propensity_status[
          !is.na(table$propensity_status) & nzchar(table$propensity_status)
        ]
      )), collapse = ","),
      max_cluster_score_norm = if (any(is.finite(table$cluster_score_norm)))
        max(table$cluster_score_norm, na.rm = TRUE) else NA_real_,
      propensity_score_min = if (any(is.finite(table$propensity_score_min)))
        min(table$propensity_score_min, na.rm = TRUE) else NA_real_,
      propensity_score_max = if (any(is.finite(table$propensity_score_max)))
        max(table$propensity_score_max, na.rm = TRUE) else NA_real_,
      effective_treated_weight = table$effective_treated_weight[1L],
      effective_comparison_weight = table$effective_comparison_weight[1L],
      max_abs_refit_change = if (any(is.finite(table$refit_delete_estimate)))
        max(abs(table$refit_delete_estimate - table$full_estimate),
            na.rm = TRUE) else NA_real_,
      n_refit_failures = if (any(!is.na(table$refit_valid)))
        sum(!table$refit_valid, na.rm = TRUE) else NA_integer_,
      refit_failure_reasons = paste(sort(unique(
        table$refit_reason[
          !is.na(table$refit_reason) & nzchar(table$refit_reason) &
            table$refit_valid %in% FALSE
        ]
      )), collapse = "; "),
      max_abs_linearized_refit_difference =
        if (any(is.finite(table$linearized_refit_difference)))
          max(abs(table$linearized_refit_difference), na.rm = TRUE) else
            NA_real_,
      max_discrepancy_cluster = if (any(is.finite(
        table$linearized_refit_difference
      ))) table$cluster[which.max(abs(table$linearized_refit_difference))] else
        NA_character_,
      linearized_refit_sign_change = if (any(!is.na(
        table$linearized_refit_sign_change
      ))) any(table$linearized_refit_sign_change, na.rm = TRUE) else NA,
      linearized_refit_significance_change = if (any(!is.na(
        table$linearized_refit_significance_change
      ))) any(table$linearized_refit_significance_change, na.rm = TRUE) else NA,
      failure_reasons = paste(failures, collapse = "; "),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  attr(out, "relative_threshold") <- relative_threshold
  attr(out, "level") <- level
  attr(out, "significance_diagnostic") <-
    "Full and delete-one estimates are compared using the same full-sample pairwise-CV3 standard error."
  out
}
