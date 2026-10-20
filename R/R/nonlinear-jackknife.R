.deleted_propensity_update <- function(block, cluster) {
  information <- block$fit$information -
    block$cluster_information[, , cluster]
  solved <- .solve_logit_information(information)
  blank <- list(
    coefficients = rep(NA_real_, length(block$coefficients)),
    change = rep(NA_real_, length(block$coefficients)),
    information = information, status = "singular_deleted_information",
    reason = solved$reason, rank = solved$rank,
    reciprocal_condition = solved$condition
  )
  if (!solved$valid) return(blank)
  change <- tryCatch({
    R <- chol(information)
    -backsolve(R, forwardsolve(t(R), block$cluster_score[cluster, ]))
  }, error = function(e) NULL)
  if (is.null(change) || any(!is.finite(change))) {
    blank$reason <- "the deleted propensity Newton system could not be solved"
    return(blank)
  }
  list(
    coefficients = block$coefficients + change,
    change = as.numeric(change), information = information,
    status = "valid", reason = "", rank = solved$rank,
    reciprocal_condition = solved$condition
  )
}

.deleted_outcome_coefficients <- function(k, deleted_control,
                                          deleted_treated,
                                          control_nuisance,
                                          treated_nuisance) {
  ids <- c(
    control_base = unname(control_nuisance$cell_blocks[k, "base_block"]),
    control_post = unname(control_nuisance$cell_blocks[k, "post_block"]),
    treated_base = unname(treated_nuisance$cell_blocks[k, "base_block"]),
    treated_post = unname(treated_nuisance$cell_blocks[k, "post_block"])
  )
  blocks <- list(
    control_base = deleted_control[[as.integer(ids["control_base"])]],
    control_post = deleted_control[[as.integer(ids["control_post"])]],
    treated_base = deleted_treated[[as.integer(ids["treated_base"])]],
    treated_post = deleted_treated[[as.integer(ids["treated_post"])]]
  )
  invalid <- vapply(blocks, `[[`, character(1), "status") != "valid"
  list(
    coefficients = lapply(blocks, `[[`, "coefficients"),
    valid = !any(invalid),
    rank_changed = any(vapply(blocks, `[[`, logical(1), "rank_changed")),
    failed = names(blocks)[invalid],
    reason = paste(unique(vapply(blocks[invalid], `[[`, character(1),
                                 "reason")), collapse = "; ")
  )
}

.deleted_propensity_score_diagnostics <- function(block, regression_data,
                                                  coefficients, keep) {
  if (!is.null(block$X) && !is.null(block$w)) {
    X <- block$X[keep, , drop = FALSE]
    w <- block$w[keep]
  } else {
    X <- regression_data$X[block$rows[keep], , drop = FALSE]
    w <- regression_data$weight[block$rows[keep]]
  }
  state <- .weighted_logit_state(X, block$D[keep], w, coefficients)
  list(score_norm = max(abs(state$score)) / sum(w),
       loglik = state$loglik)
}

.deleted_propensity_refit <- function(block, cluster, coefficients,
                                      regression_data) {
  if (!is.null(block$X) && !is.null(block$w)) {
    keep <- block$cluster != cluster
    return(.fit_weighted_logit(
      block$X[keep, , drop = FALSE], block$D[keep], block$w[keep],
      start = coefficients
    ))
  }
  keep <- block$cluster != cluster
  .fit_weighted_logit(
    regression_data$X[block$rows[keep], , drop = FALSE],
    block$D[keep],
    regression_data$weight[block$rows[keep]],
    start = coefficients
  )
}

.deleted_overlap_status <- function(block, regression_data, keep, cluster) {
  if (!is.null(block$X)) {
    start <- block$overlap$coefficients
    if (is.null(start) || any(!is.finite(start))) {
      return(.propensity_overlap_status(
        block$X[keep, , drop = FALSE], block$D[keep]
      ))
    }
    if (!is.null(block$overlap_cluster_information) &&
        all(is.finite(start))) {
      deleted_information <- block$overlap$information -
        block$overlap_cluster_information[, , cluster]
      solved <- .solve_logit_information(deleted_information)
      if (solved$valid) {
        change <- tryCatch({
          R <- chol(deleted_information)
          -backsolve(
            R,
            forwardsolve(t(R), block$overlap_cluster_score[cluster, ])
          )
        }, error = function(e) NULL)
        if (!is.null(change) && all(is.finite(change))) {
          start <- start + as.numeric(change)
        }
      }
    }
    fit <- .fit_weighted_logit(
      block$X[keep, , drop = FALSE], block$D[keep], block$unit_w[keep],
      start = start, tolerance = 1e-8, max_iterations = 100L
    )
    maximum <- if (all(is.finite(fit$probability))) max(fit$probability) else
      NA_real_
    return(list(
      valid = fit$status == "valid" && is.finite(maximum) && maximum < 0.999,
      maximum = maximum,
      status = if (fit$status == "valid" && is.finite(maximum) &&
                   maximum >= 0.999) {
        "overlap_violation"
      } else fit$status,
      reason = if (fit$status == "valid" && is.finite(maximum) &&
                   maximum >= 0.999) {
        paste0("unweighted overlap prefit has maximum propensity ",
               format(maximum, digits = 6), " (threshold 0.999)")
      } else fit$reason,
      iterations = fit$iterations
    ))
  }
  .propensity_overlap_status(
    regression_data$X[block$rows[keep], , drop = FALSE], block$D[keep]
  )
}

.delete_nonlinear_att_matrix <- function(moments, att_map, evaluation,
                                         method = c("ipw", "dr"),
                                         jackknife_method = c("linearized",
                                                              "refit")) {
  method <- match.arg(method)
  jackknife_method <- match.arg(jackknife_method)
  deletion_start <- proc.time()[["elapsed"]]
  propensity <- evaluation$propensity
  H <- length(moments$regions)
  K <- nrow(att_map)
  full_att <- evaluation$att
  estimates <- matrix(rep(full_att, each = H), H, K)
  status <- matrix(0L, H, K)
  status_label <- matrix("unaffected", H, K)
  reason <- matrix("unaffected", H, K)
  coefficient_change <- matrix(0, H, K)
  deleted_score_norm <- matrix(0, H, K)
  deleted_condition <- matrix(NA_real_, H, K)
  failed_propensity_block <- matrix("", H, K)
  propensity_status <- matrix("unaffected", H, K)
  direct_change <- matrix(0, H, K)
  outcome_change <- matrix(0, H, K)
  propensity_change <- matrix(0, H, K)
  refit_iterations <- matrix(0L, H, K)

  labels <- list(as.character(moments$regions), .att_names(att_map))
  for (object in c("estimates", "status", "status_label", "reason",
                   "coefficient_change", "deleted_score_norm",
                   "deleted_condition", "failed_propensity_block",
                   "propensity_status", "direct_change", "outcome_change",
                   "propensity_change", "refit_iterations")) {
    target <- get(object)
    dimnames(target) <- labels
    assign(object, target)
  }

  unsupported <- !is.finite(full_att)
  if (any(unsupported)) {
    status[, unsupported] <- 2L
    status_label[, unsupported] <- matrix(
      evaluation$support$status[unsupported], H, sum(unsupported), byrow = TRUE
    )
    reason[, unsupported] <- matrix(
      paste0("full-sample ATT is unsupported: ",
             evaluation$support$reason[unsupported]),
      H, sum(unsupported), byrow = TRUE
    )
  }

  for (h in seq_len(H)) {
    deleted_control <- deleted_treated <- NULL
    if (method == "dr") {
      deleted_control <- .delete_nuisance_blocks(
        evaluation$control_nuisance, h
      )
      deleted_treated <- .delete_nuisance_blocks(
        evaluation$treated_nuisance, h
      )
    }
    for (k in seq_len(K)) {
      block <- propensity$blocks[[k]]
      if (!is.finite(full_att[k]) || block$cluster_weight[h] <= 0) next
      keep <- block$cluster != h
      status[h, k] <- 1L
      status_label[h, k] <- "updated"
      reason[h, k] <- "updated"

      propensity_update <- .deleted_propensity_update(block, h)
      deleted_condition[h, k] <- propensity_update$reciprocal_condition
      if (propensity_update$status != "valid") {
        estimates[h, k] <- NA_real_
        status[h, k] <- 2L
        status_label[h, k] <- propensity_update$status
        reason[h, k] <- propensity_update$reason
        failed_propensity_block[h, k] <- paste0("block:", k)
        propensity_status[h, k] <- propensity_update$status
        next
      }
      coefficient_change[h, k] <- max(abs(propensity_update$change))

      full_outcome <- deleted_outcome <- NULL
      if (method == "dr") {
        full_outcome <- .nonlinear_outcome_coefficients(
          k, evaluation$control_nuisance, evaluation$treated_nuisance
        )
        deleted <- .deleted_outcome_coefficients(
          k, deleted_control, deleted_treated,
          evaluation$control_nuisance, evaluation$treated_nuisance
        )
        if (!deleted$valid) {
          estimates[h, k] <- NA_real_
          status[h, k] <- 2L
          status_label[h, k] <- if (deleted$rank_changed) "rank_changed" else
            "singular_outcome_regression"
          reason[h, k] <- deleted$reason
          next
        }
        deleted_outcome <- deleted$coefficients
      }

      direct <- .evaluate_nonlinear_sample(
        block, moments$regression_data, block$coefficients, method,
        full_outcome, keep
      )
      with_outcome <- if (method == "dr") .evaluate_nonlinear_sample(
        block, moments$regression_data, block$coefficients, method,
        deleted_outcome, keep
      ) else direct
      if (direct$status != "valid" || with_outcome$status != "valid") {
        failed <- if (direct$status != "valid") direct else with_outcome
        estimates[h, k] <- NA_real_
        status[h, k] <- 2L
        status_label[h, k] <- failed$status
        reason[h, k] <- failed$reason
        next
      }
      direct_change[h, k] <- direct$att - full_att[k]
      outcome_change[h, k] <- with_outcome$att - direct$att

      if (jackknife_method == "linearized") {
        X_deleted <- moments$regression_data$X[
          block$rows[keep], , drop = FALSE
        ]
        linearized_probability <- stats::plogis(as.numeric(
          X_deleted %*% propensity_update$coefficients
        ))
        if (max(linearized_probability) >= 0.999) {
          estimates[h, k] <- NA_real_
          status[h, k] <- 2L
          status_label[h, k] <- "overlap_violation"
          reason[h, k] <- paste0(
            "linearized deleted propensity reaches ",
            format(max(linearized_probability), digits = 6),
            " (threshold 0.999)"
          )
          failed_propensity_block[h, k] <- paste0("block:", k)
          propensity_status[h, k] <- "overlap_violation"
          next
        }
        propensity_effect <- as.numeric(crossprod(
          with_outcome$gradient, propensity_update$change
        ))
        estimate <- with_outcome$att + propensity_effect
        propensity_status[h, k] <- "linearized"
        score_diagnostic <- .deleted_propensity_score_diagnostics(
          block, moments$regression_data, propensity_update$coefficients, keep
        )
      } else {
        refit <- .deleted_propensity_refit(
          block, h, propensity_update$coefficients, moments$regression_data
        )
        overlap <- .deleted_overlap_status(
          block, moments$regression_data, keep, h
        )
        if (refit$status == "valid" && !overlap$valid) {
          refit$converged <- FALSE
          refit$status <- overlap$status
          refit$reason <- overlap$reason
        }
        propensity_status[h, k] <- refit$status
        refit_iterations[h, k] <- refit$iterations
        if (refit$status != "valid") {
          estimates[h, k] <- NA_real_
          status[h, k] <- 2L
          status_label[h, k] <- refit$status
          reason[h, k] <- refit$reason
          failed_propensity_block[h, k] <- paste0("block:", k)
          deleted_score_norm[h, k] <- refit$score_norm
          next
        }
        exact <- .evaluate_nonlinear_sample(
          block, moments$regression_data, refit$coefficients, method,
          deleted_outcome, keep
        )
        if (exact$status != "valid") {
          estimates[h, k] <- NA_real_
          status[h, k] <- 2L
          status_label[h, k] <- exact$status
          reason[h, k] <- exact$reason
          next
        }
        estimate <- exact$att
        coefficient_change[h, k] <- max(abs(
          refit$coefficients - block$coefficients
        ))
        score_diagnostic <- list(score_norm = refit$score_norm,
                                 loglik = refit$loglik)
      }
      estimates[h, k] <- estimate
      deleted_score_norm[h, k] <- score_diagnostic$score_norm
      propensity_change[h, k] <- estimate - with_outcome$att
    }
  }
  list(
    estimates = estimates,
    changes = sweep(estimates, 2L, full_att, `-`),
    status = status, status_label = status_label, reason = reason,
    status_legend = c(unaffected = 0L, updated = 1L, unsupported = 2L),
    influence = sweep(matrix(full_att, H, K, byrow = TRUE), 2L,
                      rep(0, K), `+`) - estimates,
    nuisance_diagnostics = list(
      jackknife_method = jackknife_method,
      propensity_coefficient_change = coefficient_change,
      deleted_score_norm = deleted_score_norm,
      deleted_information_condition = deleted_condition,
      failed_propensity_block = failed_propensity_block,
      propensity_status = propensity_status,
      direct_moment_change = direct_change,
      outcome_regression_change = outcome_change,
      propensity_change = propensity_change,
      refit_iterations = refit_iterations,
      affected_cluster_blocks = sum(status != 0L),
      elapsed_seconds = proc.time()[["elapsed"]] - deletion_start
    )
  )
}
