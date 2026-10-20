.log1pexp <- function(x) {
  pmax(x, 0) + log1p(exp(-abs(x)))
}

.weighted_logit_state <- function(X, D, w, coefficients,
                                  keep_score_rows = TRUE) {
  eta <- as.numeric(X %*% coefficients)
  probability <- stats::plogis(eta)
  weighted_residual <- w * (D - probability)
  score <- as.numeric(crossprod(X, weighted_residual))
  score_rows <- if (keep_score_rows) X * weighted_residual else NULL
  information_weight <- w * probability * (1 - probability)
  list(
    eta = eta,
    probability = probability,
    loglik = sum(w * (D * eta - .log1pexp(eta))),
    score_rows = score_rows,
    score = score,
    information_weight = information_weight,
    information = crossprod(X, information_weight * X)
  )
}

.solve_logit_information <- function(information,
                                     tolerance = .Machine$double.eps) {
  p <- nrow(information)
  singular_values <- tryCatch(
    svd(information, nu = 0L, nv = 0L)$d,
    error = function(e) numeric()
  )
  if (!length(singular_values)) {
    return(list(valid = FALSE, rank = 0L, condition = 0,
                reason = "propensity information could not be decomposed"))
  }
  scale <- max(singular_values)
  rank <- if (scale > 0) sum(singular_values > scale * tolerance) else 0L
  condition <- tryCatch(rcond(information), error = function(e) 0)
  if (rank < p || !is.finite(condition) || condition < tolerance) {
    return(list(
      valid = FALSE, rank = as.integer(rank), condition = condition,
      reason = paste0("propensity information rank ", rank, " of ", p,
                      "; reciprocal condition number ",
                      format(condition, digits = 5))
    ))
  }
  list(valid = TRUE, rank = as.integer(rank), condition = condition,
       reason = "")
}

.solve_logit_information_fast <- function(information,
                                          tolerance = .Machine$double.eps) {
  condition <- tryCatch(rcond(information), error = function(e) 0)
  if (is.finite(condition) && condition >= tolerance) {
    R <- tryCatch(chol(information), error = function(e) NULL)
    if (!is.null(R)) {
      return(list(
        valid = TRUE, rank = as.integer(nrow(information)),
        condition = condition, reason = "", chol = R
      ))
    }
  }
  solved <- .solve_logit_information(information, tolerance)
  solved$chol <- if (isTRUE(solved$valid)) {
    tryCatch(chol(information), error = function(e) NULL)
  } else NULL
  solved
}

.fit_weighted_logit <- function(X, D, w, start = NULL, tolerance = 1e-10,
                                max_iterations = 100L,
                                information_tolerance = .Machine$double.eps) {
  X <- as.matrix(X)
  D <- as.numeric(D)
  w <- as.numeric(w)
  p <- ncol(X)
  blank <- list(
    coefficients = rep(NA_real_, p), converged = FALSE, iterations = 0L,
    score = rep(NA_real_, p), score_norm = Inf,
    information = matrix(NA_real_, p, p), loglik = NA_real_,
    probability = rep(NA_real_, nrow(X)), status = "nonconvergence",
    reason = "", rank = 0L, reciprocal_condition = 0,
    step_halvings = 0L, separation = FALSE
  )
  if (!nrow(X) || length(D) != nrow(X) || length(w) != nrow(X)) {
    blank$status <- "missing_propensity_sample"
    blank$reason <- "the propensity sample is empty or misaligned"
    return(blank)
  }
  if (any(!is.finite(X)) || any(!is.finite(D)) || any(!is.finite(w)) ||
      any(w < 0) || sum(w) <= 0) {
    blank$status <- "invalid_propensity_sample"
    blank$reason <- "the propensity sample contains invalid values or zero weight"
    return(blank)
  }
  if (!any(D == 1 & w > 0) || !any(D == 0 & w > 0)) {
    blank$status <- "missing_propensity_class"
    blank$reason <- "positive-weight treated and control observations are required"
    return(blank)
  }
  coefficients <- if (is.null(start)) rep(0, p) else as.numeric(start)
  if (length(coefficients) != p || any(!is.finite(coefficients))) {
    blank$status <- "invalid_propensity_start"
    blank$reason <- "the propensity starting vector is invalid"
    return(blank)
  }

  state <- .weighted_logit_state(X, D, w, coefficients,
                                 keep_score_rows = FALSE)
  total_halvings <- 0L
  converged <- FALSE
  last_step <- rep(Inf, p)
  for (iteration in seq_len(max_iterations)) {
    solved <- .solve_logit_information_fast(state$information,
                                            information_tolerance)
    if (!solved$valid) {
      blank$coefficients <- coefficients
      blank$iterations <- iteration - 1L
      blank$score <- state$score
      blank$score_norm <- max(abs(state$score)) / sum(w)
      blank$information <- state$information
      blank$loglik <- state$loglik
      blank$probability <- state$probability
      blank$status <- "singular_propensity_information"
      blank$reason <- solved$reason
      blank$rank <- solved$rank
      blank$reciprocal_condition <- solved$condition
      blank$step_halvings <- total_halvings
      return(blank)
    }
    step <- tryCatch({
      R <- solved$chol
      if (is.null(R)) R <- chol(state$information)
      backsolve(R, forwardsolve(t(R), state$score))
    }, error = function(e) NULL)
    if (is.null(step) || any(!is.finite(step))) {
      blank$status <- "singular_propensity_information"
      blank$reason <- "the propensity Newton system could not be solved"
      return(blank)
    }
    multiplier <- 1
    accepted <- FALSE
    for (halving in 0:30) {
      candidate <- coefficients + multiplier * step
      candidate_state <- .weighted_logit_state(
        X, D, w, candidate, keep_score_rows = FALSE
      )
      if (is.finite(candidate_state$loglik) &&
          candidate_state$loglik >= state$loglik -
            1e-12 * (1 + abs(state$loglik))) {
        accepted <- TRUE
        total_halvings <- total_halvings + halving
        break
      }
      multiplier <- multiplier / 2
    }
    if (!accepted) {
      blank$coefficients <- coefficients
      blank$iterations <- iteration - 1L
      blank$score <- state$score
      blank$score_norm <- max(abs(state$score)) / sum(w)
      blank$information <- state$information
      blank$loglik <- state$loglik
      blank$probability <- state$probability
      blank$status <- "propensity_line_search_failed"
      blank$reason <- "step halving could not increase the weighted log likelihood"
      blank$rank <- solved$rank
      blank$reciprocal_condition <- solved$condition
      blank$step_halvings <- total_halvings
      return(blank)
    }
    last_step <- multiplier * step
    coefficients <- candidate
    state <- candidate_state
    score_norm <- max(abs(state$score)) / sum(w)
    if (score_norm <= tolerance &&
        max(abs(last_step)) <= sqrt(tolerance) *
          (1 + max(abs(coefficients)))) {
      converged <- TRUE
      break
    }
  }
  solved <- .solve_logit_information_fast(state$information,
                                          information_tolerance)
  score_norm <- max(abs(state$score)) / sum(w)
  status <- if (converged) "valid" else "propensity_nonconvergence"
  reason <- if (converged) "" else
    paste0("weighted logit did not converge in ", max_iterations,
           " iterations; normalized score norm ",
           format(score_norm, digits = 5))
  separation <- max(abs(coefficients)) > 30 ||
    max(state$probability) > 1 - 1e-12 ||
    min(state$probability) < 1e-12
  list(
    coefficients = coefficients, converged = converged,
    iterations = as.integer(iteration), score = state$score,
    score_norm = score_norm, information = state$information,
    loglik = state$loglik, probability = state$probability,
    status = status, reason = reason, rank = solved$rank,
    reciprocal_condition = solved$condition,
    step_halvings = total_halvings, separation = separation
  )
}

.one_logit_newton_step <- function(X, D, w, coefficients) {
  state <- .weighted_logit_state(X, D, w, coefficients,
                                 keep_score_rows = FALSE)
  solved <- .solve_logit_information_fast(state$information)
  if (!solved$valid) return(coefficients)
  step <- tryCatch({
    R <- solved$chol
    if (is.null(R)) R <- chol(state$information)
    backsolve(R, forwardsolve(t(R), state$score))
  }, error = function(e) NULL)
  if (is.null(step) || any(!is.finite(step))) coefficients else
    coefficients + step
}

.propensity_overlap_status <- function(X, D) {
  fit <- .fit_weighted_logit(
    X, D, rep(1, length(D)), tolerance = 1e-8, max_iterations = 100L
  )
  maximum <- if (all(is.finite(fit$probability))) max(fit$probability) else
    NA_real_
  list(
    valid = fit$status == "valid" && is.finite(maximum) && maximum < 0.999,
    maximum = maximum,
    status = if (fit$status == "valid" && is.finite(maximum) && maximum >= 0.999)
      "overlap_violation" else fit$status,
    reason = if (fit$status == "valid" && is.finite(maximum) && maximum >= 0.999)
      paste0("unweighted overlap prefit has maximum propensity ",
             format(maximum, digits = 6), " (threshold 0.999)") else
        fit$reason,
    coefficients = fit$coefficients,
    information = fit$information,
    score = fit$score,
    iterations = fit$iterations
  )
}
