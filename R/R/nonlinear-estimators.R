.fit_treated_regression_blocks <- function(moments, att_map) {
  compilation_start <- proc.time()[["elapsed"]]
  regression_data <- moments$regression_data
  H <- length(moments$regions)
  keys <- unique(data.frame(
    group = rep(att_map$group, 2L),
    period = c(att_map$base, att_map$post)
  ))
  keys <- keys[order(keys$group, keys$period), , drop = FALSE]
  key_name <- paste(keys$group, keys$period, sep = "|")
  blocks <- lapply(seq_len(nrow(keys)), function(id) {
    regions <- which(moments$assignment == keys$group[id])
    period_index <- match(keys$period[id], moments$times)
    rows <- which(regression_data$region_index %in% regions &
                    regression_data$time_index == period_index)
    block <- .fit_one_regression_block(
      regression_data, rows, id, keys$period[id], regions, H
    )
    block$group <- keys$group[id]
    block
  })
  lookup <- stats::setNames(seq_along(blocks), key_name)
  cell_blocks <- cbind(
    base_block = unname(lookup[paste(att_map$group, att_map$base, sep = "|")]),
    post_block = unname(lookup[paste(att_map$group, att_map$post, sep = "|")])
  )
  metadata <- data.frame(
    block_id = seq_along(blocks), group = keys$group, period = keys$period,
    weight = vapply(blocks, `[[`, numeric(1), "total_weight"),
    rank = vapply(blocks, `[[`, integer(1), "rank"),
    status = vapply(blocks, `[[`, character(1), "status"),
    reason = vapply(blocks, `[[`, character(1), "reason"),
    stringsAsFactors = FALSE
  )
  list(
    blocks = blocks, metadata = metadata, cell_blocks = cell_blocks,
    timing = c(
      block_compilation = max(0, proc.time()[["elapsed"]] - compilation_start -
        sum(vapply(blocks, function(block) sum(block$timing), numeric(1)))),
      crossproducts = sum(vapply(blocks, function(block) {
        unname(block$timing["crossproducts"])
      }, numeric(1))),
      full_regression_solves = sum(vapply(blocks, function(block) {
        unname(block$timing["solve"])
      }, numeric(1)))
    )
  )
}

.weighted_ratio_with_gradient <- function(weight, value, X = NULL) {
  denominator <- sum(weight)
  if (!is.finite(denominator) || denominator <= 0) {
    return(list(value = NA_real_, denominator = denominator,
                gradient = if (is.null(X)) numeric() else rep(NA_real_, ncol(X))))
  }
  estimate <- sum(weight * value) / denominator
  gradient <- if (is.null(X)) numeric() else
    colSums(X * (weight * (value - estimate))) / denominator
  list(value = estimate, denominator = denominator, gradient = gradient)
}

.evaluate_nonlinear_sample <- function(block, regression_data, coefficients,
                                       method = c("ipw", "dr"),
                                       outcome_coefficients = NULL,
                                       keep_rows = rep(TRUE, length(block$rows))) {
  method <- match.arg(method)
  rows <- block$rows[keep_rows]
  if (!is.null(block$X) && !is.null(block$y) && !is.null(block$w)) {
    X <- block$X[keep_rows, , drop = FALSE]
    y <- block$y[keep_rows]
    w <- block$w[keep_rows]
  } else {
    X <- regression_data$X[rows, , drop = FALSE]
    y <- regression_data$y[rows]
    w <- regression_data$weight[rows]
  }
  D <- block$D[keep_rows]
  post <- block$post[keep_rows]
  if (!length(y)) {
    return(list(att = NA_real_, status = "missing_propensity_sample",
                reason = "the deleted ATT sample is empty"))
  }
  probability <- pmin(stats::plogis(as.numeric(X %*% coefficients)), 1-1e-6)
  trim <- D == 1 | probability < 0.995
  odds <- probability / (1 - probability)
  a10 <- trim * w * D * (1-post)
  a11 <- trim * w * D * post
  a00 <- trim * w * odds * (1-D) * (1-post)
  a01 <- trim * w * odds * (1-D) * post
  z <- y
  if (method == "dr") {
    required <- c("control_base", "control_post", "treated_base", "treated_post")
    if (is.null(outcome_coefficients) ||
        any(!vapply(outcome_coefficients[required], function(x) all(is.finite(x)),
                    logical(1)))) {
      return(list(att = NA_real_, status = "singular_outcome_regression",
                  reason = "a required DR outcome-regression block is invalid"))
    }
    m00 <- as.numeric(X %*% outcome_coefficients$control_base)
    m01 <- as.numeric(X %*% outcome_coefficients$control_post)
    m10 <- as.numeric(X %*% outcome_coefficients$treated_base)
    m11 <- as.numeric(X %*% outcome_coefficients$treated_post)
    z <- y - (post * m01 + (1-post) * m00)
  }
  r10 <- .weighted_ratio_with_gradient(a10, z)
  r11 <- .weighted_ratio_with_gradient(a11, z)
  r00 <- .weighted_ratio_with_gradient(a00, z, X)
  r01 <- .weighted_ratio_with_gradient(a01, z, X)
  ratios <- list(r10, r11, r00, r01)
  if (any(!vapply(ratios, function(x) is.finite(x$value), logical(1)))) {
    return(list(
      att = NA_real_, status = "zero_propensity_normalizer",
      reason = "a treated or propensity-weighted period normalizer is zero",
      normalizers = vapply(ratios, `[[`, numeric(1), "denominator")
    ))
  }
  att <- (r11$value-r10$value) - (r01$value-r00$value)
  gradient <- r00$gradient-r01$gradient
  correction <- 0
  if (method == "dr") {
    ad <- trim*w*D
    e_d_post <- .weighted_ratio_with_gradient(ad, m11-m01)
    e_1_post <- .weighted_ratio_with_gradient(a11, m11-m01)
    e_d_pre <- .weighted_ratio_with_gradient(ad, m10-m00)
    e_0_pre <- .weighted_ratio_with_gradient(a10, m10-m00)
    correction_parts <- list(e_d_post, e_1_post, e_d_pre, e_0_pre)
    if (any(!vapply(correction_parts, function(x) is.finite(x$value),
                    logical(1)))) {
      return(list(att = NA_real_, status = "zero_propensity_normalizer",
                  reason = "a locally efficient DR normalizer is zero"))
    }
    correction <- e_d_post$value-e_1_post$value-e_d_pre$value+e_0_pre$value
    att <- att + correction
  }
  list(
    att = att, status = "valid", reason = "", gradient = gradient,
    normalizers = c(treated_base=r10$denominator,
                    treated_post=r11$denominator,
                    control_base=r00$denominator,
                    control_post=r01$denominator),
    probability = probability, trim = trim, correction = correction
  )
}

.nonlinear_outcome_coefficients <- function(k, control_nuisance,
                                            treated_nuisance) {
  cb <- control_nuisance$cell_blocks[k, "base_block"]
  cp <- control_nuisance$cell_blocks[k, "post_block"]
  tb <- treated_nuisance$cell_blocks[k, "base_block"]
  tp <- treated_nuisance$cell_blocks[k, "post_block"]
  list(
    control_base = control_nuisance$blocks[[cb]]$coefficients,
    control_post = control_nuisance$blocks[[cp]]$coefficients,
    treated_base = treated_nuisance$blocks[[tb]]$coefficients,
    treated_post = treated_nuisance$blocks[[tp]]$coefficients
  )
}

.evaluate_nonlinear_blocks <- function(moments, att_map, propensity,
                                       method = c("ipw", "dr"),
                                       control_nuisance = NULL,
                                       treated_nuisance = NULL) {
  method <- match.arg(method)
  K <- nrow(att_map)
  att <- rep(NA_real_, K)
  evaluations <- vector("list", K)
  support <- data.frame(
    supported = FALSE, status = "unsupported", reason = "",
    propensity_block = seq_len(K), stringsAsFactors = FALSE
  )
  for (k in seq_len(K)) {
    block <- propensity$blocks[[k]]
    if (block$fit$status != "valid") {
      support$status[k] <- block$fit$status
      support$reason[k] <- block$fit$reason
      next
    }
    outcome <- if (method == "dr") .nonlinear_outcome_coefficients(
      k, control_nuisance, treated_nuisance
    ) else NULL
    evaluated <- .evaluate_nonlinear_sample(
      block, moments$regression_data, block$coefficients, method, outcome
    )
    evaluations[[k]] <- evaluated
    att[k] <- evaluated$att
    support$supported[k] <- evaluated$status == "valid"
    support$status[k] <- evaluated$status
    support$reason[k] <- evaluated$reason
  }
  names(att) <- .att_names(att_map)
  groups <- sort(unique(att_map$group))
  group_size <- stats::setNames(vapply(groups, function(group) {
    sum(moments$region_total_weight[moments$assignment == group])
  }, numeric(1)), as.character(groups))
  list(att = att, support = support, evaluations = evaluations,
       group_size = group_size, propensity = propensity,
       control_nuisance = control_nuisance,
       treated_nuisance = treated_nuisance)
}
