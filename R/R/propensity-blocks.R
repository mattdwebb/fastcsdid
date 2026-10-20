.propensity_sample_rows <- function(regression_data, moments, cell,
                                    control_group) {
  control_regions <- .regression_control_regions(
    moments, cell$post, control_group
  )
  treated_regions <- which(moments$assignment == cell$group)
  period_indices <- match(c(cell$base, cell$post), moments$times)
  rows <- which(
    regression_data$region_index %in% c(treated_regions, control_regions) &
      regression_data$time_index %in% period_indices
  )
  list(rows = rows, control_regions = control_regions,
       treated_regions = treated_regions, period_indices = period_indices)
}

.propensity_block_data <- function(regression_data, moments, cell,
                                   control_group) {
  definition <- .propensity_sample_rows(
    regression_data, moments, cell, control_group
  )
  rows <- definition$rows
  list(
    rows = rows,
    D = as.numeric(regression_data$assignment[rows] == cell$group),
    post = as.numeric(regression_data$time_index[rows] ==
                        definition$period_indices[2L]),
    cluster = regression_data$region_index[rows],
    control_regions = definition$control_regions,
    treated_regions = definition$treated_regions
  )
}

.fit_propensity_master_blocks <- function(moments, att_map, control_group) {
  regression_data <- moments$regression_data
  definitions <- list()
  key_to_id <- new.env(parent = emptyenv())
  cell_master <- integer(nrow(att_map))
  for (k in seq_len(nrow(att_map))) {
    controls <- .regression_control_regions(
      moments, att_map$post[k], control_group
    )
    key <- paste(att_map$group[k], paste(controls, collapse = ","), sep = "|")
    if (!exists(key, key_to_id, inherits = FALSE)) {
      id <- length(definitions) + 1L
      assign(key, id, key_to_id)
      definitions[[id]] <- list(
        master_id = id, key = key, group = att_map$group[k],
        control_regions = controls, cell_ids = integer()
      )
    }
    id <- get(key, key_to_id, inherits = FALSE)
    cell_master[k] <- id
    definitions[[id]]$cell_ids <- c(definitions[[id]]$cell_ids, k)
  }
  blocks <- lapply(definitions, function(definition) {
    treated_regions <- which(moments$assignment == definition$group)
    descendant_periods <- sort(unique(c(
      att_map$base[definition$cell_ids], att_map$post[definition$cell_ids]
    )))
    period_indices <- match(descendant_periods, moments$times)
    rows <- which(
      regression_data$region_index %in%
        c(treated_regions, definition$control_regions) &
        regression_data$time_index %in% period_indices
    )
    X <- regression_data$X[rows, , drop = FALSE]
    D <- as.numeric(regression_data$assignment[rows] == definition$group)
    w <- regression_data$weight[rows]
    fit_start <- proc.time()[["elapsed"]]
    fit <- .fit_weighted_logit(X, D, w)
    elapsed <- proc.time()[["elapsed"]] - fit_start
    c(definition, list(
      rows = rows, periods = descendant_periods,
      coefficients = fit$coefficients, fit = fit,
      elapsed_seconds = elapsed
    ))
  })
  list(blocks = blocks, cell_master = cell_master)
}

.fit_propensity_blocks <- function(moments, att_map, control_group,
                                   ps_fit = c("warm_start", "standard")) {
  ps_fit <- match.arg(ps_fit)
  regression_data <- moments$regression_data
  H <- length(moments$regions)
  p <- ncol(regression_data$X)
  masters <- .fit_propensity_master_blocks(moments, att_map, control_group)
  blocks <- vector("list", nrow(att_map))
  for (k in seq_len(nrow(att_map))) {
    sample <- .propensity_block_data(
      regression_data, moments, att_map[k, ], control_group
    )
    X <- regression_data$X[sample$rows, , drop = FALSE]
    y <- regression_data$y[sample$rows]
    w <- regression_data$weight[sample$rows]
    master <- masters$blocks[[masters$cell_master[k]]]
    start <- rep(0, p)
    start_source <- "zero"
    master_start <- master$coefficients
    one_step_start <- master_start
    master_score_norm <- one_step_score_norm <- NA_real_
    if (all(is.finite(master_start))) {
      master_score_norm <- max(abs(.weighted_logit_state(
        X, sample$D, w, master_start, keep_score_rows = FALSE
      )$score)) / sum(w)
    }
    if (ps_fit == "warm_start" && master$fit$status == "valid") {
      one_step_start <- .one_logit_newton_step(
        X, sample$D, w, master_start
      )
      start <- one_step_start
      start_source <- "master_one_step"
    }
    if (all(is.finite(one_step_start))) {
      one_step_score_norm <- max(abs(.weighted_logit_state(
        X, sample$D, w, one_step_start, keep_score_rows = FALSE
      )$score)) / sum(w)
    }
    fit_start <- proc.time()[["elapsed"]]
    fit <- .fit_weighted_logit(X, sample$D, w, start = start)
    fit_elapsed <- proc.time()[["elapsed"]] - fit_start
    overlap <- .propensity_overlap_status(X, sample$D)
    if (fit$status == "valid" && !overlap$valid) {
      fit$converged <- FALSE
      fit$status <- overlap$status
      fit$reason <- overlap$reason
    }
    state <- if (all(is.finite(fit$coefficients)))
      .weighted_logit_state(
        X, sample$D, w, fit$coefficients, keep_score_rows = FALSE
      ) else
      NULL
    overlap_state <- if (!is.null(overlap$coefficients) &&
                         all(is.finite(overlap$coefficients)))
      .weighted_logit_state(
        X, sample$D, rep(1, nrow(X)), overlap$coefficients,
        keep_score_rows = FALSE
      ) else
      NULL
    cluster_score <- matrix(0, nrow = H, ncol = p)
    cluster_information <- array(0, dim = c(p, p, H))
    cluster_weight <- numeric(H)
    overlap_cluster_score <- matrix(0, nrow = H, ncol = p)
    overlap_cluster_information <- array(0, dim = c(p, p, H))
    overlap_cluster_weight <- numeric(H)
    if (!is.null(state)) {
      for (h in unique(sample$cluster)) {
        local <- which(sample$cluster == h)
        cluster_score[h, ] <- colSums(
          X[local, , drop = FALSE] *
            (w[local] * (sample$D[local] - state$probability[local]))
        )
        Xh <- X[local, , drop = FALSE]
        cluster_information[, , h] <- crossprod(
          Xh, state$information_weight[local] * Xh
        )
        cluster_weight[h] <- sum(w[local])
      }
    }
    if (!is.null(overlap_state)) {
      unit_w <- rep(1, nrow(X))
      for (h in unique(sample$cluster)) {
        local <- which(sample$cluster == h)
        overlap_cluster_score[h, ] <- colSums(
          X[local, , drop = FALSE] *
            (unit_w[local] *
               (sample$D[local] - overlap_state$probability[local]))
        )
        Xh <- X[local, , drop = FALSE]
        overlap_cluster_information[, , h] <- crossprod(
          Xh, overlap_state$information_weight[local] * Xh
        )
        overlap_cluster_weight[h] <- sum(unit_w[local])
      }
    }
    blocks[[k]] <- c(sample, list(
      block_id = k, group = att_map$group[k], base = att_map$base[k],
      post_period = att_map$post[k], master_id = masters$cell_master[k],
      X = X, y = y, w = w, unit_w = rep(1, nrow(X)),
      start_source = start_source, master_start = master_start,
      one_step_start = one_step_start, fit = fit,
      master_score_norm = master_score_norm,
      one_step_score_norm = one_step_score_norm,
      elapsed_seconds = fit_elapsed,
      overlap = overlap,
      coefficients = fit$coefficients,
      probability = if (is.null(state)) fit$probability else state$probability,
      cluster_score = cluster_score,
      cluster_information = cluster_information,
      cluster_weight = cluster_weight,
      overlap_cluster_score = overlap_cluster_score,
      overlap_cluster_information = overlap_cluster_information,
      overlap_cluster_weight = overlap_cluster_weight
    ))
  }
  metadata <- data.frame(
    block_id = seq_along(blocks),
    group = att_map$group,
    base = att_map$base,
    post = att_map$post,
    master_id = vapply(blocks, `[[`, integer(1), "master_id"),
    control_regions = vapply(blocks, function(block) {
      paste(as.character(moments$regions[block$control_regions]), collapse = ",")
    }, character(1)),
    treated_weight = vapply(blocks, function(block) {
      sum(regression_data$weight[block$rows] * block$D)
    }, numeric(1)),
    control_weight = vapply(blocks, function(block) {
      sum(regression_data$weight[block$rows] * (1-block$D))
    }, numeric(1)),
    converged = vapply(blocks, function(block) block$fit$converged, logical(1)),
    iterations = vapply(blocks, function(block) block$fit$iterations, integer(1)),
    score_norm = vapply(blocks, function(block) block$fit$score_norm, numeric(1)),
    master_score_norm = vapply(blocks, `[[`, numeric(1), "master_score_norm"),
    one_step_score_norm = vapply(blocks, `[[`, numeric(1), "one_step_score_norm"),
    elapsed_seconds = vapply(blocks, `[[`, numeric(1), "elapsed_seconds"),
    loglik = vapply(blocks, function(block) block$fit$loglik, numeric(1)),
    rank = vapply(blocks, function(block) block$fit$rank, integer(1)),
    reciprocal_condition = vapply(
      blocks, function(block) block$fit$reciprocal_condition, numeric(1)
    ),
    status = vapply(blocks, function(block) block$fit$status, character(1)),
    reason = vapply(blocks, function(block) block$fit$reason, character(1)),
    start_source = vapply(blocks, `[[`, character(1), "start_source"),
    stringsAsFactors = FALSE
  )
  master_metadata <- data.frame(
    master_id = seq_along(masters$blocks),
    group = vapply(masters$blocks, `[[`, numeric(1), "group"),
    periods = vapply(masters$blocks, function(block) {
      paste(block$periods, collapse = ",")
    }, character(1)),
    control_regions = vapply(masters$blocks, function(block) {
      paste(as.character(moments$regions[block$control_regions]), collapse = ",")
    }, character(1)),
    descendants = vapply(masters$blocks, function(block) {
      paste(block$cell_ids, collapse = ",")
    }, character(1)),
    iterations = vapply(masters$blocks, function(block) block$fit$iterations,
                         integer(1)),
    elapsed_seconds = vapply(masters$blocks, `[[`, numeric(1),
                             "elapsed_seconds"),
    status = vapply(masters$blocks, function(block) block$fit$status,
                    character(1)),
    stringsAsFactors = FALSE
  )
  list(blocks = blocks, metadata = metadata, masters = masters$blocks,
       master_metadata = master_metadata, ps_fit = ps_fit)
}
