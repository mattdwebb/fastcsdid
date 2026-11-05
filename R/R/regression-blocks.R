.solve_regression_crossproducts <- function(Q, q, total_weight,
                                            condition_tolerance = .Machine$double.eps) {
  p <- nrow(Q)
  blank <- list(
    coefficients = rep(NA_real_, p), rank = 0L, condition = 0,
    status = "singular_outcome_regression", reason = ""
  )
  if (!is.finite(total_weight) || total_weight <= 0) {
    blank$status <- "zero_weight_regression_sample"
    blank$reason <- "the regression sample has zero positive sampling-weight mass"
    return(blank)
  }
  if (any(!is.finite(Q)) || any(!is.finite(q))) {
    blank$reason <- "the regression cross-products are non-finite"
    return(blank)
  }
  singular_values <- tryCatch(svd(Q, nu = 0L, nv = 0L)$d,
                              error = function(e) numeric())
  if (!length(singular_values)) {
    blank$reason <- "the regression cross-products could not be decomposed"
    return(blank)
  }
  scale <- max(singular_values)
  rank <- if (scale > 0) sum(singular_values > scale * condition_tolerance) else 0L
  condition <- tryCatch(rcond(Q), error = function(e) 0)
  blank$rank <- as.integer(rank)
  blank$condition <- condition
  if (rank < p || !is.finite(condition) || condition < condition_tolerance) {
    blank$reason <- paste0(
      "outcome-regression design rank ", rank, " of ", p,
      "; reciprocal condition number ", format(condition, digits = 5)
    )
    return(blank)
  }
  coefficients <- tryCatch({
    R <- chol(Q)
    backsolve(R, forwardsolve(t(R), q))
  }, error = function(e) NULL)
  if (is.null(coefficients) || any(!is.finite(coefficients))) {
    blank$reason <- "the full-rank regression system could not be solved"
    return(blank)
  }
  list(
    coefficients = as.numeric(coefficients), rank = as.integer(rank),
    condition = condition, status = "valid", reason = ""
  )
}

.regression_control_regions <- function(moments, post, control_group) {
  if (control_group == "never") {
    which(is.infinite(moments$assignment))
  } else {
    which(is.infinite(moments$assignment) | moments$assignment > post)
  }
}

.fit_one_regression_block <- function(regression_data, rows, block_id,
                                      period, control_regions, H) {
  X <- regression_data$X[rows, , drop = FALSE]
  y <- regression_data$y[rows]
  w <- regression_data$weight[rows]
  p <- ncol(X)
  crossproduct_start <- proc.time()[["elapsed"]]
  Q <- crossprod(X, w * X)
  q <- as.numeric(crossprod(X, w * y))
  total_weight <- sum(w)
  crossproduct_time <- proc.time()[["elapsed"]] - crossproduct_start
  solve_start <- proc.time()[["elapsed"]]
  solved <- .solve_regression_crossproducts(Q, q, total_weight)
  solve_time <- proc.time()[["elapsed"]] - solve_start

  cluster_crossproduct_start <- proc.time()[["elapsed"]]
  cluster_Q <- array(0, dim = c(p, p, H))
  cluster_q <- matrix(0, nrow = H, ncol = p)
  cluster_weight <- numeric(H)
  if (length(rows)) {
    row_clusters <- regression_data$region_index[rows]
    for (h in unique(row_clusters)) {
      local <- which(row_clusters == h)
      Xh <- X[local, , drop = FALSE]
      wh <- w[local]
      yh <- y[local]
      cluster_Q[, , h] <- crossprod(Xh, wh * Xh)
      cluster_q[h, ] <- as.numeric(crossprod(Xh, wh * yh))
      cluster_weight[h] <- sum(wh)
    }
  }
  cluster_crossproduct_time <-
    proc.time()[["elapsed"]] - cluster_crossproduct_start
  c(list(
    block_id = block_id,
    period = period,
    control_regions = control_regions,
    rows = rows,
    total_weight = total_weight,
    Q = Q,
    q = q,
    cluster_Q = cluster_Q,
    cluster_q = cluster_q,
    cluster_weight = cluster_weight,
    timing = c(
      crossproducts = crossproduct_time + cluster_crossproduct_time,
      solve = solve_time
    )
  ), solved)
}

.fit_nuisance_blocks <- function(moments, att_map,
                                 control_group = c("never", "not_yet")) {
  control_group <- match.arg(control_group)
  regression_data <- moments$regression_data
  if (is.null(regression_data)) {
    stop("Regression moments are absent; call collapse_rcs() with xformla.",
         call. = FALSE)
  }
  K <- nrow(att_map)
  H <- length(moments$regions)
  cell_blocks <- matrix(NA_integer_, nrow = K, ncol = 2L,
                        dimnames = list(.att_names(att_map),
                                        c("base_block", "post_block")))
  definitions <- list()
  key_to_id <- new.env(parent = emptyenv())
  compilation_start <- proc.time()[["elapsed"]]

  for (k in seq_len(K)) {
    control_regions <- .regression_control_regions(
      moments, att_map$post[k], control_group
    )
    for (role in c("base", "post")) {
      period <- att_map[[role]][k]
      period_index <- match(period, moments$times)
      key <- paste0(period_index, "|", paste(control_regions, collapse = ","))
      if (!.key_cache_has(key_to_id, key)) {
        block_id <- length(definitions) + 1L
        .key_cache_set(key_to_id, key, block_id)
        definitions[[block_id]] <- list(
          block_id = block_id,
          key = key,
          period = period,
          period_index = period_index,
          control_regions = control_regions,
          cell_ids = integer(),
          roles = character()
        )
      }
      block_id <- .key_cache_get(key_to_id, key)
      cell_blocks[k, paste0(role, "_block")] <- block_id
      definitions[[block_id]]$cell_ids <-
        c(definitions[[block_id]]$cell_ids, k)
      definitions[[block_id]]$roles <-
        c(definitions[[block_id]]$roles, role)
    }
  }
  compilation_time <- proc.time()[["elapsed"]] - compilation_start

  blocks <- lapply(definitions, function(definition) {
    rows <- which(
      regression_data$time_index == definition$period_index &
        regression_data$region_index %in% definition$control_regions
    )
    block <- .fit_one_regression_block(
      regression_data, rows, definition$block_id, definition$period,
      definition$control_regions, H
    )
    block$cell_ids <- definition$cell_ids
    block$roles <- definition$roles
    block$key <- definition$key
    block
  })

  metadata <- data.frame(
    block_id = seq_along(blocks),
    period = vapply(blocks, `[[`, numeric(1), "period"),
    control_regions = vapply(blocks, function(block) {
      paste(as.character(moments$regions[block$control_regions]), collapse = ",")
    }, character(1)),
    att_ids = vapply(blocks, function(block) {
      paste(sort(unique(block$cell_ids)), collapse = ",")
    }, character(1)),
    roles = vapply(blocks, function(block) {
      paste(sort(unique(block$roles)), collapse = ",")
    }, character(1)),
    weight = vapply(blocks, `[[`, numeric(1), "total_weight"),
    rank = vapply(blocks, `[[`, integer(1), "rank"),
    columns = ncol(regression_data$X),
    reciprocal_condition = vapply(blocks, `[[`, numeric(1), "condition"),
    status = vapply(blocks, `[[`, character(1), "status"),
    reason = vapply(blocks, `[[`, character(1), "reason"),
    stringsAsFactors = FALSE
  )
  list(
    blocks = blocks, metadata = metadata, cell_blocks = cell_blocks,
    timing = c(
      block_compilation = compilation_time,
      crossproducts = sum(vapply(blocks, function(block) {
        unname(block$timing["crossproducts"])
      }, numeric(1))),
      full_regression_solves = sum(vapply(blocks, function(block) {
        unname(block$timing["solve"])
      }, numeric(1)))
    )
  )
}

.compile_treated_regression_moments <- function(moments, att_map) {
  regression_data <- moments$regression_data
  H <- length(moments$regions)
  p <- ncol(regression_data$X)
  row_keys <- paste(regression_data$assignment, regression_data$time_index,
                    sep = "|")
  row_catalog <- split(seq_along(row_keys), row_keys)
  lookup <- function(group, period) {
    key <- paste(group, match(period, moments$times), sep = "|")
    rows <- row_catalog[[key]]
    if (is.null(rows)) integer() else rows
  }

  lapply(seq_len(nrow(att_map)), function(k) {
    base_rows <- lookup(att_map$group[k], att_map$base[k])
    post_rows <- lookup(att_map$group[k], att_map$post[k])
    all_rows <- c(base_rows, post_rows)
    base_weight_by_cluster <- post_weight_by_cluster <- numeric(H)
    base_outcome_by_cluster <- post_outcome_by_cluster <- numeric(H)
    x_by_cluster <- matrix(0, nrow = H, ncol = p)

    if (length(base_rows)) {
      base_weight_by_cluster <- rowsum(
        regression_data$weight[base_rows], regression_data$region_index[base_rows],
        reorder = FALSE
      )
      base_outcome_by_cluster <- rowsum(
        regression_data$weight[base_rows] * regression_data$y[base_rows],
        regression_data$region_index[base_rows], reorder = FALSE
      )
      bw <- numeric(H); by <- numeric(H)
      bw[as.integer(rownames(base_weight_by_cluster))] <- base_weight_by_cluster[, 1L]
      by[as.integer(rownames(base_outcome_by_cluster))] <- base_outcome_by_cluster[, 1L]
      base_weight_by_cluster <- bw
      base_outcome_by_cluster <- by
    }
    if (length(post_rows)) {
      post_weight_table <- rowsum(
        regression_data$weight[post_rows], regression_data$region_index[post_rows],
        reorder = FALSE
      )
      post_outcome_table <- rowsum(
        regression_data$weight[post_rows] * regression_data$y[post_rows],
        regression_data$region_index[post_rows], reorder = FALSE
      )
      post_weight_by_cluster[as.integer(rownames(post_weight_table))] <-
        post_weight_table[, 1L]
      post_outcome_by_cluster[as.integer(rownames(post_outcome_table))] <-
        post_outcome_table[, 1L]
    }
    if (length(all_rows)) {
      weighted_X <- regression_data$X[all_rows, , drop = FALSE] *
        regression_data$weight[all_rows]
      x_table <- rowsum(weighted_X, regression_data$region_index[all_rows],
                        reorder = FALSE)
      x_by_cluster[as.integer(rownames(x_table)), ] <- x_table
    }
    list(
      base_weight = sum(base_weight_by_cluster),
      post_weight = sum(post_weight_by_cluster),
      base_outcome = sum(base_outcome_by_cluster),
      post_outcome = sum(post_outcome_by_cluster),
      pooled_x = colSums(x_by_cluster),
      base_weight_by_cluster = base_weight_by_cluster,
      post_weight_by_cluster = post_weight_by_cluster,
      base_outcome_by_cluster = base_outcome_by_cluster,
      post_outcome_by_cluster = post_outcome_by_cluster,
      pooled_x_by_cluster = x_by_cluster
    )
  })
}

.evaluate_att_from_blocks <- function(moments, att_map, nuisance) {
  treated_start <- proc.time()[["elapsed"]]
  treated_moments <- .compile_treated_regression_moments(moments, att_map)
  treated_time <- proc.time()[["elapsed"]] - treated_start
  K <- nrow(att_map)
  att <- rep(NA_real_, K)
  support <- data.frame(
    treated_base_weight = numeric(K),
    treated_post_weight = numeric(K),
    control_base_weight = numeric(K),
    control_post_weight = numeric(K),
    supported = FALSE,
    status = "unsupported",
    reason = "",
    base_block = nuisance$cell_blocks[, "base_block"],
    post_block = nuisance$cell_blocks[, "post_block"],
    stringsAsFactors = FALSE
  )

  att_start <- proc.time()[["elapsed"]]
  for (k in seq_len(K)) {
    treated <- treated_moments[[k]]
    base_block <- nuisance$blocks[[nuisance$cell_blocks[k, "base_block"]]]
    post_block <- nuisance$blocks[[nuisance$cell_blocks[k, "post_block"]]]
    support[k, 1:4] <- c(
      treated$base_weight, treated$post_weight,
      base_block$total_weight, post_block$total_weight
    )
    if (treated$base_weight <= 0 || treated$post_weight <= 0) {
      support$status[k] <- "zero_weight_treated_cell"
      support$reason[k] <- "a required treated period has zero sampling-weight mass"
      next
    }
    failed_blocks <- c(base_block$status, post_block$status) != "valid"
    if (any(failed_blocks)) {
      failed <- list(base_block, post_block)[failed_blocks]
      support$status[k] <- failed[[1L]]$status
      support$reason[k] <- paste(unique(vapply(failed, `[[`, character(1), "reason")),
                                 collapse = "; ")
      next
    }
    pooled_weight <- treated$base_weight + treated$post_weight
    pooled_x_mean <- treated$pooled_x / pooled_weight
    att[k] <- treated$post_outcome / treated$post_weight -
      treated$base_outcome / treated$base_weight -
      as.numeric(crossprod(
        pooled_x_mean, post_block$coefficients - base_block$coefficients
      ))
    support$supported[k] <- TRUE
    support$status[k] <- "valid"
  }
  names(att) <- .att_names(att_map)
  groups <- sort(unique(att_map$group))
  group_size <- stats::setNames(vapply(groups, function(group) {
    sum(moments$region_total_weight[moments$assignment == group])
  }, numeric(1)), as.character(groups))
  list(
    att = att, support = support, group_size = group_size,
    treated_moments = treated_moments, nuisance = nuisance,
    timing = c(
      treated_moment_compilation = treated_time,
      att_calculation = proc.time()[["elapsed"]] - att_start
    )
  )
}
