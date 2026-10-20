.delete_nuisance_blocks <- function(nuisance, cluster) {
  B <- length(nuisance$blocks)
  deleted <- vector("list", B)
  for (block_id in seq_len(B)) {
    block <- nuisance$blocks[[block_id]]
    changed <- block$cluster_weight[cluster] > 0
    if (!changed) {
      deleted[[block_id]] <- list(
        coefficients = block$coefficients,
        rank = block$rank,
        condition = block$condition,
        status = block$status,
        reason = block$reason,
        changed = FALSE,
        rank_changed = FALSE,
        max_coefficient_change = 0
      )
      next
    }
    solved <- .solve_regression_crossproducts(
      block$Q - block$cluster_Q[, , cluster],
      block$q - block$cluster_q[cluster, ],
      block$total_weight - block$cluster_weight[cluster]
    )
    rank_changed <- solved$rank != block$rank
    if (block$status == "valid" && rank_changed) {
      solved$status <- "rank_changed"
      solved$reason <- paste0(
        "deleting the cluster changes block ", block_id,
        " rank from ", block$rank, " to ", solved$rank
      )
      solved$coefficients[] <- NA_real_
    }
    solved$changed <- TRUE
    solved$rank_changed <- rank_changed
    solved$max_coefficient_change <- if (
      solved$status == "valid" && block$status == "valid"
    ) {
      max(abs(solved$coefficients - block$coefficients))
    } else {
      NA_real_
    }
    deleted[[block_id]] <- solved
  }
  deleted
}

.delete_regression_att_matrix <- function(moments, att_map, evaluation) {
  nuisance <- evaluation$nuisance
  treated_moments <- evaluation$treated_moments
  full_att <- evaluation$att
  H <- length(moments$regions)
  K <- nrow(att_map)
  J <- matrix(rep(full_att, each = H), nrow = H, ncol = K)
  status <- matrix(0L, H, K)
  status_label <- matrix("unaffected", H, K)
  reason <- matrix("unaffected", H, K)
  rank_changed <- matrix(FALSE, H, K)
  failed_block <- matrix("", H, K)
  max_coefficient_change <- matrix(0, H, K)
  influence_source <- matrix("unaffected", H, K)

  B <- length(nuisance$blocks)
  block_status <- matrix("unaffected", H, B)
  block_reason <- matrix("", H, B)
  block_rank <- matrix(NA_integer_, H, B)
  block_max_coefficient_change <- matrix(0, H, B)

  row_labels <- as.character(moments$regions)
  att_labels <- .att_names(att_map)
  dimnames(J) <- dimnames(status) <- dimnames(status_label) <- dimnames(reason) <-
    dimnames(rank_changed) <- dimnames(failed_block) <-
    dimnames(max_coefficient_change) <- dimnames(influence_source) <-
    list(row_labels, att_labels)
  dimnames(block_status) <- dimnames(block_reason) <- dimnames(block_rank) <-
    dimnames(block_max_coefficient_change) <-
    list(row_labels, paste0("block:", seq_len(B)))

  unsupported_full <- !is.finite(full_att)
  if (any(unsupported_full)) {
    status[, unsupported_full] <- 2L
    status_label[, unsupported_full] <- matrix(
      evaluation$support$status[unsupported_full], H, sum(unsupported_full),
      byrow = TRUE
    )
    reason[, unsupported_full] <- matrix(
      paste0("full-sample ATT is unsupported: ",
             evaluation$support$reason[unsupported_full]),
      H, sum(unsupported_full), byrow = TRUE
    )
    influence_source[, unsupported_full] <- "unsupported_full_sample"
  }

  for (h in seq_len(H)) {
    deleted_blocks <- .delete_nuisance_blocks(nuisance, h)
    for (block_id in seq_len(B)) {
      block_status[h, block_id] <- if (deleted_blocks[[block_id]]$changed)
        deleted_blocks[[block_id]]$status else "unaffected"
      block_reason[h, block_id] <- deleted_blocks[[block_id]]$reason
      block_rank[h, block_id] <- deleted_blocks[[block_id]]$rank
      block_max_coefficient_change[h, block_id] <-
        deleted_blocks[[block_id]]$max_coefficient_change
    }

    affected_blocks <- which(vapply(
      nuisance$blocks, function(block) block$cluster_weight[h] > 0, logical(1)
    ))
    regression_cells <- which(
      nuisance$cell_blocks[, "base_block"] %in% affected_blocks |
        nuisance$cell_blocks[, "post_block"] %in% affected_blocks
    )
    treated_cells <- which(vapply(treated_moments, function(treated) {
      treated$base_weight_by_cluster[h] > 0 ||
        treated$post_weight_by_cluster[h] > 0
    }, logical(1)))
    affected_cells <- union(regression_cells, treated_cells)
    for (k in affected_cells) {
      treated <- treated_moments[[k]]
      base_block_id <- nuisance$cell_blocks[k, "base_block"]
      post_block_id <- nuisance$cell_blocks[k, "post_block"]
      base_deleted <- deleted_blocks[[base_block_id]]
      post_deleted <- deleted_blocks[[post_block_id]]
      treated_changed <- treated$base_weight_by_cluster[h] > 0 ||
        treated$post_weight_by_cluster[h] > 0
      regression_changed <- base_deleted$changed || post_deleted$changed
      if (!is.finite(full_att[k])) next
      if (!treated_changed && !regression_changed) next

      status[h, k] <- 1L
      status_label[h, k] <- "updated"
      reason[h, k] <- "updated"
      influence_source[h, k] <- if (treated_changed && regression_changed) {
        "treated_moments_and_control_regression"
      } else if (treated_changed) {
        "treated_moments"
      } else {
        "control_regression"
      }

      base_weight <- treated$base_weight - treated$base_weight_by_cluster[h]
      post_weight <- treated$post_weight - treated$post_weight_by_cluster[h]
      base_outcome <- treated$base_outcome - treated$base_outcome_by_cluster[h]
      post_outcome <- treated$post_outcome - treated$post_outcome_by_cluster[h]
      pooled_x <- treated$pooled_x - treated$pooled_x_by_cluster[h, ]
      if (base_weight <= 0 || post_weight <= 0) {
        J[h, k] <- NA_real_
        status[h, k] <- 2L
        cohort_mass <- sum(moments$region_total_weight[
          moments$assignment == att_map$group[k]
        ]) - if (moments$assignment[h] == att_map$group[k])
          moments$region_total_weight[h] else 0
        status_label[h, k] <- if (cohort_mass <= 0) "cohort_removed" else
          "zero_weight_treated_cell"
        reason[h, k] <- "deletion removes positive treated support"
        next
      }
      block_pair <- list(base_deleted, post_deleted)
      block_ids <- c(base_block_id, post_block_id)
      invalid <- vapply(block_pair, `[[`, character(1), "status") != "valid"
      if (any(invalid)) {
        J[h, k] <- NA_real_
        status[h, k] <- 2L
        invalid_ids <- unique(block_ids[invalid])
        invalid_blocks <- block_pair[invalid]
        failed_block[h, k] <- paste0("block:", invalid_ids, collapse = ",")
        rank_changed[h, k] <- any(vapply(invalid_blocks, `[[`, logical(1),
                                          "rank_changed"))
        status_label[h, k] <- if (rank_changed[h, k]) "rank_changed" else
          invalid_blocks[[1L]]$status
        reason[h, k] <- paste(unique(vapply(invalid_blocks, `[[`, character(1),
                                             "reason")), collapse = "; ")
        max_coefficient_change[h, k] <- NA_real_
        next
      }
      changes <- c(base_deleted$max_coefficient_change,
                   post_deleted$max_coefficient_change)
      max_coefficient_change[h, k] <- max(changes, na.rm = TRUE)
      pooled_weight <- base_weight + post_weight
      J[h, k] <- post_outcome / post_weight - base_outcome / base_weight -
        as.numeric(crossprod(
          pooled_x / pooled_weight,
          post_deleted$coefficients - base_deleted$coefficients
        ))
    }
  }

  list(
    estimates = J,
    changes = sweep(J, 2L, full_att, `-`),
    status = status,
    status_label = status_label,
    reason = reason,
    status_legend = c(unchanged = 0L, updated = 1L, unsupported = 2L),
    nuisance_diagnostics = list(
      block_status = block_status,
      block_reason = block_reason,
      block_rank = block_rank,
      block_max_coefficient_change = block_max_coefficient_change,
      att_rank_changed = rank_changed,
      att_failed_block = failed_block,
      att_max_coefficient_change = max_coefficient_change,
      att_influence_source = influence_source
    )
  )
}
