.delete_att_matrix <- function(moments, att_map, full_att,
                               control_group = c("never", "not_yet")) {
  control_group <- match.arg(control_group)
  H <- length(moments$regions)
  K <- nrow(att_map)
  A <- moments$assignment
  J <- matrix(rep(full_att, each = H), nrow = H, ncol = K)
  status <- matrix(0L, nrow = H, ncol = K) # 0 unchanged, 1 updated, 2 unsupported
  reason <- matrix("unaffected", nrow = H, ncol = K)
  status_label <- matrix("unaffected", nrow = H, ncol = K)
  rownames(J) <- rownames(status) <- rownames(reason) <- rownames(status_label) <-
    as.character(moments$regions)
  colnames(J) <- colnames(status) <- colnames(reason) <- colnames(status_label) <-
    .att_names(att_map)

  for (k in seq_len(K)) {
    g <- att_map$group[k]
    ib <- match(att_map$base[k], moments$times)
    it <- match(att_map$post[k], moments$times)
    treated <- A == g
    control <- if (control_group == "never") {
      is.infinite(A)
    } else {
      is.infinite(A) | A > att_map$post[k]
    }
    affected <- treated | control

    ntb <- sum(moments$weight_sum[treated, ib]) -
      as.numeric(treated) * moments$weight_sum[, ib]
    ntp <- sum(moments$weight_sum[treated, it]) -
      as.numeric(treated) * moments$weight_sum[, it]
    ncb <- sum(moments$weight_sum[control, ib]) -
      as.numeric(control) * moments$weight_sum[, ib]
    ncp <- sum(moments$weight_sum[control, it]) -
      as.numeric(control) * moments$weight_sum[, it]

    stb <- sum(moments$outcome_sum[treated, ib]) -
      as.numeric(treated) * moments$outcome_sum[, ib]
    stp <- sum(moments$outcome_sum[treated, it]) -
      as.numeric(treated) * moments$outcome_sum[, it]
    scb <- sum(moments$outcome_sum[control, ib]) -
      as.numeric(control) * moments$outcome_sum[, ib]
    scp <- sum(moments$outcome_sum[control, it]) -
      as.numeric(control) * moments$outcome_sum[, it]

    valid <- ntb > 0 & ntp > 0 & ncb > 0 & ncp > 0
    deleted <- (.safe_ratio(stp, ntp) - .safe_ratio(stb, ntb)) -
      (.safe_ratio(scp, ncp) - .safe_ratio(scb, ncb))
    J[affected & valid, k] <- deleted[affected & valid]
    J[affected & !valid, k] <- NA_real_
    status[affected & valid, k] <- 1L
    status[affected & !valid, k] <- 2L
    reason[affected & valid, k] <- "updated"
    status_label[affected & valid, k] <- "updated"
    invalid_rows <- which(affected & !valid)
    for (h in invalid_rows) {
      missing_cells <- c(
        treated_base = ntb[h], treated_post = ntp[h],
        control_base = ncb[h], control_post = ncp[h]
      ) <= 0
      reason[h, k] <- paste0(
        "zero weight in ", paste(names(missing_cells)[missing_cells], collapse = ", ")
      )
      treated_mass <- sum(moments$region_total_weight[treated]) -
        as.numeric(treated[h]) * moments$region_total_weight[h]
      control_mass <- sum(moments$region_total_weight[control]) -
        as.numeric(control[h]) * moments$region_total_weight[h]
      status_label[h, k] <- if (treated_mass <= 0) {
        "cohort_removed"
      } else if (control_mass <= 0) {
        "control_group_removed"
      } else if (ntb[h] <= 0 || ntp[h] <= 0) {
        "zero_weight_treated_cell"
      } else {
        "zero_weight_control_cell"
      }
    }
  }
  list(
    estimates = J,
    changes = sweep(J, 2L, full_att, `-`),
    status = status,
    status_label = status_label,
    reason = reason,
    status_legend = c(unchanged = 0L, updated = 1L, unsupported = 2L)
  )
}

.delete_aggregate_matrix <- function(moments, att_map, delete_att,
                                     aggregation_types,
                                     support = c("changing", "fixed")) {
  support <- match.arg(support)
  H <- nrow(delete_att)
  full_group_size <- stats::setNames(vapply(sort(unique(att_map$group)), function(g) {
    sum(moments$region_total_weight[moments$assignment == g])
  }, numeric(1)), as.character(sort(unique(att_map$group))))
  # Determine required cells from full-sample masses, not an already-deleted row.
  template <- aggregate_att(rep(0, nrow(att_map)), att_map,
                            full_group_size, aggregation_types)
  labels <- template$metadata$label
  out <- matrix(NA_real_, nrow = H, ncol = length(labels),
                dimnames = list(as.character(moments$regions), labels))
  changed <- matrix(FALSE, nrow = H, ncol = length(labels), dimnames = dimnames(out))

  for (h in seq_len(H)) {
    gs <- full_group_size
    gh <- moments$assignment[h]
    if (is.finite(gh) && as.character(gh) %in% names(gs)) {
      gs[as.character(gh)] <- gs[as.character(gh)] - moments$region_total_weight[h]
    }
    agg <- aggregate_att(delete_att[h, ], att_map, gs, aggregation_types)
    out[h, ] <- agg$estimates[labels]
    removed <- gs[as.character(att_map$group)] <= 0
    for (j in seq_along(labels)) {
      needed <- is.finite(template$weights[j, ]) & template$weights[j, ] > 0
      changed[h, j] <- any(needed & (removed | !is.finite(delete_att[h, ])))
      if (support == "fixed" && changed[h, j]) out[h, j] <- NA_real_
    }
  }
  status <- ifelse(is.na(out), 2L, 1L)
  dimnames(status) <- dimnames(out)
  reason <- ifelse(
    is.na(out),
    "one or more required ATT cells are unsupported, or aggregate weight is zero",
    "updated"
  )
  dimnames(reason) <- dimnames(out)
  reason[changed] <- "deletion changes full-sample aggregate support"
  status_label <- ifelse(is.na(out), "undefined_aggregation", "updated")
  dimnames(status_label) <- dimnames(out)
  list(
    estimates = out,
    status = status,
    status_label = status_label,
    reason = reason,
    support = support,
    support_changed = changed,
    status_legend = c(updated = 1L, unsupported = 2L)
  )
}

.jackknife_psd_check <- function(covariance, keep = NULL,
                                 tolerance = 1e-10) {
  if (is.null(keep)) keep <- rep(TRUE, nrow(covariance))
  block <- covariance[keep, keep, drop = FALSE]
  if (!nrow(block)) return(NA)
  if (any(!is.finite(block))) return(NA)
  eigenvalues <- eigen((block + t(block)) / 2, symmetric = TRUE,
                       only.values = TRUE)$values
  scale <- max(1, max(abs(eigenvalues)))
  min(eigenvalues) >= -tolerance * scale
}

.jackknife_covariance_block <- function(deleted, full) {
  G <- nrow(deleted)
  delta_full <- sweep(deleted, 2L, full, `-`)
  delta_mean <- sweep(deleted, 2L, colMeans(deleted), `-`)
  list(
    cv3 = (G - 1) / G * crossprod(delta_full),
    cv3j = (G - 1) / G * crossprod(delta_mean)
  )
}

#' Compute jackknife covariance matrices
#'
#' @param deleted Matrix whose rows are delete-one-cluster estimates.
#' @param full Full-sample estimate vector.
#' @param strict If `TRUE`, any unsupported deletion makes the corresponding
#'   covariance matrix unavailable. If `FALSE`, complete rows are used.
#' @return Conventional fully-supported, common-support, and explicitly
#'   pairwise CV3 and CV3J covariance results with deletion-count diagnostics.
#' @export
jackknife_covariance <- function(deleted, full, strict = TRUE) {
  deleted <- as.matrix(deleted)
  full <- as.numeric(full)
  if (ncol(deleted) != length(full)) {
    stop("deleted and full have incompatible dimensions.", call. = FALSE)
  }
  H <- nrow(deleted)
  p <- ncol(deleted)
  valid_deletions <- is.finite(deleted) &
    matrix(is.finite(full), H, p, byrow = TRUE)
  complete <- apply(valid_deletions, 1L, all)
  complete_by_stat <- is.finite(full) & colSums(valid_deletions) == H
  names_full <- colnames(deleted)
  blank <- matrix(NA_real_, p, p, dimnames = list(names_full, names_full))

  if (strict) {
    keep <- complete_by_stat
    D <- deleted[, keep, drop = FALSE]
    full_used <- full[keep]
  } else {
    keep <- rep(TRUE, ncol(deleted))
    D <- deleted[complete, , drop = FALSE]
    full_used <- full
  }
  G <- nrow(D)
  if (G < 2L) {
    stop("At least two valid deletions are required.", call. = FALSE)
  }
  cv3 <- cv3j <- blank
  if (any(keep)) {
    conventional <- .jackknife_covariance_block(D, full_used)
    cv3[keep, keep] <- conventional$cv3
    cv3j[keep, keep] <- conventional$cv3j
  }

  finite_estimands <- is.finite(full)
  common_support_rows <- if (any(finite_estimands)) {
    apply(valid_deletions[, finite_estimands, drop = FALSE], 1L, all)
  } else {
    rep(FALSE, H)
  }
  common_repetitions <- sum(common_support_rows)
  common_cv3 <- common_cv3j <- blank
  if (common_repetitions >= 2L && any(finite_estimands)) {
    common <- .jackknife_covariance_block(
      deleted[common_support_rows, finite_estimands, drop = FALSE],
      full[finite_estimands]
    )
    common_cv3[finite_estimands, finite_estimands] <- common$cv3
    common_cv3j[finite_estimands, finite_estimands] <- common$cv3j
  }

  n_valid_cov <- crossprod(valid_deletions * 1L)
  dimnames(n_valid_cov) <- dimnames(blank)
  pairwise_cv3 <- pairwise_cv3j <- blank
  for (i in seq_len(p)) {
    if (!is.finite(full[i])) next
    for (j in i:p) {
      if (!is.finite(full[j])) next
      pair_rows <- valid_deletions[, i] & valid_deletions[, j]
      pair_G <- sum(pair_rows)
      if (pair_G < 2L) next
      di <- deleted[pair_rows, i]
      dj <- deleted[pair_rows, j]
      factor <- (pair_G - 1) / pair_G
      pairwise_cv3[i, j] <- pairwise_cv3[j, i] <-
        factor * sum((di - full[i]) * (dj - full[j]))
      pairwise_cv3j[i, j] <- pairwise_cv3j[j, i] <-
        factor * sum((di - mean(di)) * (dj - mean(dj)))
    }
  }

  repetitions_by_stat <- if (strict) {
    ifelse(complete_by_stat, H, 0L)
  } else {
    rep.int(G, p)
  }
  list(
    cv3 = cv3,
    cv3j = cv3j,
    se_cv3 = sqrt(diag(cv3)),
    se_cv3j = sqrt(diag(cv3j)),
    repetitions = G,
    repetitions_by_stat = stats::setNames(repetitions_by_stat, names_full),
    complete = all(complete_by_stat),
    complete_by_stat = stats::setNames(complete_by_stat, names_full),
    valid_deletions = `colnames<-`(valid_deletions, names_full),
    n_valid_by_stat = stats::setNames(colSums(valid_deletions), names_full),
    n_valid_cov = n_valid_cov,
    common_support_rows = common_support_rows,
    common_repetitions = common_repetitions,
    common_cv3 = common_cv3,
    common_cv3j = common_cv3j,
    common_se_cv3 = sqrt(diag(common_cv3)),
    common_se_cv3j = sqrt(diag(common_cv3j)),
    pairwise_cv3 = pairwise_cv3,
    pairwise_cv3j = pairwise_cv3j,
    pairwise_se_cv3 = sqrt(diag(pairwise_cv3)),
    pairwise_se_cv3j = sqrt(diag(pairwise_cv3j)),
    covariance_type = c(
      cv3 = "all deletions for fully supported estimands",
      common_cv3 = "one common deletion set across finite estimands",
      pairwise_cv3 = "entry-specific pairwise deletion sets"
    ),
    positive_semidefinite = c(
      cv3 = .jackknife_psd_check(cv3, keep),
      cv3j = .jackknife_psd_check(cv3j, keep),
      common_cv3 = .jackknife_psd_check(common_cv3, finite_estimands),
      common_cv3j = .jackknife_psd_check(common_cv3j, finite_estimands),
      pairwise_cv3 = .jackknife_psd_check(pairwise_cv3),
      pairwise_cv3j = .jackknife_psd_check(pairwise_cv3j)
    )
  )
}
