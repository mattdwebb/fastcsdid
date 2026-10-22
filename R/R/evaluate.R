#' Evaluate all ATT(g,t) effects for a region assignment
#'
#' @param moments Collapsed region-period moments.
#' @param assignment Region-level treatment cohorts, internally using `Inf` for
#'   never treated. Defaults to the observed assignment.
#' @param att_map Compiled ATT map.
#' @param control_group Either `"never"` or `"not_yet"`.
#' @param active Optional logical vector indicating regions retained in the
#'   sample. Set an entry to `FALSE` to exclude that region.
#' @return A list containing the ATT vector, support diagnostics, and cohort sizes.
#' @export
evaluate_assignment <- function(moments, assignment = moments$assignment,
                                att_map = compile_att_map(moments),
                                control_group = c("never", "not_yet"),
                                active = rep.int(TRUE, length(moments$regions))) {
  control_group <- match.arg(control_group)
  H <- length(moments$regions)
  if (length(assignment) != H) {
    stop("assignment must have one value per region.", call. = FALSE)
  }
  assignment <- .as_numeric_design(assignment, "assignment")
  if (anyNA(assignment)) {
    stop("assignment must not contain missing or nonnumeric values.", call. = FALSE)
  }
  if (length(active) != H || anyNA(active)) {
    stop("active must be a nonmissing logical vector with one value per region.", call. = FALSE)
  }
  active <- as.logical(active)

  K <- nrow(att_map)
  att <- rep(NA_real_, K)
  support <- data.frame(
    treated_base_weight = rep(0, K),
    treated_post_weight = rep(0, K),
    control_base_weight = rep(0, K),
    control_post_weight = rep(0, K),
    supported = rep(FALSE, K),
    status = rep("unsupported", K),
    reason = rep("", K),
    stringsAsFactors = FALSE
  )

  for (k in seq_len(K)) {
    g <- att_map$group[k]
    b <- match(att_map$base[k], moments$times)
    t <- match(att_map$post[k], moments$times)
    treated <- active & assignment == g
    control <- if (control_group == "never") {
      active & is.infinite(assignment)
    } else {
      active & (is.infinite(assignment) | assignment > att_map$post[k])
    }

    ntb <- sum(moments$weight_sum[treated, b])
    ntp <- sum(moments$weight_sum[treated, t])
    ncb <- sum(moments$weight_sum[control, b])
    ncp <- sum(moments$weight_sum[control, t])
    support[k, 1:4] <- c(ntb, ntp, ncb, ncp)

    if (all(c(ntb, ntp, ncb, ncp) > 0)) {
      mtb <- sum(moments$outcome_sum[treated, b]) / ntb
      mtp <- sum(moments$outcome_sum[treated, t]) / ntp
      mcb <- sum(moments$outcome_sum[control, b]) / ncb
      mcp <- sum(moments$outcome_sum[control, t]) / ncp
      att[k] <- (mtp - mtb) - (mcp - mcb)
      support$supported[k] <- TRUE
      support$status[k] <- "valid"
    } else {
      missing_cells <- c(
        treated_base = ntb, treated_post = ntp,
        control_base = ncb, control_post = ncp
      ) <= 0
      support$status[k] <- if (ntb <= 0 || ntp <= 0) {
        "zero_weight_treated_cell"
      } else {
        "zero_weight_control_cell"
      }
      support$reason[k] <- paste0(
        "zero weight in ", paste(names(missing_cells)[missing_cells], collapse = ", ")
      )
    }
  }

  names(att) <- .att_names(att_map)
  finite_groups <- sort(unique(att_map$group))
  group_size <- stats::setNames(vapply(finite_groups, function(g) {
    sum(moments$region_total_weight[active & assignment == g])
  }, numeric(1)), as.character(finite_groups))

  list(att = att, support = support, group_size = group_size, active = active)
}
