#' Aggregate group-time effects
#'
#' Computes simple, cohort, calendar-time, and event-time summaries. Cell weights
#' for simple, calendar, and dynamic summaries are proportional to the treated
#' cohort's total sample weight. Cohort-specific summaries average the cohort's
#' post-treatment effects equally.
#'
#' @param att Numeric ATT vector.
#' @param att_map ATT specification table.
#' @param group_size Named treated-cohort sample sizes.
#' @param types Requested aggregation types.
#' @return A list with estimates, a sparse-in-spirit weight matrix, and metadata.
#' @export
aggregate_att <- function(att, att_map, group_size,
                          types = c("simple", "group", "calendar", "dynamic")) {
  types <- .match_aggregation_types(types)
  K <- nrow(att_map)
  if (length(att) != K) {
    stop("att and att_map have incompatible lengths.", call. = FALSE)
  }
  gs <- as.numeric(group_size[as.character(att_map$group)])
  if (anyNA(gs)) {
    stop("group_size is missing one or more treated cohorts.", call. = FALSE)
  }

  rows <- list()
  add_row <- function(type, index, label, selector, raw_weight) {
    w <- rep(0, K)
    w[selector] <- raw_weight
    active <- selector & is.finite(w) & w > 0
    if (!any(active) || sum(w[active]) <= 0 || anyNA(att[active])) {
      estimate <- NA_real_
      w[] <- NA_real_
    } else {
      w[!active] <- 0
      w <- w / sum(w)
      estimate <- sum(w[active] * att[active])
    }
    rows[[length(rows) + 1L]] <<- list(
      type = type, index = index, label = label, estimate = estimate, weight = w
    )
  }

  if ("simple" %in% types) {
    add_row("simple", NA_real_, "simple", rep(TRUE, K), gs)
  }
  if ("group" %in% types) {
    for (g in sort(unique(att_map$group))) {
      sel <- att_map$group == g
      add_row("group", g, paste0("group:", g), sel, rep(1, sum(sel)))
    }
  }
  if ("calendar" %in% types) {
    for (t in sort(unique(att_map$post))) {
      sel <- att_map$post == t
      add_row("calendar", t, paste0("calendar:", t), sel, gs[sel])
    }
  }
  if ("dynamic" %in% types) {
    for (e in sort(unique(att_map$event))) {
      sel <- att_map$event == e
      add_row("dynamic", e, paste0("event:", e), sel, gs[sel])
    }
  }

  metadata <- data.frame(
    type = vapply(rows, `[[`, character(1), "type"),
    index = vapply(rows, function(x) if (is.na(x$index)) NA_real_ else as.numeric(x$index), numeric(1)),
    label = vapply(rows, `[[`, character(1), "label"),
    estimate = vapply(rows, `[[`, numeric(1), "estimate"),
    stringsAsFactors = FALSE
  )
  weight_matrix <- do.call(rbind, lapply(rows, `[[`, "weight"))
  rownames(weight_matrix) <- metadata$label
  colnames(weight_matrix) <- .att_names(att_map)
  estimates <- metadata$estimate
  names(estimates) <- metadata$label
  list(estimates = estimates, weights = weight_matrix, metadata = metadata)
}
