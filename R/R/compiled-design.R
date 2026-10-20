.compile_csdid_design <- function(moments, att_map, control_group,
                                  aggregation_types) {
  control_group <- match.arg(control_group, c("never", "not_yet"))
  aggregation_types <- .match_aggregation_types(aggregation_types)

  finite_cohorts <- sort(unique(moments$assignment[is.finite(moments$assignment)]))
  cohort_labels <- c(finite_cohorts, Inf)
  group_labels <- sort(unique(att_map$group))

  aggregate_specs <- list()
  add_spec <- function(type, index, label, cells, weighting) {
    aggregate_specs[[length(aggregate_specs) + 1L]] <<- list(
      type = type,
      index = index,
      label = label,
      cells = which(cells),
      weighting = weighting
    )
  }
  if ("simple" %in% aggregation_types) {
    add_spec("simple", NA_real_, "simple", rep(TRUE, nrow(att_map)), "group_size")
  }
  if ("group" %in% aggregation_types) {
    for (g in group_labels) {
      add_spec("group", g, paste0("group:", g), att_map$group == g, "equal")
    }
  }
  if ("calendar" %in% aggregation_types) {
    for (tt in sort(unique(att_map$post))) {
      add_spec("calendar", tt, paste0("calendar:", tt),
               att_map$post == tt, "group_size")
    }
  }
  if ("dynamic" %in% aggregation_types) {
    for (event in sort(unique(att_map$event))) {
      add_spec("dynamic", event, paste0("event:", event),
               att_map$event == event, "group_size")
    }
  }

  control_indices <- lapply(att_map$post, function(tt) {
    if (control_group == "never") {
      which(is.infinite(cohort_labels))
    } else {
      which(is.infinite(cohort_labels) | cohort_labels > tt)
    }
  })

  aggregate_labels <- vapply(aggregate_specs, `[[`, character(1), "label")
  list(
    att_map = att_map,
    control_group = control_group,
    aggregation_types = aggregation_types,
    H = length(moments$regions),
    TT = length(moments$times),
    K = nrow(att_map),
    A = length(aggregate_specs),
    cluster_order = moments$regions,
    times = moments$times,
    cohort_labels = cohort_labels,
    cohort_count_target = vapply(cohort_labels, function(label) {
      sum(moments$assignment == label)
    }, integer(1)),
    group_labels = group_labels,
    att_cohort_index = match(att_map$group, cohort_labels),
    att_group_index = match(att_map$group, group_labels),
    base_index = match(att_map$base, moments$times),
    post_index = match(att_map$post, moments$times),
    control_indices = control_indices,
    aggregate_specs = aggregate_specs,
    att_names = .att_names(att_map),
    aggregate_labels = aggregate_labels,
    statistic_names = c(
      paste0("att:", .att_names(att_map)),
      paste0("aggregate:", aggregate_labels)
    ),
    status_legend = c(
      valid = 0L,
      missing_treated_base = 1L,
      missing_treated_post = 2L,
      missing_control_base = 4L,
      missing_control_post = 8L,
      undefined_aggregation = 16L,
      zero_aggregate_weight = 32L
    )
  )
}

.compiled_aggregate_metadata <- function(compiled, estimates) {
  data.frame(
    type = vapply(compiled$aggregate_specs, `[[`, character(1), "type"),
    index = vapply(compiled$aggregate_specs, function(x) {
      if (is.na(x$index)) NA_real_ else as.numeric(x$index)
    }, numeric(1)),
    label = compiled$aggregate_labels,
    estimate = as.numeric(estimates),
    stringsAsFactors = FALSE
  )
}
