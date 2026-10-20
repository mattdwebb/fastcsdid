make_test_data <- function(unbalanced = FALSE) {
  regions <- 1:8
  times <- 1:6
  cohort <- c(3, 3, 5, 5, 0, 0, 0, 0)
  out <- vector("list", 0L)
  z <- 0L
  for (h in regions) {
    for (t in times) {
      n <- if (unbalanced) 3L + ((h + 2L * t) %% 5L) else 5L
      for (j in seq_len(n)) {
        z <- z + 1L
        treated <- cohort[h] > 0 && t >= cohort[h]
        out[[z]] <- data.frame(
          id = z,
          region = h,
          time = t,
          cohort = cohort[h],
          y = 0.4 * h + 0.2 * t + 1.5 * treated + 0.03 * j
        )
      }
    }
  }
  do.call(rbind, out)
}

manual_att <- function(d, g, base, post, control = c("never", "not_yet")) {
  control <- match.arg(control)
  treated <- d$cohort == g
  comparison <- if (control == "never") d$cohort == 0 else d$cohort == 0 | d$cohort > post
  mean(d$y[treated & d$time == post]) - mean(d$y[treated & d$time == base]) -
    mean(d$y[comparison & d$time == post]) + mean(d$y[comparison & d$time == base])
}

numeric_labels <- function(x) {
  if (is.factor(x)) x <- as.character(x)
  suppressWarnings(as.numeric(x))
}

reference_att <- function(d, att_map, control = c("never", "not_yet"),
                          weights = NULL, never_treated = 0) {
  control <- match.arg(control)
  y <- numeric_labels(d$y)
  time <- numeric_labels(d$time)
  cohort <- numeric_labels(d$cohort)
  never <- numeric_labels(never_treated)
  cohort[cohort == never] <- Inf
  w <- if (is.null(weights)) rep(1, nrow(d)) else numeric_labels(d[[weights]])
  weighted_cell_mean <- function(selector) {
    denominator <- sum(w[selector])
    if (!is.finite(denominator) || denominator <= 0) return(NA_real_)
    sum(w[selector] * y[selector]) / denominator
  }

  out <- vapply(seq_len(nrow(att_map)), function(k) {
    g <- att_map$group[k]
    b <- att_map$base[k]
    t <- att_map$post[k]
    treated <- cohort == g
    comparison <- if (control == "never") {
      is.infinite(cohort)
    } else {
      is.infinite(cohort) | cohort > t
    }
    cells <- c(
      weighted_cell_mean(treated & time == b),
      weighted_cell_mean(treated & time == t),
      weighted_cell_mean(comparison & time == b),
      weighted_cell_mean(comparison & time == t)
    )
    if (anyNA(cells)) return(NA_real_)
    (cells[2L] - cells[1L]) - (cells[4L] - cells[3L])
  }, numeric(1))
  names(out) <- paste0("ATT(", att_map$group, ",", att_map$post, ")")
  out
}

reference_group_size <- function(d, groups, weights = NULL, never_treated = 0) {
  cohort <- numeric_labels(d$cohort)
  never <- numeric_labels(never_treated)
  cohort[cohort == never] <- Inf
  w <- if (is.null(weights)) rep(1, nrow(d)) else numeric_labels(d[[weights]])
  stats::setNames(vapply(groups, function(g) sum(w[cohort == g]), numeric(1)),
                  as.character(groups))
}

reference_aggregate <- function(att, att_map, group_size,
                                types = c("simple", "group", "calendar", "dynamic")) {
  rows <- list()
  add <- function(label, selector, raw_weight) {
    if (!any(selector) || anyNA(att[selector]) || sum(raw_weight) <= 0) {
      value <- NA_real_
    } else {
      value <- sum(att[selector] * raw_weight / sum(raw_weight))
    }
    rows[[label]] <<- value
  }
  gs <- as.numeric(group_size[as.character(att_map$group)])
  if ("simple" %in% types) add("simple", rep(TRUE, length(att)), gs)
  if ("group" %in% types) {
    for (g in sort(unique(att_map$group))) {
      selector <- att_map$group == g
      add(paste0("group:", g), selector, rep(1, sum(selector)))
    }
  }
  if ("calendar" %in% types) {
    for (tt in sort(unique(att_map$post))) {
      selector <- att_map$post == tt
      add(paste0("calendar:", tt), selector, gs[selector])
    }
  }
  if ("dynamic" %in% types) {
    for (event in sort(unique(att_map$event))) {
      selector <- att_map$event == event
      add(paste0("event:", event), selector, gs[selector])
    }
  }
  unlist(rows)
}

reference_fit <- function(d, att_map, control = "never", weights = NULL,
                          never_treated = 0) {
  att <- reference_att(d, att_map, control, weights, never_treated)
  sizes <- reference_group_size(d, sort(unique(att_map$group)), weights,
                                never_treated)
  list(
    att = att,
    aggregate = reference_aggregate(att, att_map, sizes),
    group_size = sizes
  )
}

relabel_data <- function(d, regions, assignment, never_treated = 0) {
  out <- d
  raw_assignment <- assignment
  raw_assignment[is.infinite(raw_assignment)] <- numeric_labels(never_treated)
  out$cohort <- raw_assignment[match(as.character(out$region), as.character(regions))]
  out
}
