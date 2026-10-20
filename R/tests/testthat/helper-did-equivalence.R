did_fixture <- function(unbalanced = FALSE, weighted = FALSE,
                        missing_cells = FALSE, first_period_treated = FALSE) {
  regions <- 1:10
  times <- c(2001, 2003, 2006, 2010, 2015)
  cohort_by_region <- c(2006, 2006, 2010, 2010, 2015, 2015, 0, 0, 0, 0)
  if (first_period_treated) cohort_by_region[1L] <- times[1L]
  rows <- vector("list", length(regions) * length(times))
  z <- 0L
  id <- 0L
  for (h in regions) {
    for (tt in times) {
      z <- z + 1L
      n <- if (unbalanced) 2L + ((3L * h + match(tt, times)) %% 5L) else 4L
      obs <- seq_len(n)
      ids <- id + obs
      id <- id + n
      cohort <- cohort_by_region[h]
      treated <- cohort > 0 && tt >= cohort
      rows[[z]] <- data.frame(
        id = ids,
        region = h,
        time = tt,
        cohort = cohort,
        y = h / 7 + match(tt, times) / 5 + obs / 31 +
          treated * (0.6 + h / 20 + match(tt, times) / 10),
        w = if (weighted) 0.15 + ((ids * 7L + h) %% 19L) / 6 else 1
      )
    }
  }
  out <- do.call(rbind, rows)
  if (missing_cells) {
    out <- out[!(out$region == 8L & out$time == 2003) &
      !(out$region == 4L & out$time == 2015 & out$id %% 2L == 0L), , drop = FALSE]
  }
  rownames(out) <- NULL
  out
}

fit_did_reference <- function(data, control_group, weighted = FALSE) {
  did::att_gt(
    yname = "y", tname = "time", idname = "id", gname = "cohort",
    xformla = ~1, data = data, panel = FALSE,
    control_group = if (control_group == "never") "nevertreated" else "notyettreated",
    anticipation = 0, weightsname = if (weighted) "w" else NULL,
    bstrap = FALSE, cband = FALSE, est_method = "dr",
    base_period = "universal", faster_mode = TRUE, compute_inffunc = TRUE
  )
}

compare_did_post_att <- function(fast, reference, tolerance = 1e-10) {
  external <- data.frame(group = reference$group, post = reference$t,
                         did_att = reference$att)
  external <- external[external$post >= external$group, , drop = FALSE]
  internal <- data.frame(group = fast$att_map$group, post = fast$att_map$post,
                         fast_att = unname(fast$att_gt))
  joined <- merge(internal, external, by = c("group", "post"), all = TRUE,
                  sort = TRUE)
  testthat::expect_equal(nrow(joined), nrow(internal))
  testthat::expect_false(anyNA(joined[c("fast_att", "did_att")]))
  testthat::expect_equal(joined$fast_att, joined$did_att, tolerance = tolerance)
  invisible(joined)
}

compare_did_aggregates <- function(fast, reference, tolerance = 1e-10) {
  simple <- did::aggte(reference, type = "simple", bstrap = FALSE, cband = FALSE)
  testthat::expect_equal(unname(fast$aggregation$estimates["simple"]),
                         unname(simple$overall.att), tolerance = tolerance)

  specifications <- list(
    group = c(prefix = "group:", fast_type = "group"),
    calendar = c(prefix = "calendar:", fast_type = "calendar"),
    dynamic = c(prefix = "event:", fast_type = "dynamic")
  )
  for (type in names(specifications)) {
    external <- did::aggte(reference, type = type, bstrap = FALSE, cband = FALSE)
    keep <- if (type == "dynamic") external$egt >= 0 else rep(TRUE, length(external$egt))
    external_estimate <- stats::setNames(
      external$att.egt[keep],
      paste0(specifications[[type]]["prefix"], external$egt[keep])
    )
    fast_rows <- fast$aggregation$metadata$type == specifications[[type]]["fast_type"]
    internal_estimate <- fast$aggregation$estimates[fast_rows]
    testthat::expect_setequal(names(internal_estimate), names(external_estimate))
    testthat::expect_equal(
      unname(internal_estimate[sort(names(internal_estimate))]),
      unname(external_estimate[sort(names(external_estimate))]),
      tolerance = tolerance,
      info = paste(type, "aggregation")
    )
  }
  invisible(TRUE)
}
