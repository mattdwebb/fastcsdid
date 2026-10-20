test_that("cached exact nonlinear deletion matches old subset refit oracle", {
  data <- make_nonlinear_fixture(H = 12L, observations_per_cell = 6L)
  moments <- collapse_rcs(
    data, "y", "time", "cohort", "region", weights = "w",
    xformla = ~x1 + x2 + female
  )
  att_map <- compile_att_map(moments)
  propensity <- .fit_propensity_blocks(moments, att_map, "not_yet",
                                       "warm_start")

  for (block in propensity$blocks) {
    for (cluster in unique(block$cluster)) {
      if (block$cluster_weight[cluster] <= 0) next
      update <- .deleted_propensity_update(block, cluster)
      if (update$status != "valid") next
      keep <- block$cluster != cluster

      old_refit <- .fit_weighted_logit(
        moments$regression_data$X[block$rows[keep], , drop = FALSE],
        block$D[keep],
        moments$regression_data$weight[block$rows[keep]],
        start = update$coefficients
      )
      cached_refit <- .deleted_propensity_refit(
        block, cluster, update$coefficients, moments$regression_data
      )
      expect_identical(cached_refit$status, old_refit$status)
      expect_equal(cached_refit$coefficients, old_refit$coefficients,
                   tolerance = 1e-10)

      old_overlap <- .propensity_overlap_status(
        moments$regression_data$X[block$rows[keep], , drop = FALSE],
        block$D[keep]
      )
      cached_overlap <- .deleted_overlap_status(
        block, moments$regression_data, keep, cluster
      )
      expect_identical(cached_overlap$status, old_overlap$status)
      expect_equal(cached_overlap$maximum, old_overlap$maximum,
                   tolerance = 1e-8)
    }
  }
})

test_that("cached exact nonlinear deletion is order invariant", {
  data <- make_nonlinear_fixture(H = 12L, observations_per_cell = 6L)
  moments <- collapse_rcs(
    data, "y", "time", "cohort", "region", weights = "w",
    xformla = ~x1 + x2 + female
  )
  att_map <- compile_att_map(moments)
  propensity <- .fit_propensity_blocks(moments, att_map, "not_yet",
                                       "warm_start")
  pairs <- list()
  for (k in seq_along(propensity$blocks)) {
    block <- propensity$blocks[[k]]
    for (cluster in unique(block$cluster)) {
      if (block$cluster_weight[cluster] > 0) {
        update <- .deleted_propensity_update(block, cluster)
        if (update$status == "valid") {
          pairs[[length(pairs) + 1L]] <- c(k, cluster)
        }
      }
    }
  }

  run_order <- function(order) {
    out <- vector("list", length(order))
    for (i in seq_along(order)) {
      pair <- pairs[[order[i]]]
      block <- propensity$blocks[[pair[1L]]]
      update <- .deleted_propensity_update(block, pair[2L])
      out[[i]] <- .deleted_propensity_refit(
        block, pair[2L], update$coefficients, moments$regression_data
      )$coefficients
    }
    names(out) <- vapply(order, function(i) paste(pairs[[i]], collapse = ":"),
                         character(1))
    out[sort(names(out))]
  }

  natural <- run_order(seq_along(pairs))
  reverse <- run_order(rev(seq_along(pairs)))
  set.seed(20260825)
  random <- run_order(sample(seq_along(pairs)))
  expect_equal(natural, reverse, tolerance = 0)
  expect_equal(natural, random, tolerance = 0)
})
