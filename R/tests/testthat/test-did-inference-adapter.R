test_that("common adapter retains ordinary cluster t inference", {
  for (unequal in c(FALSE, TRUE)) {
    d <- make_test_data(unequal)
    d$y <- d$y + .08 * cos(d$region * d$time)
    design <- did_design(d, "y", "region", "time", "cohort")
    fast <- did_jackknife(design)
    slow <- did_jackknife(design, engine = "reference")
    expect_equal(fast$delete_estimates, slow$delete_estimates, tolerance = 1e-10)
    expect_equal(fast$jackknife_df, 7)
    expect_equal(fast$p_value, 2 * pt(-abs(fast$statistic), 7))
  }
})
