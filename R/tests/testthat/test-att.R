test_that("ATT cells equal direct four-cell means", {
  d <- make_test_data(unbalanced = TRUE)
  fit <- fast_csdid(d, "y", "time", "cohort", "region", inference = "none")
  expected <- mapply(
    manual_att,
    g = fit$att_map$group,
    base = fit$att_map$base,
    post = fit$att_map$post,
    MoreArgs = list(d = d, control = "never")
  )
  expect_equal(unname(fit$att_gt), unname(expected), tolerance = 1e-12)
})

test_that("not-yet-treated comparisons use late cohorts before adoption", {
  d <- make_test_data(unbalanced = TRUE)
  fit <- fast_csdid(d, "y", "time", "cohort", "region",
                    control_group = "not_yet", inference = "none")
  expected <- mapply(
    manual_att,
    g = fit$att_map$group,
    base = fit$att_map$base,
    post = fit$att_map$post,
    MoreArgs = list(d = d, control = "not_yet")
  )
  expect_equal(unname(fit$att_gt), unname(expected), tolerance = 1e-12)
})
