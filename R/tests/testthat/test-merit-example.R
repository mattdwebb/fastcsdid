test_that("MERIT paper example is stable under changing-support jackknife", {
  path <- system.file("extdata", "merit.csv", package = "fastrcsdid")
  d <- utils::read.csv(path)
  d$cohort[d$cohort == 0] <- Inf

  fit <- fast_csdid(
    d, "outcome", "year", "cohort", "state", weights = "weight",
    est_method = "none", inference = "jackknife"
  )

  expect_equal(unname(fit$aggregation$estimates["simple"]),
               0.0477819967004955, tolerance = 1e-12)
  expect_equal(unname(fit$jackknife$aggregate_covariance$se_cv3j["simple"]),
               0.0146168254236831, tolerance = 1e-12)
  expect_equal(nrow(fit$jackknife$aggregate), 51L)
  expect_equal(sum(fit$jackknife$support_changed[, "simple"]), 4L)
  expect_true(fit$jackknife$aggregate_covariance$complete_by_stat["simple"])
})
