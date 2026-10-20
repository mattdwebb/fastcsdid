test_that("singleton deletions invalidate only aggregates requiring the removed cohort", {
  d <- make_test_data(TRUE)
  d <- d[d$region != 2, , drop = FALSE]
  fit <- fast_csdid(d, "y", "time", "cohort", "region", inference = "jackknife",
                    jackknife_support = "fixed")
  expect_identical(fit$settings$jackknife_support, "fixed")
  expect_true(is.na(fit$jackknife$aggregate["1", "simple"]))
  expect_true(fit$jackknife$support_changed["1", "simple"])
  expect_match(fit$jackknife$aggregate_reason["1", "simple"], "support")
  expect_true(is.na(fit$jackknife$aggregate_covariance$se_cv3j["simple"]))
  expect_false(fit$jackknife$aggregate_covariance$complete_by_stat["simple"])
  expect_true(is.finite(fit$jackknife$aggregate["1", "group:5"]))
  expect_true(is.na(fit$jackknife$aggregate["1", "group:3"]))

  historical <- update(fit, jackknife_support = "changing")
  expect_equal(fit$aggregation$estimates, historical$aggregation$estimates)
  expect_equal(fit$vcov$analytical, historical$vcov$analytical)
  expect_true(is.finite(historical$jackknife$aggregate["1", "simple"]))
  expect_identical(historical$jackknife$support, "changing")
  expect_true(historical$jackknife$support_changed["1", "simple"])
})

test_that("all four methods enforce fixed support and retain the full-sample fit", {
  d <- make_test_data(TRUE)
  d <- d[d$region != 2, , drop = FALSE]
  d$x <- sin(seq_len(nrow(d)) / 3) + d$region / 11
  for (method in c("none", "reg", "ipw", "dr")) {
    fit <- fast_csdid(d, "y", "time", "cohort", "region", est_method = method,
      xformla = if (method == "none") NULL else ~x,
      aggregation = "simple", inference = "jackknife", jackknife_method = "refit",
      jackknife_support = "fixed")
    expect_true(is.finite(fit$aggregation$estimates["simple"]), info = method)
    expect_true(is.na(fit$jackknife$aggregate["1", "simple"]), info = method)
    expect_false(fit$jackknife$aggregate_covariance$complete_by_stat["simple"], info = method)
  }
})

test_that("support-preserving results do not change under the policy switch", {
  d <- make_test_data(TRUE)
  fixed <- fast_csdid(d, "y", "time", "cohort", "region", inference = "jackknife",
                      jackknife_support = "fixed")
  changing <- update(fixed, jackknife_support = "changing")
  expect_equal(fixed$jackknife$aggregate, changing$jackknife$aggregate)
  expect_equal(fixed$jackknife$aggregate_covariance, changing$jackknife$aggregate_covariance)
  expect_true(all(fixed$jackknife$aggregate_covariance$complete_by_stat))
  expect_error(update(fixed, jackknife_support = "silent"), "arg")
})

test_that("default reaggregates after singleton deletion for all four methods", {
  d <- make_test_data(TRUE)
  d <- d[d$region != 2, , drop = FALSE]
  d$x <- sin(seq_len(nrow(d)) / 3) + d$region / 11
  for (method in c("none", "reg", "ipw", "dr")) {
    fit <- fast_csdid(d, "y", "time", "cohort", "region", est_method = method,
      xformla = if (method == "none") NULL else ~x,
      aggregation = "simple", inference = "jackknife", jackknife_method = "refit")
    explicit <- update(fit, jackknife_support = "changing")
    expect_identical(fit$settings$jackknife_support, "changing")
    expect_equal(fit$jackknife$aggregate, explicit$jackknife$aggregate)
    expect_true(fit$jackknife$support_changed["1", "simple"])
    expect_true(is.finite(fit$jackknife$aggregate["1", "simple"]))
    expect_true(fit$jackknife$aggregate_covariance$complete_by_stat["simple"])
    expect_true(is.finite(fit$jackknife$aggregate_covariance$se_cv3j["simple"]))
    refit <- fast_csdid(d[d$region != 1, ], "y", "time", "cohort", "region",
      est_method = method, xformla = if (method == "none") NULL else ~x,
      aggregation = "simple", inference = "none")
    expect_equal(unname(fit$jackknife$aggregate["1", "simple"]),
                 unname(refit$aggregation$estimates["simple"]), tolerance = 1e-7)
  }
})
