test_that("v1 exposes only retained inference methods", {
  d <- fastrcsdid_example
  for (method in c("ri", "both", "rijack", "randjack")) {
    expect_error(fast_csdid(d, "y", "time", "cohort", "region",
                          inference = method), "arg")
  }
  for (option in c("ri_reps", "ri_seed", "ri_assignments", "ri_strata",
                   "ri_subtract", "ri_batch_size", "keep_ri_distribution",
                   "store_ri_distribution")) {
    args <- list(d, "y", "time", "cohort", "region", inference = "none")
    args[[option]] <- 1
    expect_error(do.call(fast_csdid, args), "unused argument")
  }
  ns <- asNamespace("fastrcsdid")
  for (symbol in c("randomization_inference", "prepare_ri_assignments",
                   "did_randomization_inference", "did_assignments", "did_placebo")) {
    expect_false(exists(symbol, envir = ns, inherits = FALSE))
  }
  fit <- fast_csdid(d, "y", "time", "cohort", "region", inference = "jackknife")
  expect_false("ri" %in% names(fit))
  expect_true(all(is.finite(fit$jackknife$aggregate_covariance$se_cv3j)))
})
