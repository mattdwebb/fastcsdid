test_that("production weighted logit matches transparent glm and DRDID", {
  data <- make_regression_fixture(unbalanced = TRUE, weighted = TRUE)
  prepared <- reference_nonlinear_prepare(data, ~x1 + x2 + female, "w")
  sample <- reference_att_sample(prepared, reference_reg_map(data)[2L, ],
                                 "never")
  production <- .fit_weighted_logit(sample$X, sample$D, sample$w)
  reference <- reference_logit(sample$X, sample$D, sample$w)
  expect_true(production$converged)
  expect_identical(production$status, "valid")
  expect_lt(production$score_norm, 1e-10)
  expect_equal(production$coefficients, reference$coefficients,
               tolerance = 1e-8)
  expect_equal(production$probability, reference$probabilities,
               tolerance = 1e-9)

  external <- suppressWarnings(fastglm::fastglmPure(
    x = sample$X, y = sample$D, family = stats::binomial(),
    weights = sample$w / mean(sample$w), method = 3L,
    tol = 1e-8, maxit = 100L
  ))
  expect_equal(production$coefficients, unname(external$coefficients),
               tolerance = 1e-7)
  expect_equal(production$probability, unname(external$fitted.values),
               tolerance = 1e-8)
})

test_that("master one-step starts preserve roots and reduce child iterations", {
  data <- make_regression_fixture(unbalanced = TRUE, weighted = TRUE)
  warm <- fast_csdid(
    data, "y", "time", "cohort", "region", weights = "w",
    xformla = ~x1 + x2 + female, est_method = "ipw",
    control_group = "never", inference = "none", ps_fit = "warm_start"
  )
  standard <- fast_csdid(
    data, "y", "time", "cohort", "region", weights = "w",
    xformla = ~x1 + x2 + female, est_method = "ipw",
    control_group = "never", inference = "none", ps_fit = "standard"
  )
  expect_equal(warm$att_gt, standard$att_gt, tolerance = 1e-10)
  expect_equal(warm$propensity$coefficients,
               standard$propensity$coefficients, tolerance = 1e-8)
  expect_true(all(warm$propensity$blocks$status == "valid"))
  expect_true(all(standard$propensity$blocks$status == "valid"))
  expect_lt(sum(warm$propensity$blocks$iterations),
            sum(standard$propensity$blocks$iterations))
  expect_true(all(warm$propensity$blocks$start_source == "master_one_step"))
})

test_that("logit failures are explicit and warm starts do not hide separation", {
  X <- cbind(1, c(-3, -2, -1, 1, 2, 3))
  D <- c(0, 0, 0, 1, 1, 1)
  standard <- .fit_weighted_logit(X, D, rep(1, 6))
  warm <- .fit_weighted_logit(X, D, rep(1, 6), start = c(0, 2))
  expect_true(standard$separation || standard$status != "valid")
  expect_true(warm$separation || warm$status != "valid")

  singular <- .fit_weighted_logit(cbind(1, 1:6, 2*(1:6)),
                                  c(0, 1, 0, 1, 0, 1), rep(1, 6))
  expect_identical(singular$status, "singular_propensity_information")
})

test_that("one-step sign reduces the exact child score", {
  data <- make_regression_fixture(unbalanced = TRUE, weighted = TRUE)
  moments <- collapse_rcs(
    data, "y", "time", "cohort", "region", weights = "w",
    xformla = ~x1 + x2 + female
  )
  map <- compile_att_map(moments)
  masters <- .fit_propensity_master_blocks(moments, map, "never")
  sample <- .propensity_block_data(
    moments$regression_data, moments, map[1L, ], "never"
  )
  X <- moments$regression_data$X[sample$rows, , drop = FALSE]
  w <- moments$regression_data$weight[sample$rows]
  start <- masters$blocks[[masters$cell_master[1L]]]$coefficients
  before <- .weighted_logit_state(X, sample$D, w, start)
  one_step <- .one_logit_newton_step(X, sample$D, w, start)
  after <- .weighted_logit_state(X, sample$D, w, one_step)
  expect_lt(max(abs(after$score)), max(abs(before$score)))
})
