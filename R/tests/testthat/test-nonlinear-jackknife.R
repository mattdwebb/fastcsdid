test_that("sparse exact nonlinear refits match physical deletion", {
  designs <- list(
    unweighted = make_regression_fixture(unbalanced = TRUE,
                                         missing_cells = TRUE),
    weighted = make_regression_fixture(unbalanced = TRUE, weighted = TRUE,
                                       missing_cells = TRUE)
  )
  designs$weighted$w[designs$weighted$id %% 19L == 0L] <- 0
  designs$weighted$w[
    designs$weighted$region == 2 &
      designs$weighted$time == max(designs$weighted$time)
  ] <- 100
  for (design_name in names(designs)) {
    data <- designs[[design_name]]
    weights <- if (design_name == "weighted") "w" else NULL
    for (method in c("ipw", "dr")) {
      for (control in c("never", "not_yet")) {
        fit <- fast_csdid(
          data, "y", "time", "cohort", "region", weights = weights,
          xformla = ~x1 + female, est_method = method,
          control_group = control, inference = "jackknife",
          jackknife_method = "refit"
        )
        for (cluster in unique(data$region)) {
          deleted <- data[data$region != cluster, , drop = FALSE]
          expected <- reference_nonlinear_fit(
            deleted, ~x1 + female, method, control, weights,
            att_map = fit$att_map
          )
          expect_equal(
            fit$jackknife$att[as.character(cluster), ], expected$att,
            tolerance = 1e-9,
            info = paste(design_name, method, control, cluster, "ATT")
          )
          expect_equal(
            fit$jackknife$aggregate[as.character(cluster), ],
            expected$aggregate, tolerance = 1e-9,
            info = paste(design_name, method, control, cluster, "aggregate")
          )
          expect_identical(
            is.na(fit$jackknife$att[as.character(cluster), ]),
            is.na(expected$att)
          )
        }
      }
    }
  }
})

test_that("deleted-Jacobian sign and block recursion equal a dense solve", {
  data <- make_regression_fixture(unbalanced = TRUE, weighted = TRUE)
  moments <- collapse_rcs(
    data, "y", "time", "cohort", "region", weights = "w",
    xformla = ~x1 + x2 + female
  )
  map <- compile_att_map(moments)
  propensity <- .fit_propensity_blocks(moments, map, "never", "warm_start")
  block <- propensity$blocks[[1L]]
  cluster <- block$treated_regions[1L]
  update <- .deleted_propensity_update(block, cluster)
  expect_identical(update$status, "valid")
  keep <- block$cluster != cluster
  deleted_state_at_full <- .weighted_logit_state(
    moments$regression_data$X[block$rows[keep], , drop = FALSE],
    block$D[keep], moments$regression_data$weight[block$rows[keep]],
    block$coefficients
  )
  expect_equal(deleted_state_at_full$score,
               -block$cluster_score[cluster, ], tolerance = 1e-10)

  direct <- .evaluate_nonlinear_sample(
    block, moments$regression_data, block$coefficients, "ipw",
    keep_rows = keep
  )
  p <- length(block$coefficients)
  jacobian <- matrix(0, p + 1L, p + 1L)
  jacobian[seq_len(p), seq_len(p)] <- -update$information
  jacobian[p + 1L, seq_len(p)] <- -direct$gradient
  jacobian[p + 1L, p + 1L] <- 1
  equations <- c(deleted_state_at_full$score, 0)
  dense_change <- solve(jacobian, -equations)
  expect_equal(dense_change[seq_len(p)], update$change, tolerance = 1e-10)
  expect_equal(dense_change[p + 1L],
               as.numeric(crossprod(direct$gradient, update$change)),
               tolerance = 1e-10)
})

test_that("linearized and exact-refit deletions agree closely in regular designs", {
  data <- make_nonlinear_fixture()
  for (method in c("ipw", "dr")) {
    linearized <- fast_csdid(
      data, "y", "time", "cohort", "region", weights = "w",
      xformla = ~x1 + female, est_method = method,
      control_group = "never", inference = "jackknife",
      jackknife_method = "linearized"
    )
    refit <- fast_csdid(
      data, "y", "time", "cohort", "region", weights = "w",
      xformla = ~x1 + female, est_method = method,
      control_group = "never", inference = "jackknife",
      jackknife_method = "refit"
    )
    common <- is.finite(linearized$jackknife$att) &
      is.finite(refit$jackknife$att)
    difference <- abs(linearized$jackknife$att[common] -
                        refit$jackknife$att[common])
    expect_gt(length(difference), 0L)
    expect_lt(stats::median(difference), 0.02)
    linear_change <- linearized$jackknife$att_changes[common]
    exact_change <- refit$jackknife$att_changes[common]
    expect_gt(stats::cor(linear_change, exact_change), 0.95)

    common_se <- is.finite(linearized$jackknife$att_covariance$pairwise_se_cv3) &
      is.finite(refit$jackknife$att_covariance$pairwise_se_cv3)
    expect_lt(max(abs(
      linearized$jackknife$att_covariance$pairwise_se_cv3[common_se] -
        refit$jackknife$att_covariance$pairwise_se_cv3[common_se]
    )), 0.1)
  }
})

test_that("nonlinear dependency graph leaves unrelated entries exact", {
  data <- make_regression_fixture(unbalanced = TRUE, weighted = TRUE)
  fit <- fast_csdid(
    data, "y", "time", "cohort", "region", weights = "w",
    xformla = ~x1 + female, est_method = "ipw",
    control_group = "not_yet", inference = "jackknife",
    jackknife_method = "linearized"
  )
  unaffected <- fit$jackknife$att_status_label == "unaffected"
  full <- matrix(fit$att_gt, nrow(unaffected), ncol(unaffected), byrow = TRUE)
  expect_true(any(unaffected))
  expect_identical(unname(fit$jackknife$att[unaffected]),
                   unname(full[unaffected]))
  expect_identical(unname(fit$jackknife$att_influence[unaffected]),
                   rep(0, sum(unaffected)))
})

test_that("nonlinear diagnostics expose method, scores, conditions, and sources", {
  data <- make_regression_fixture(unbalanced = TRUE, weighted = TRUE)
  fit <- fast_csdid(
    data, "y", "time", "cohort", "region", weights = "w",
    xformla = ~x1 + female, est_method = "dr",
    control_group = "not_yet", inference = "jackknife",
    jackknife_method = "linearized"
  )
  table <- jackknife_cluster_table(fit, names(fit$att_gt)[1L])
  expect_true(all(table$jackknife_method == "linearized"))
  expect_true(any(is.finite(table$max_propensity_coefficient_change)))
  expect_true(any(is.finite(table$deleted_score_norm)))
  expect_true(any(is.finite(table$deleted_information_condition)))
  expect_true(any(is.finite(table$direct_moment_change)))
  expect_true(any(is.finite(table$outcome_regression_change)))
  expect_true(any(is.finite(table$propensity_change)))
  summary <- jackknife_influence_summary(fit)
  expect_true(all(summary$jackknife_method == "linearized"))
  expect_true(any(is.finite(summary$max_propensity_coefficient_change)))
})
