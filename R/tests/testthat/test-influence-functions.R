test_that("observation ATT influence functions match did after exact alignment", {
  skip_if_not_installed("did", minimum_version = "2.5.1")
  data <- make_nonlinear_fixture(H = 20L, observations_per_cell = 8L)
  specifications <- list(
    none = list(local = "none", external = "reg", formula = ~1),
    reg = list(local = "reg", external = "reg", formula = ~x1 + female),
    ipw = list(local = "ipw", external = "ipw", formula = ~x1 + female),
    dr = list(local = "dr", external = "dr", formula = ~x1 + female)
  )
  for (name in names(specifications)) {
    specification <- specifications[[name]]
    for (control in c("never", "not_yet")) {
      fit <- fast_csdid(
        data, "y", "time", "cohort", "region", weights = "w",
        xformla = if (name == "none") NULL else specification$formula,
        est_method = specification$local, control_group = control,
        inference = "none"
      )
      external <- did_nonlinear_reference(
        data, specification$formula, specification$external, control,
        weighted = TRUE, compute_inffunc = TRUE, clustered = TRUE
      )
      columns <- match(
        paste(fit$att_map$group, fit$att_map$post, sep = "|"),
        paste(external$group, external$t, sep = "|")
      )
      internal <- as.data.frame(external$DIDparams$data)
      expected <- matrix(0, nrow(data), length(columns))
      expected[internal$.rowid, ] <- as.matrix(
        external$inffunc[, columns, drop = FALSE]
      )
      observed <- influence_function(fit, type = "att")
      tolerance <- if (name %in% c("ipw", "dr")) 1e-8 else 1e-11
      expect_equal(observed, expected, tolerance = tolerance,
                   ignore_attr = TRUE, info = paste(name, control))
      expect_lt(max(abs(colSums(observed))), 1e-7)
      cluster_expected <- rowsum(expected, data$region, reorder = TRUE)
      expect_equal(
        influence_function(fit, type = "att", cluster_sum = TRUE),
        cluster_expected, tolerance = tolerance, ignore_attr = TRUE
      )
      expect_equal(
        fit$vcov$analytical_raw[names(fit$att_gt), names(fit$att_gt)],
        crossprod(cluster_expected) / nrow(data)^2,
        tolerance = 1e-10, ignore_attr = TRUE
      )
    }
  }
})

test_that("aggregate influence includes cohort-mass weight estimation", {
  data <- make_nonlinear_fixture(H = 20L, observations_per_cell = 8L)
  fit <- fast_csdid(
    data, "y", "time", "cohort", "region", weights = "w",
    xformla = ~x1 + female, est_method = "dr", control_group = "not_yet",
    inference = "none"
  )
  att_if <- influence_function(fit, type = "att")
  aggregate_if <- influence_function(fit, type = "aggregate")
  for (a in seq_along(fit$aggregation$estimates)) {
    fixed_part <- as.numeric(att_if %*% fit$aggregation$weights[a, ])
    if (fit$aggregation$metadata$type[a] == "group") {
      expect_equal(aggregate_if[, a], fixed_part, tolerance = 1e-11)
    } else {
      expect_equal(
        aggregate_if[, a] - fixed_part,
        .aggregate_weight_observation_if(fit$influence$factorized, a),
        tolerance = 1e-11
      )
    }
  }
  cluster <- influence_function(fit, type = "aggregate", cluster_sum = TRUE)
  expect_equal(cluster, rowsum(aggregate_if, data$region, reorder = TRUE),
               tolerance = 1e-9, ignore_attr = TRUE)
  expect_lt(max(abs(colSums(aggregate_if))), 1e-7)
})

test_that("unused unsupported ATT cells do not contaminate aggregate IFs", {
  data <- did_fixture(unbalanced = TRUE, weighted = TRUE)
  data <- data[!(data$cohort == 2006 & data$time == 2015), , drop = FALSE]
  fit <- fast_csdid(
    data, "y", "time", "cohort", "region", weights = "w",
    inference = "none"
  )
  expect_false(fit$att_support$supported[
    fit$att_map$group == 2006 & fit$att_map$post == 2015
  ])
  expect_true(is.finite(fit$aggregation$estimates["calendar:2010"]))
  observation <- influence_function(
    fit, "calendar:2010", "aggregate"
  )
  cluster <- influence_function(
    fit, "calendar:2010", "aggregate", cluster_sum = TRUE
  )
  expect_true(all(is.finite(cluster)))
  expect_equal(
    cluster, rowsum(observation, data$region, reorder = TRUE),
    tolerance = 1e-11, ignore_attr = TRUE
  )
})

test_that("dense stacked DR observation solve equals block recursion", {
  data <- make_nonlinear_fixture(H = 12L, observations_per_cell = 6L)
  fit <- fast_csdid(
    data, "y", "time", "cohort", "region", weights = "w",
    xformla = ~x1 + female, est_method = "dr", inference = "none"
  )
  factorized <- fit$influence$factorized
  k <- 1L
  block <- factorized$evaluation$propensity$blocks[[k]]
  regression <- fit$moments$regression_data
  rows <- block$rows
  X <- regression$X[rows, , drop = FALSE]
  y <- regression$y[rows]
  D <- block$D
  post <- block$post
  w <- regression$weight[rows]
  w <- w / mean(w)
  pscore <- pmin(stats::plogis(as.numeric(X %*% block$coefficients)), 1-1e-6)
  coefficients <- .nonlinear_outcome_coefficients(
    k, factorized$evaluation$control_nuisance,
    factorized$evaluation$treated_nuisance
  )
  m00 <- as.numeric(X %*% coefficients$control_base)
  m01 <- as.numeric(X %*% coefficients$control_post)
  m10 <- as.numeric(X %*% coefficients$treated_base)
  m11 <- as.numeric(X %*% coefficients$treated_post)
  m0 <- post*m01 + (1-post)*m00
  trim <- D == 1 | pscore < .995
  odds <- pscore/(1-pscore)
  w10 <- trim*w*D*(1-post); w11 <- trim*w*D*post
  w00 <- trim*w*odds*(1-D)*(1-post); w01 <- trim*w*odds*(1-D)*post
  wd <- trim*w*D; wd1 <- wd*post; wd0 <- wd*(1-post)
  r10 <- .if_ratio(w10, y-m0); r11 <- .if_ratio(w11, y-m0)
  r00 <- .if_ratio(w00, y-m0); r01 <- .if_ratio(w01, y-m0)
  edp <- .if_ratio(wd, m11-m01); e1p <- .if_ratio(wd1, m11-m01)
  ed0 <- .if_ratio(wd, m10-m00); e00 <- .if_ratio(wd0, m10-m00)
  direct <- (r11$influence-r10$influence) -
    (r01$influence-r00$influence) +
    (edp$influence-e1p$influence) - (ed0$influence-e00$influence)
  n <- length(y); q <- ncol(X)
  score <- cbind(
    X*(w*(D-pscore)),
    X*(w*(1-D)*(1-post)*(y-m00)),
    X*(w*(1-D)*post*(y-m01)),
    X*(w*D*(1-post)*(y-m10)),
    X*(w*D*post*(y-m11))
  )
  information <- matrix(0, 5*q, 5*q)
  weights <- list(
    w*pscore*(1-pscore), w*(1-D)*(1-post), w*(1-D)*post,
    w*D*(1-post), w*D*post
  )
  for (j in seq_len(5L)) {
    index <- ((j-1L)*q+1L):(j*q)
    information[index, index] <- crossprod(X, weights[[j]]*X)/n
  }
  dense_representation <- .if_solve_rows(information, score)
  M1post <- -as.numeric(crossprod(w11, X))/n/mean(w11)
  M1pre <- -as.numeric(crossprod(w10, X))/n/mean(w10)
  M2pre <- as.numeric(crossprod(w00*(y-m0-r00$estimate), X))/n/mean(w00)
  M2post <- as.numeric(crossprod(w01*(y-m0-r01$estimate), X))/n/mean(w01)
  M3post <- -as.numeric(crossprod(w01, X))/n/mean(w01)
  M3pre <- -as.numeric(crossprod(w00, X))/n/mean(w00)
  moment_post <- as.numeric(crossprod(wd/mean(wd)-wd1/mean(wd1), X))/n
  moment_pre <- as.numeric(crossprod(wd/mean(wd)-wd0/mean(wd0), X))/n
  loading <- c(
    -(M2post-M2pre),
    M1pre-M3pre+moment_pre,
    M1post-M3post-moment_post,
    -moment_pre,
    moment_post
  )
  dense_if <- direct + as.numeric(dense_representation %*% loading)
  recursive_if <- .local_dr_if(
    y, post, D, X, regression$weight[rows], pscore, coefficients
  )
  expect_equal(recursive_if, dense_if, tolerance = 1e-10)
})

test_that("RIF generation is lazy, batched, and preserves retained row order", {
  data <- make_nonlinear_fixture(H = 20L, observations_per_cell = 8L)
  fit <- fast_csdid(
    data, "y", "time", "cohort", "region", weights = "w",
    xformla = ~x1 + female, est_method = "ipw", inference = "none"
  )
  expect_false(any(vapply(fit$influence$factorized, function(x) {
    is.matrix(x) && identical(dim(x), c(nrow(data), length(fit$att_gt)))
  }, logical(1))))
  one <- influence_function(fit, names(fit$att_gt)[1L], "att")
  expect_identical(dim(one), c(nrow(data), 1L))
  all_1 <- influence_function(fit, type = "att", batch_size = 1L)
  all_5 <- influence_function(fit, type = "att", batch_size = 5L)
  expect_equal(all_1, all_5)
  recentered <- rif(fit, names(fit$att_gt)[1L])
  expect_equal(colMeans(recentered), fit$att_gt[1L], tolerance = 1e-10,
               ignore_attr = TRUE)
  metadata <- attr(recentered, "observation_metadata")
  expect_identical(metadata$original_row, seq_len(nrow(data)))
  expect_identical(metadata$cluster, data$region)
  expect_equal(metadata$period, data$time)
  simple <- rif(fit, "simple")
  expect_equal(colMeans(simple), fit$aggregation$estimates["simple"],
               tolerance = 1e-10, ignore_attr = TRUE)
})

test_that("moments-only no-covariate fits retain cluster influence only", {
  data <- make_regression_fixture()
  moments <- collapse_rcs(data, "y", "time", "cohort", "region")
  fit <- fit_csdid_moments(moments, inference = "none")
  expect_identical(dim(fit$influence$cluster_att),
                   c(length(moments$regions), nrow(fit$att_map)))
  expect_error(influence_function(fit, type = "att"),
               "without retained microdata")
  expect_silent(influence_function(fit, type = "att", cluster_sum = TRUE))
})

test_that("analytical and both nonlinear jackknives are compared together", {
  data <- make_nonlinear_fixture(H = 20L, observations_per_cell = 8L)
  fit <- fast_csdid(
    data, "y", "time", "cohort", "region", weights = "w",
    xformla = ~x1 + female, est_method = "dr", inference = "jackknife",
    jackknife_method = "both"
  )
  expect_named(fit$jackknife$methods, c("linearized", "refit"))
  comparison <- inference_comparison(fit)
  expect_true(all(is.finite(comparison$analytical_se)))
  expect_true(any(is.finite(comparison$linearized_cv3_se)))
  expect_true(any(is.finite(comparison$refit_cv3_se)))
  expect_true(any(is.finite(comparison$linearized_refit_ratio)))
  table <- jackknife_cluster_table(fit, names(fit$att_gt)[1L])
  expect_true(any(is.finite(table$linearized_refit_difference)))
  expect_true(any(is.finite(table$cluster_score_norm)))
  expect_true(all(table$refit_valid))
  expect_true(all(table$refit_status %in% c("updated", "unaffected")))
  summary <- jackknife_influence_summary(fit)
  expect_true(any(is.finite(summary$max_abs_linearized_refit_difference)))
  expect_true(all(summary$n_refit_failures == 0L))
  expect_true(all(summary$refit_failure_reasons == ""))
})

test_that("exact linear deletion populates primary comparison columns", {
  data <- make_regression_fixture()
  fit <- fast_csdid(
    data, "y", "time", "cohort", "region", inference = "jackknife"
  )
  comparison <- inference_comparison(fit)
  expect_true(all(is.finite(comparison$analytical_se)))
  expect_true(all(is.finite(comparison$linearized_cv3_se)))
  expect_true(all(is.na(comparison$refit_cv3_se)))
  expect_identical(
    fit$vcov$analytical_metadata,
    fit$vcov$metadata
  )
})
