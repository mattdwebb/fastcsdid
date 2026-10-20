compare_fast_regression_to_did <- function(data, xformla, control,
                                           weighted = FALSE,
                                           tolerance = 1e-10) {
  fast <- suppressWarnings(fast_csdid(
    data, "y", "time", "cohort", "region",
    weights = if (weighted) "w" else NULL,
    xformla = xformla, est_method = "reg", control_group = control,
    inference = "none"
  ))
  external <- suppressMessages(suppressWarnings(
    did_reg_reference(data, xformla, control, weighted, compute_inffunc = TRUE)
  ))
  external_att <- data.frame(group = external$group, post = external$t,
                             estimate = external$att)
  external_att <- external_att[external_att$post >= external_att$group, ]
  joined <- merge(
    data.frame(group = fast$att_map$group, post = fast$att_map$post,
               estimate_fast = unname(fast$att_gt)),
    external_att, by = c("group", "post"), all = TRUE
  )
  expect_equal(joined$estimate_fast, joined$estimate, tolerance = tolerance)

  simple <- did::aggte(external, type = "simple", bstrap = FALSE, cband = FALSE)
  expect_equal(unname(fast$aggregation$estimates["simple"]), simple$overall.att,
               tolerance = tolerance)
  for (type in c("group", "calendar", "dynamic")) {
    result <- did::aggte(external, type = type, bstrap = FALSE, cband = FALSE)
    keep <- if (type == "dynamic") result$egt >= 0 else rep(TRUE, length(result$egt))
    prefix <- switch(type, group = "group:", calendar = "calendar:",
                     dynamic = "event:")
    labels <- paste0(prefix, result$egt[keep])
    expect_equal(unname(fast$aggregation$estimates[labels]), result$att.egt[keep],
                 tolerance = tolerance)
  }
  invisible(list(fast = fast, did = external, joined = joined))
}

test_that("optimized regression estimator directly matches did", {
  skip_if_not_installed("did", minimum_version = "2.5.1")
  data <- make_regression_fixture(unbalanced = TRUE, weighted = TRUE,
                                  missing_cells = TRUE)
  data$x_almost <- data$x1 + 1e-3 * sin(data$id)
  formulas <- list(
    intercept = ~1,
    numeric_binary = ~x1 + x2 + female,
    factor_interaction = ~x1 * sector + female,
    high_collinearity = ~x1 + x_almost + x2
  )
  for (xformla in formulas) {
    for (control in c("never", "not_yet")) {
      compare_fast_regression_to_did(data, xformla, control, weighted = TRUE)
    }
  }
})

test_that("complete-case and zero-weight semantics match did", {
  skip_if_not_installed("did", minimum_version = "2.5.1")
  data <- make_regression_fixture(unbalanced = TRUE, weighted = TRUE)
  data$x1[1L] <- NA_real_
  data$y[2L] <- NA_real_
  data$w[3L] <- 0
  comparison <- compare_fast_regression_to_did(
    data, ~x1 + x2 + female, "not_yet", weighted = TRUE
  )
  expect_equal(length(comparison$fast$moments$regression_data$dropped_rows), 2L)
  expect_equal(nrow(comparison$fast$moments$regression_data$X), nrow(data) - 2L)
})

test_that("rank-deficient full-sample regressions remain explicit", {
  skip_if_not_installed("did", minimum_version = "2.5.1")
  data <- make_regression_fixture()
  data$x_duplicate <- 2 * data$x1
  fast <- fast_csdid(
    data, "y", "time", "cohort", "region",
    xformla = ~x1 + x_duplicate, est_method = "reg", inference = "none"
  )
  external <- suppressMessages(suppressWarnings(
    did_reg_reference(data, ~x1 + x_duplicate, "never", FALSE)
  ))
  post <- external$t >= external$group
  expect_true(all(is.na(fast$att_gt)))
  expect_true(all(is.na(external$att[post])))
  expect_true(all(fast$att_support$status == "singular_outcome_regression"))
  expect_true(all(fast$regression$blocks$status == "singular_outcome_regression"))
})

test_that("moments-level regression fitting retains one aligned model matrix", {
  data <- make_regression_fixture(unbalanced = TRUE, weighted = TRUE)
  moments <- collapse_rcs(
    data, "y", "time", "cohort", "region", weights = "w",
    xformla = ~x1 * sector + female
  )
  fit <- fit_csdid_moments(moments, est_method = "reg", inference = "none")
  expect_equal(nrow(moments$regression_data$X), nrow(data))
  expect_equal(ncol(moments$regression_data$X),
               length(moments$regression_data$column_names))
  expect_equal(fit$regression$column_names,
               moments$regression_data$column_names)
  expect_error(
    fit_csdid_moments(collapse_rcs(data, "y", "time", "cohort", "region"),
                      est_method = "reg", inference = "none"),
    "requires moments created with xformla"
  )
})
