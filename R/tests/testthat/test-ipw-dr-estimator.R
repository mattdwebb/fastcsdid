test_that("optimized IPW and DR point estimates match the independent reference", {
    designs <- list(balanced = make_regression_fixture(), weighted = make_regression_fixture(unbalanced = TRUE,
        weighted = TRUE, missing_cells = TRUE))
    formulas <- list(intercept = ~1, numeric = ~x1 + x2 + female, interaction = ~x1 * sector + female)
    for (design_name in names(designs)) {
        data <- designs[[design_name]]
        weighted <- design_name == "weighted"
        for (method in c("ipw", "dr")) {
            for (formula_name in names(formulas)) {
                if (weighted && method == "dr" && formula_name == "interaction")
                  next
                for (control in c("never", "not_yet")) {
                  fit <- fast_csdid(data, "y", "time", "cohort", "region", weights = if (weighted)
                    "w"
                  else NULL, xformla = formulas[[formula_name]], est_method = method, control_group = control,
                    inference = "none")
                  expected <- reference_nonlinear_fit(data, formulas[[formula_name]], method, control,
                    weights = if (weighted)
                      "w"
                    else NULL, att_map = fit$att_map)
                  expect_equal(fit$att_gt, expected$att, tolerance = 1e-09, info = paste(design_name,
                    method, formula_name, control))
                  expect_equal(fit$aggregation$estimates, expected$aggregate, tolerance = 1e-09, info = paste(design_name,
                    method, formula_name, control, "aggregate"))
                }
            }
        }
    }
})
test_that("optimized IPW and DR directly match did point estimates and aggregations", {
    skip_if_not_installed("did", minimum_version = "2.5.1")
    data <- make_regression_fixture(unbalanced = TRUE, weighted = TRUE, missing_cells = TRUE)
    for (method in c("ipw", "dr")) {
        for (control in c("never", "not_yet")) {
            fit <- fast_csdid(data, "y", "time", "cohort", "region", weights = "w", xformla = ~x1 + x2 +
                female, est_method = method, control_group = control, inference = "none")
            external <- suppressMessages(suppressWarnings(did_nonlinear_reference(data, ~x1 + x2 + female,
                method, control, weighted = TRUE, compute_inffunc = TRUE)))
            external_att <- data.frame(group = external$group, post = external$t, att = external$att)
            external_att <- external_att[external_att$post >= external_att$group, , drop = FALSE]
            joined <- merge(data.frame(group = fit$att_map$group, post = fit$att_map$post, att = unname(fit$att_gt)),
                external_att, by = c("group", "post"), suffixes = c("_fast", "_did"))
            expect_equal(joined$att_fast, joined$att_did, tolerance = 2e-09)
            expect_equal(fit$aggregation$estimates, reference_did_aggregates(external), tolerance = 2e-09)
        }
    }
})
test_that("intercept-only IPW and DR reduce to the validated estimator", {
    data <- make_regression_fixture(unbalanced = TRUE, weighted = TRUE, missing_cells = TRUE)
    for (control in c("never", "not_yet")) {
        baseline <- fast_csdid(data, "y", "time", "cohort", "region", weights = "w", control_group = control,
            inference = "none")
        for (method in c("ipw", "dr")) {
            fit <- fast_csdid(data, "y", "time", "cohort", "region", weights = "w", xformla = ~1, est_method = method,
                control_group = control, inference = "none")
            expect_equal(fit$att_gt, baseline$att_gt, tolerance = 1e-10)
            expect_equal(fit$aggregation$estimates, baseline$aggregation$estimates, tolerance = 1e-10)
        }
    }
})
test_that("propensity and DR outcome block maps expose genuine sharing", {
    data <- make_regression_fixture(unbalanced = TRUE, weighted = TRUE)
    fit <- fast_csdid(data, "y", "time", "cohort", "region", weights = "w", xformla = ~x1 + x2 + female,
        est_method = "dr", control_group = "never", inference = "none")
    expect_equal(nrow(fit$propensity$blocks), nrow(fit$att_map))
    expect_lt(nrow(fit$propensity$masters), nrow(fit$propensity$blocks))
    expect_equal(nrow(fit$regression$control_blocks), length(unique(c(fit$att_map$base, fit$att_map$post))))
    expected_treated <- nrow(unique(data.frame(group = rep(fit$att_map$group, 2L), period = c(fit$att_map$base,
        fit$att_map$post))))
    expect_equal(nrow(fit$regression$treated_blocks), expected_treated)
})
