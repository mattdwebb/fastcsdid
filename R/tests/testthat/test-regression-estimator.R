test_that("optimized full-sample regression estimates match the independent reference", {
    designs <- list(balanced = make_regression_fixture(), unbalanced = make_regression_fixture(unbalanced = TRUE),
        missing_cells = make_regression_fixture(unbalanced = TRUE, missing_cells = TRUE), weighted = make_regression_fixture(unbalanced = TRUE,
            weighted = TRUE))
    formulas <- list(intercept = ~1, numeric = ~x1 + x2 + female, factor = ~x1 * sector + female)
    for (design_name in names(designs)) {
        data <- designs[[design_name]]
        weighted <- identical(design_name, "weighted")
        for (formula_name in names(formulas)) {
            xformla <- formulas[[formula_name]]
            for (control in c("never", "not_yet")) {
                fit <- fast_csdid(data, "y", "time", "cohort", "region", weights = if (weighted)
                  "w"
                else NULL, control_group = control, inference = "none", xformla = xformla, est_method = "reg")
                expected <- reference_reg_fit(data, xformla, control, weights = if (weighted)
                  "w"
                else NULL, att_map = fit$att_map)
                expect_equal(fit$att_gt, expected$att, tolerance = 1e-10, info = paste(design_name, formula_name,
                  control, "ATT"))
                expect_equal(fit$aggregation$estimates, expected$aggregate, tolerance = 1e-10, info = paste(design_name,
                  formula_name, control, "aggregate"))
            }
        }
    }
})
test_that("intercept-only regression reduces to the no-covariate estimator", {
    data <- make_regression_fixture(unbalanced = TRUE, weighted = TRUE, missing_cells = TRUE)
    for (control in c("never", "not_yet")) {
        unadjusted <- fast_csdid(data, "y", "time", "cohort", "region", weights = "w", control_group = control,
            inference = "none")
        adjusted <- fast_csdid(data, "y", "time", "cohort", "region", weights = "w", control_group = control,
            inference = "none", xformla = ~1, est_method = "reg")
        expect_equal(adjusted$att_gt, unadjusted$att_gt, tolerance = 1e-12)
        expect_equal(adjusted$aggregation$estimates, unadjusted$aggregation$estimates, tolerance = 1e-12)
    }
})
test_that("regression blocks are shared exactly by period and control risk set", {
    data <- make_regression_fixture(unbalanced = TRUE)
    never <- fast_csdid(data, "y", "time", "cohort", "region", xformla = ~x1 + x2, est_method = "reg",
        control_group = "never", inference = "none")
    used_periods <- sort(unique(c(never$att_map$base, never$att_map$post)))
    expect_equal(nrow(never$regression$blocks), length(used_periods))
    expect_equal(sort(never$regression$blocks$period), used_periods)
    expect_true(all(table(never$regression$blocks$period) == 1L))
    not_yet <- fast_csdid(data, "y", "time", "cohort", "region", xformla = ~x1 + x2, est_method = "reg",
        control_group = "not_yet", inference = "none")
    for (post in unique(not_yet$att_map$post)) {
        cells <- which(not_yet$att_map$post == post)
        expect_equal(length(unique(not_yet$regression$cell_blocks[cells, "post_block"])), 1L)
    }
    keys <- paste(not_yet$regression$blocks$period, not_yet$regression$blocks$control_regions)
    expect_equal(anyDuplicated(keys), 0L)
})
