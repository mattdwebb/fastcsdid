test_that("all ATT and aggregate deletions match an independent physical reference", {
    designs <- list(balanced = make_test_data(FALSE), unequal_cells = make_test_data(TRUE), missing_cells = subset(make_test_data(TRUE),
        !(region == 2 & time %in% c(2, 6))), unequal_clusters = subset(make_test_data(TRUE), !(region ==
        7 & time %in% c(1, 2, 3))))
    controls <- c("never", "not_yet")
    for (design_name in names(designs)) {
        d <- designs[[design_name]]
        for (control in controls) {
            fit <- fast_csdid(d, "y", "time", "cohort", "region", control_group = control, inference = "jackknife")
            full_reference <- reference_fit(d, fit$att_map, control)
            expect_equal(fit$att_gt, full_reference$att, tolerance = 1e-12, info = paste(design_name,
                control, "full ATT"))
            expect_equal(fit$aggregation$estimates, full_reference$aggregate, tolerance = 1e-12, info = paste(design_name,
                control, "full aggregate"))
            for (h in unique(d$region)) {
                deleted_data <- d[d$region != h, , drop = FALSE]
                deleted_reference <- reference_fit(deleted_data, fit$att_map, control)
                expect_equal(fit$jackknife$att[as.character(h), ], deleted_reference$att, tolerance = 1e-12,
                  info = paste(design_name, control, "ATT deletion", h))
                expect_equal(fit$jackknife$aggregate[as.character(h), ], deleted_reference$aggregate,
                  tolerance = 1e-12, info = paste(design_name, control, "aggregate deletion", h))
            }
        }
    }
})
test_that("weighted moments and direct moments fitting match the reference", {
    d <- make_test_data(TRUE)
    d$w <- 0.5 + (d$id%%7)/3
    moments <- collapse_rcs(d, "y", "time", "cohort", "region", weights = "w")
    fit <- fit_csdid_moments(moments, control_group = "not_yet", inference = "jackknife")
    expected <- reference_fit(d, fit$att_map, "not_yet", weights = "w")
    expect_equal(fit$att_gt, expected$att, tolerance = 1e-12)
    expect_equal(fit$aggregation$estimates, expected$aggregate, tolerance = 1e-12)
})
test_that("several treated cohorts with multiple clusters match the reference", {
    d <- expand.grid(region = 1:10, time = 1:6, obs = 1:4)
    d$cohort <- c(3, 3, 4, 4, 5, 5, 0, 0, 0, 0)[d$region]
    d$y <- d$region/7 + d$time/9 + d$obs/20 + 1.2 * (d$cohort > 0 & d$time >= d$cohort)
    d <- d[!(d$region == 8 & d$time == 2 & d$obs > 2), , drop = FALSE]
    for (control in c("never", "not_yet")) {
        fit <- fast_csdid(d, "y", "time", "cohort", "region", control_group = control, inference = "jackknife")
        expected <- reference_fit(d, fit$att_map, control)
        expect_equal(fit$att_gt, expected$att, tolerance = 1e-12)
        expect_equal(fit$aggregation$estimates, expected$aggregate, tolerance = 1e-12)
        expect_equal(sort(unique(fit$att_map$group)), c(3, 4, 5))
    }
})
test_that("numeric-like character and factor design labels are stable", {
    baseline <- make_test_data(TRUE)
    expected <- fast_csdid(baseline, "y", "time", "cohort", "region", inference = "none")
    character_data <- baseline
    character_data$region <- paste0("region-", character_data$region * 10 + 1)
    character_data$time <- as.character(character_data$time)
    character_data$cohort <- as.character(character_data$cohort)
    character_fit <- fast_csdid(character_data, "y", "time", "cohort", "region", never_treated = "0",
        inference = "none")
    factor_data <- character_data
    factor_data$region <- factor(factor_data$region)
    factor_data$time <- factor(factor_data$time)
    factor_data$cohort <- factor(factor_data$cohort)
    factor_fit <- fast_csdid(factor_data, "y", "time", "cohort", "region", never_treated = "0", inference = "none")
    expect_equal(character_fit$att_gt, expected$att_gt, tolerance = 1e-12)
    expect_equal(factor_fit$att_gt, expected$att_gt, tolerance = 1e-12)
    expect_equal(as.character(factor_fit$moments$regions), unique(as.character(factor_data$region)))
    integer_data <- baseline
    integer_data$region <- integer_data$region * 10L + 7L
    integer_fit <- fast_csdid(integer_data, "y", "time", "cohort", "region", inference = "none")
    expect_equal(integer_fit$att_gt, expected$att_gt, tolerance = 1e-12)
    expect_error(fast_csdid(transform(baseline, cohort = paste0("g", cohort)), "y", "time", "cohort",
        "region", inference = "none"), "numeric-like")
})
test_that("support-destroying deletions remain explicit and do not change the map", {
    d <- make_test_data(TRUE)
    d <- d[!(d$region == 2), , drop = FALSE]
    fit <- fast_csdid(d, "y", "time", "cohort", "region", inference = "jackknife")
    destroyed <- fit$att_map$group == 3
    preserved <- fit$att_map$group == 5
    expect_true(all(is.na(fit$jackknife$att["1", destroyed])))
    expect_true(all(fit$jackknife$att_status["1", destroyed] == 2L))
    expect_equal(fit$jackknife$att["1", preserved], fit$att_gt[preserved], tolerance = 0)
    expect_true(all(fit$jackknife$att_status["1", preserved] == 0L))
    expect_true(is.na(fit$jackknife$aggregate["1", "group:3"]))
    expect_true(is.finite(fit$jackknife$aggregate["1", "group:5"]))
    expect_equal(fit$jackknife$aggregate_status["1", "group:3"], 2L)
    expect_false(fit$jackknife$att_covariance$complete_by_stat[1L])
    expect_true(any(fit$jackknife$att_covariance$complete_by_stat[preserved]))
    expect_true(all(is.finite(diag(fit$jackknife$att_covariance$cv3)[preserved])))
})
test_that("deleting the sole never-treated control reports every failure", {
    d <- make_test_data(TRUE)
    d <- d[d$region %in% 1:5, , drop = FALSE]
    fit <- fast_csdid(d, "y", "time", "cohort", "region", inference = "jackknife")
    expect_true(all(is.na(fit$jackknife$att["5", ])))
    expect_true(all(fit$jackknife$att_status["5", ] == 2L))
    expect_true(all(is.na(fit$jackknife$aggregate["5", ])))
    expect_true(all(fit$jackknife$aggregate_status["5", ] == 2L))
    expect_equal(nrow(fit$diagnostics$jackknife_att_failures), nrow(fit$att_map))
    expect_equal(fit$jackknife$att["5", ], reference_att(d[d$region != 5, , drop = FALSE], fit$att_map),
        tolerance = 1e-12)
})
test_that("CV3 and CV3J covariance matrices equal their definitions", {
    d <- make_test_data(TRUE)
    fit <- fast_csdid(d, "y", "time", "cohort", "region", inference = "jackknife")
    deleted <- fit$jackknife$att
    full <- fit$att_gt
    H <- nrow(deleted)
    expected_cv3 <- (H - 1)/H * crossprod(sweep(deleted, 2L, full, `-`))
    expected_cv3j <- (H - 1)/H * crossprod(sweep(deleted, 2L, colMeans(deleted), `-`))
    expect_equal(fit$jackknife$att_covariance$cv3, expected_cv3, tolerance = 1e-14)
    expect_equal(fit$jackknife$att_covariance$cv3j, expected_cv3j, tolerance = 1e-14)
})
test_that("sparse dependency status is exact for both comparison rules", {
    d <- make_test_data(FALSE)
    never_fit <- fast_csdid(d, "y", "time", "cohort", "region", control_group = "never", inference = "jackknife")
    expect_true(all(never_fit$jackknife$att_status["1", never_fit$att_map$group == 5] == 0L))
    expect_true(all(never_fit$jackknife$att["1", never_fit$att_map$group == 5] == never_fit$att_gt[never_fit$att_map$group ==
        5]))
    not_yet_fit <- fast_csdid(d, "y", "time", "cohort", "region", control_group = "not_yet", inference = "jackknife")
    affected <- not_yet_fit$att_map$group == 5 | (not_yet_fit$att_map$group == 3 & not_yet_fit$att_map$post <
        5)
    expect_true(all(not_yet_fit$jackknife$att_status["3", affected] == 1L))
    expect_true(all(not_yet_fit$jackknife$att_status["3", !affected] == 0L))
    expect_equal(not_yet_fit$jackknife$att["3", !affected], not_yet_fit$att_gt[!affected], tolerance = 0)
})
