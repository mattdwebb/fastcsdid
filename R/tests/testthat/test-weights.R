replicate_integer_weights <- function(data, weight) {
    index <- rep(seq_len(nrow(data)), times = weight)
    out <- data[index, setdiff(names(data), "w"), drop = FALSE]
    out$id <- seq_len(nrow(out))
    rownames(out) <- NULL
    out
}
test_that("moments carry explicit sampling-weight semantics", {
    data <- make_test_data(TRUE)
    unweighted <- collapse_rcs(data, "y", "time", "cohort", "region")
    data$w <- 0.5 + (data$id%%7L)/3
    weighted <- collapse_rcs(data, "y", "time", "cohort", "region", weights = "w")
    expect_identical(unweighted$weighting$weighted, FALSE)
    expect_null(unweighted$weighting$column)
    expect_identical(weighted$weighting$weighted, TRUE)
    expect_identical(weighted$weighting$column, "w")
    expect_identical(weighted$weighting$type, "sampling")
    expect_match(weighted$weighting$cell_mean, "sum\\(weight \\* outcome\\)")
    selector <- data$region == data$region[1L] & data$time == data$time[1L]
    expect_equal(weighted$outcome_sq_sum[1L, 1L], sum(data$w[selector] * data$y[selector]^2), tolerance = 1e-12)
    fit <- fit_csdid_moments(weighted, inference = "none")
    expect_identical(fit$settings$weighting, weighted$weighting)
})
test_that("sampling weights reject missing, nonfinite, and negative values", {
    data <- make_test_data(FALSE)
    data$w <- 1
    for (bad in list(NA_real_, Inf, -0.01)) {
        invalid <- data
        invalid$w[1L] <- bad
        expect_error(collapse_rcs(invalid, "y", "time", "cohort", "region", weights = "w"), "weights must be nonnegative")
    }
    expect_error(collapse_rcs(data, "y", "time", "cohort", "region", weights = "missing"), "Missing required columns")
})
test_that("integer sampling weights equal literal row replication", {
    data <- make_test_data(TRUE)
    data$w <- 1L + data$id%%4L
    replicated <- replicate_integer_weights(data, data$w)
    for (control in c("never", "not_yet")) {
        weighted <- fast_csdid(data, "y", "time", "cohort", "region", weights = "w", control_group = control,
            inference = "jackknife")
        expanded <- fast_csdid(replicated, "y", "time", "cohort", "region", control_group = control,
            inference = "jackknife")
        expect_equal(weighted$att_gt, expanded$att_gt, tolerance = 1e-12)
        expect_equal(weighted$aggregation$estimates, expanded$aggregation$estimates, tolerance = 1e-12)
        expect_equal(weighted$jackknife$att, expanded$jackknife$att, tolerance = 1e-12)
        expect_equal(weighted$jackknife$aggregate, expanded$jackknife$aggregate, tolerance = 1e-12)
    }
})
test_that("zero-weight observations equal removal and zero-mass cells are unsupported", {
    data <- make_test_data(TRUE)
    data$w <- 0.2 + (data$id%%13L)/4
    zero_rows <- data$region == 8L & data$time == 2L & data$id%%2L == 0L
    data$w[zero_rows] <- 0
    weighted <- fast_csdid(data, "y", "time", "cohort", "region", weights = "w", control_group = "not_yet",
        inference = "jackknife")
    removed <- fast_csdid(data[!zero_rows, ], "y", "time", "cohort", "region", weights = "w", control_group = "not_yet",
        inference = "jackknife")
    expect_equal(weighted$att_gt, removed$att_gt, tolerance = 1e-12)
    expect_equal(weighted$aggregation$estimates, removed$aggregation$estimates, tolerance = 1e-12)
    expect_equal(weighted$jackknife$att, removed$jackknife$att, tolerance = 1e-12)
    all_zero <- data
    all_zero$w[all_zero$cohort == 3 & all_zero$time == 2] <- 0
    unsupported <- fast_csdid(all_zero, "y", "time", "cohort", "region", weights = "w", inference = "none")
    affected <- unsupported$att_map$group == 3
    expect_true(all(!unsupported$att_support$supported[affected]))
    expect_true(all(is.na(unsupported$att_gt[affected])))
    expect_true(all(unsupported$att_support$treated_base_weight[affected] == 0))
    expect_true(all(unsupported$att_support$status[affected] == "zero_weight_treated_cell"))
    expect_true(all(grepl("treated_base", unsupported$att_support$reason[affected])))
})
test_that("fractional and highly unequal weighted deletions match physical deletion", {
    data <- make_test_data(TRUE)
    data$w <- 0.001 + (data$id%%17L)/13
    data$w[data$region == 2L & data$time == 6L] <- 1000
    for (control in c("never", "not_yet")) {
        fit <- fast_csdid(data, "y", "time", "cohort", "region", weights = "w", control_group = control,
            inference = "jackknife")
        expected <- reference_fit(data, fit$att_map, control, weights = "w")
        expect_equal(fit$att_gt, expected$att, tolerance = 1e-11)
        expect_equal(fit$aggregation$estimates, expected$aggregate, tolerance = 1e-11)
        for (h in unique(data$region)) {
            deleted <- data[data$region != h, , drop = FALSE]
            delete_expected <- reference_fit(deleted, fit$att_map, control, weights = "w")
            expect_equal(fit$jackknife$att[as.character(h), ], delete_expected$att, tolerance = 1e-10)
            expect_equal(fit$jackknife$aggregate[as.character(h), ], delete_expected$aggregate, tolerance = 1e-10)
        }
    }
})
test_that("weighted sparse deletions leave unrelated ATT entries bitwise unchanged", {
    data <- make_test_data(FALSE)
    data$w <- 0.1 + (data$id%%13L)/7
    never <- fast_csdid(data, "y", "time", "cohort", "region", weights = "w", control_group = "never",
        inference = "jackknife")
    unrelated <- never$att_map$group == 5
    expect_identical(unname(never$jackknife$att["1", unrelated]), unname(never$att_gt[unrelated]))
    expect_true(all(never$jackknife$att_status["1", unrelated] == 0L))
    not_yet <- fast_csdid(data, "y", "time", "cohort", "region", weights = "w", control_group = "not_yet",
        inference = "jackknife")
    affected <- not_yet$att_map$group == 5 | (not_yet$att_map$group == 3 & not_yet$att_map$post < 5)
    expect_identical(unname(not_yet$jackknife$att["3", !affected]), unname(not_yet$att_gt[!affected]))
    expect_true(all(not_yet$jackknife$att_status["3", !affected] == 0L))
})
test_that("weighted support failures distinguish cohort, treated-cell, and control failures", {
    data <- make_test_data(FALSE)
    data$w <- 1
    sole_cohort <- data[data$region != 2L, , drop = FALSE]
    cohort_fit <- fast_csdid(sole_cohort, "y", "time", "cohort", "region", weights = "w", inference = "jackknife")
    cohort_cells <- cohort_fit$att_map$group == 3
    expect_true(all(cohort_fit$jackknife$att_status_label["1", cohort_cells] == "cohort_removed"))
    treated_cell <- data
    treated_cell$w[treated_cell$region == 2L & treated_cell$time == 2L] <- 0
    treated_fit <- fast_csdid(treated_cell, "y", "time", "cohort", "region", weights = "w", inference = "jackknife")
    expect_true(all(treated_fit$jackknife$att_status_label["1", treated_fit$att_map$group == 3] == "zero_weight_treated_cell"))
    control_cell <- data
    control_cell$w[control_cell$region %in% 6:8 & control_cell$time == 2L] <- 0
    control_fit <- fast_csdid(control_cell, "y", "time", "cohort", "region", weights = "w", inference = "jackknife")
    expect_true(all(control_fit$jackknife$att_status_label["5", control_fit$att_map$base == 2] == "zero_weight_control_cell"))
})
