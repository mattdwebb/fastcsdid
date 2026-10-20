test_that("influence summaries reproduce the raw deletion matrices", {
    data <- make_test_data(TRUE)
    data$w <- 0.2 + (data$id%%9L)/5
    fit <- fast_csdid(data, "y", "time", "cohort", "region", weights = "w", control_group = "not_yet",
        inference = "jackknife")
    summary <- jackknife_influence_summary(fit)
    expect_equal(nrow(summary), length(fit$att_gt) + length(fit$aggregation$estimates))
    expect_setequal(summary$estimand, c(names(fit$att_gt), names(fit$aggregation$estimates)))
    expect_equal(summary, fit$diagnostics$jackknife_influence)
    label <- names(fit$att_gt)[1L]
    table <- jackknife_cluster_table(fit, label)
    deleted <- fit$jackknife$att[, label]
    valid <- is.finite(deleted)
    row <- summary[summary$estimand == label, ]
    expect_equal(table$delete_estimate, unname(deleted))
    expect_equal(table$change, unname(deleted - fit$att_gt[label]))
    expect_equal(row$min_delete, min(deleted[valid]))
    expect_equal(row$max_delete, max(deleted[valid]))
    expect_equal(row$max_abs_change, max(abs(deleted[valid] - fit$att_gt[label])))
    expected_cluster <- names(which.max(abs(deleted - fit$att_gt[label])))
    expect_identical(row$most_influential_cluster, expected_cluster)
})
test_that("near-zero estimates have deliberately undefined relative changes", {
    data <- make_test_data(FALSE)
    data$y <- data$region + 2 * data$time + data$id/10000
    fit <- fast_csdid(data, "y", "time", "cohort", "region", inference = "jackknife")
    label <- names(fit$att_gt)[1L]
    fit$att_gt[label] <- 0
    table <- jackknife_cluster_table(fit, label)
    summary <- jackknife_influence_summary(fit)
    expect_true(all(is.na(table$relative_change)))
    expect_true(is.na(summary$largest_relative_change[summary$estimand == label]))
})
test_that("sign and fixed-SE significance flips are flagged", {
    fit <- fast_csdid(make_test_data(TRUE), "y", "time", "cohort", "region", inference = "jackknife")
    label <- names(fit$att_gt)[1L]
    position <- match(label, names(fit$att_gt))
    fit$att_gt[position] <- 0.2
    fit$jackknife$att[, position] <- 0.21
    fit$jackknife$att[1L, position] <- -0.05
    fit$jackknife$att_covariance$pairwise_se_cv3[position] <- 0.1
    table <- jackknife_cluster_table(fit, label)
    expect_true(table$sign_flip[1L])
    expect_true(table$significance_flip[1L])
    expect_false(any(table$significance_flip[-1L]))
    summary <- jackknife_influence_summary(fit)
    row <- summary[summary$estimand == label, ]
    expect_true(row$sign_flip)
    expect_true(row$significance_flip)
    expect_identical(row$most_influential_cluster, table$cluster[1L])
})
test_that("diagnostic helpers reject missing jackknife results and bad labels", {
    fit <- fast_csdid(make_test_data(FALSE), "y", "time", "cohort", "region", inference = "none")
    expect_error(jackknife_influence_summary(fit), "does not contain jackknife")
    fit <- fast_csdid(make_test_data(FALSE), "y", "time", "cohort", "region", inference = "jackknife")
    expect_error(jackknife_cluster_table(fit, "not-an-estimand"), "unknown estimand")
    expect_error(jackknife_influence_summary(fit, relative_threshold = -1), "nonnegative")
})
