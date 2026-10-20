test_that("zero-mass unsupported cells do not invalidate an aggregate", {
  att_map <- data.frame(
    group = c(4, 4, 6),
    base = c(3, 3, 5),
    post = c(4, 5, 6),
    event = c(0, 1, 0)
  )
  att <- c(NA_real_, NA_real_, 2)
  group_size <- c(`4` = 0, `6` = 10)

  result <- aggregate_att(
    att, att_map, group_size,
    types = c("simple", "group", "calendar", "dynamic")
  )

  expect_equal(unname(result$estimates["simple"]), 2)
  expect_equal(unname(result$weights["simple", ]), c(0, 0, 1))
  expect_true(is.na(result$estimates["group:4"]))
  expect_equal(unname(result$estimates["group:6"]), 2)
  expect_true(is.na(result$estimates["calendar:4"]))
  expect_true(is.na(result$estimates["event:1"]))
  expect_equal(unname(result$estimates["event:0"]), 2)
})

test_that("positive-mass unsupported cells invalidate every affected aggregate", {
  att_map <- data.frame(
    group = c(4, 4, 6, 6),
    base = c(3, 3, 5, 5),
    post = c(4, 5, 6, 7),
    event = c(0, 1, 0, 1)
  )
  result <- aggregate_att(
    c(NA_real_, 1, 2, 3), att_map, c(`4` = 10, `6` = 20),
    types = c("simple", "group", "calendar", "dynamic")
  )

  expect_true(is.na(result$estimates["simple"]))
  expect_true(is.na(result$estimates["group:4"]))
  expect_true(is.na(result$estimates["calendar:4"]))
  expect_true(is.na(result$estimates["event:0"]))
  expect_equal(unname(result$estimates["group:6"]), 2.5)
  expect_equal(unname(result$estimates["event:1"]), 7 / 3)
})

test_that("multiple zero-mass unsupported cells are removed and renormalized", {
  att_map <- data.frame(
    group = c(4, 4, 4, 6, 6),
    base = c(3, 3, 3, 5, 5),
    post = c(4, 5, 6, 6, 7),
    event = c(0, 1, 2, 0, 1)
  )
  result <- aggregate_att(
    c(NA_real_, NA_real_, NA_real_, 2, 4), att_map,
    c(`4` = 0, `6` = 12),
    types = c("simple", "calendar", "dynamic")
  )

  expect_equal(unname(result$estimates["simple"]), 3)
  expect_equal(unname(result$weights["simple", ]), c(0, 0, 0, 0.5, 0.5))
  expect_equal(unname(result$estimates["calendar:6"]), 2)
  expect_equal(unname(result$estimates["event:0"]), 2)
  expect_equal(unname(result$estimates["event:1"]), 4)
  expect_true(is.na(result$estimates["event:2"]))
})

test_that("all zero-mass cells produce an explicitly undefined aggregate", {
  att_map <- data.frame(
    group = c(4, 4), base = c(3, 3), post = c(4, 5), event = c(0, 1)
  )
  result <- aggregate_att(
    c(NA_real_, NA_real_), att_map, c(`4` = 0),
    types = c("simple", "group", "calendar", "dynamic")
  )

  expect_true(all(is.na(result$estimates)))
  expect_true(all(is.na(result$weights)))
})

test_that("deleting a singleton treated cohort reaggregates surviving cohorts", {
  d <- make_test_data(TRUE)
  d <- d[!(d$region == 2), , drop = FALSE] # cohort 3 now has one region
  fit <- fast_csdid(d, "y", "time", "cohort", "region",
                    aggregation = "simple", inference = "jackknife",
                    jackknife_support = "changing")

  expect_true(all(is.na(fit$jackknife$att["1", fit$att_map$group == 3])))
  expect_true(is.finite(fit$jackknife$aggregate["1", "simple"]))
  expect_equal(fit$jackknife$aggregate_status["1", "simple"], 1L)

  surviving <- fit$jackknife$att["1", fit$att_map$group == 5]
  expect_equal(
    unname(fit$jackknife$aggregate["1", "simple"]),
    mean(surviving),
    tolerance = 1e-12
  )
})

singleton_fixture <- function(weighted = FALSE, one_cohort = FALSE) {
  d <- did_fixture(unbalanced = TRUE, weighted = weighted)
  d <- d[d$region %in% c(1L, 3L, 7L, 8L, 9L, 10L), , drop = FALSE]
  d$cohort <- 0
  d$cohort[d$region == 1L] <- 2006
  if (!one_cohort) d$cohort[d$region == 3L] <- 2010
  d
}

expect_singleton_reaggregation <- function(fit, deleted_cluster, removed_group,
                                           surviving_group) {
  removed <- fit$att_map$group == removed_group
  surviving <- fit$att_map$group == surviving_group
  row <- as.character(deleted_cluster)

  expect_true(all(is.na(fit$jackknife$att[row, removed])))
  expect_true(all(fit$jackknife$att_status_label[row, removed] %in% c(
    "cohort_removed", "zero_propensity_normalizer", "rank_changed"
  )))
  expect_true(is.na(fit$jackknife$aggregate[row, paste0("group:", removed_group)]))
  expect_equal(fit$jackknife$aggregate_status_label[
    row, paste0("group:", removed_group)
  ], "undefined_aggregation")
  expect_true(is.finite(fit$jackknife$aggregate[row, "simple"]))

  deleted_att <- fit$jackknife$att[row, surviving]
  expect_equal(
    unname(fit$jackknife$aggregate[row, "simple"]),
    mean(deleted_att), tolerance = 1e-10
  )
  expect_equal(
    unname(fit$jackknife$aggregate[row, paste0("group:", surviving_group)]),
    mean(deleted_att), tolerance = 1e-10
  )

  # Calendar and dynamic rows that retain supported cells remain defined.
  calendar <- fit$aggregation$metadata$label[
    fit$aggregation$metadata$type == "calendar"
  ]
  dynamic <- fit$aggregation$metadata$label[
    fit$aggregation$metadata$type == "dynamic"
  ]
  expect_true(any(is.finite(fit$jackknife$aggregate[row, calendar])))
  expect_true(any(is.finite(fit$jackknife$aggregate[row, dynamic])))
}

test_that("early and late singleton deletions work for all aggregation types", {
  fit <- fast_csdid(
    singleton_fixture(), "y", "time", "cohort", "region",
    aggregation = c("simple", "group", "calendar", "dynamic"),
    inference = "jackknife", est_method = "none", jackknife_support = "changing"
  )

  expect_singleton_reaggregation(fit, 1L, 2006, 2010)
  expect_singleton_reaggregation(fit, 3L, 2010, 2006)
  expect_true(all(is.finite(fit$att_gt)))
  expect_true(all(is.finite(fit$aggregation$estimates)))
})

test_that("weighted no-covariate, regression, IPW, and DR deletions renormalize", {
  d <- singleton_fixture(weighted = TRUE)
  d$x <- sin(d$id / 5) + d$region / 11 + match(d$time, unique(d$time)) / 7
  specifications <- list(
    none = list(xformla = NULL, jackknife_method = "linearized"),
    reg = list(xformla = ~x, jackknife_method = "linearized"),
    ipw = list(xformla = ~x, jackknife_method = "refit"),
    dr = list(xformla = ~x, jackknife_method = "refit")
  )

  for (method in names(specifications)) {
    specification <- specifications[[method]]
    fit <- fast_csdid(
      d, "y", "time", "cohort", "region", weights = "w",
      xformla = specification$xformla, est_method = method,
      aggregation = c("simple", "group", "calendar", "dynamic"),
      inference = "jackknife",
      jackknife_method = specification$jackknife_method, jackknife_support = "changing"
    )
    expect_singleton_reaggregation(fit, 1L, 2006, 2010)
    expect_singleton_reaggregation(fit, 3L, 2010, 2006)
  }
})

test_that("deleting the only treated cluster makes aggregates explicitly undefined", {
  fit <- fast_csdid(
    singleton_fixture(one_cohort = TRUE), "y", "time", "cohort", "region",
    aggregation = c("simple", "group", "calendar", "dynamic"),
    inference = "jackknife"
  )
  deleted <- fit$jackknife$aggregate["1", ]

  expect_true(all(is.na(deleted)))
  expect_true(all(fit$jackknife$aggregate_status["1", ] == 2L))
  expect_true(all(fit$jackknife$aggregate_status_label["1", ] ==
                    "undefined_aggregation"))
})
