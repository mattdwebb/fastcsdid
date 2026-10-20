test_that("simple aggregation uses cohort-size weights per ATT cell", {
  d <- make_test_data(unbalanced = TRUE)
  fit <- fast_csdid(d, "y", "time", "cohort", "region", inference = "none")
  gs <- fit$aggregation$weights["simple", ]
  raw <- fit$moments$region_total_weight
  group_size <- vapply(fit$att_map$group, function(g) {
    sum(raw[fit$moments$assignment == g])
  }, numeric(1))
  expected_weight <- group_size / sum(group_size)
  expect_equal(unname(gs), unname(expected_weight), tolerance = 1e-12)
  expect_equal(unname(fit$aggregation$estimates["simple"]),
               sum(expected_weight * fit$att_gt), tolerance = 1e-12)
})
