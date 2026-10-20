test_that("sparse delete matrix equals brute-force deletion", {
  d <- make_test_data(unbalanced = TRUE)
  fit <- fast_csdid(d, "y", "time", "cohort", "region", inference = "jackknife")

  for (h in unique(d$region)) {
    brute <- fast_csdid(d[d$region != h, ], "y", "time", "cohort", "region",
                        inference = "none")
    expect_equal(
      unname(fit$jackknife$att[as.character(h), ]),
      unname(brute$att_gt),
      tolerance = 1e-12
    )
    expect_equal(
      unname(fit$jackknife$aggregate[as.character(h), ]),
      unname(brute$aggregation$estimates),
      tolerance = 1e-12
    )
  }
})

test_that("unaffected never-control ATT cells are unchanged", {
  d <- make_test_data()
  fit <- fast_csdid(d, "y", "time", "cohort", "region", inference = "jackknife")
  # Region 1 belongs to cohort 3, so cohort-5 ATT cells do not use it.
  late <- fit$att_map$group == 5
  expect_equal(
    unname(fit$jackknife$att["1", late]),
    unname(fit$att_gt[late]),
    tolerance = 0
  )
})
