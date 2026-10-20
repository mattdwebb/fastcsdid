test_that("post-treatment effects and aggregations match did 2.5.1", {
  skip_if_not_installed("did", minimum_version = "2.5.1")

  designs <- list(
    balanced = did_fixture(),
    unbalanced = did_fixture(unbalanced = TRUE),
    missing_cells = did_fixture(unbalanced = TRUE, missing_cells = TRUE),
    weighted = did_fixture(unbalanced = TRUE, weighted = TRUE)
  )
  for (design_name in names(designs)) {
    data <- designs[[design_name]]
    weighted <- identical(design_name, "weighted")
    for (control in c("never", "not_yet")) {
      fast <- fast_csdid(
        data, "y", "time", "cohort", "region",
        weights = if (weighted) "w" else NULL,
        control_group = control, inference = "none"
      )
      reference <- suppressMessages(suppressWarnings(
        fit_did_reference(data, control, weighted)
      ))
      compare_did_post_att(fast, reference)
      compare_did_aggregates(fast, reference)
    }
  }
})

test_that("first-period treatment is explicitly omitted as in did", {
  skip_if_not_installed("did", minimum_version = "2.5.1")
  data <- did_fixture(first_period_treated = TRUE)
  fast <- fast_csdid(data, "y", "time", "cohort", "region", inference = "none")
  reference <- suppressMessages(suppressWarnings(
    fit_did_reference(data, "never", FALSE)
  ))

  expect_equal(fast$diagnostics$omitted_groups$group, min(data$time))
  expect_identical(fast$diagnostics$omitted_groups$status, "missing_base_period")
  expect_false(min(data$time) %in% reference$group)
  compare_did_post_att(fast, reference)
})
