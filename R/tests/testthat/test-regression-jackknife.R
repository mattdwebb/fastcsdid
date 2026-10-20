test_that("exact regression jackknife matches physical deletion and refitting", {
  designs <- list(
    unweighted = make_regression_fixture(unbalanced = TRUE, missing_cells = TRUE),
    weighted = make_regression_fixture(unbalanced = TRUE, weighted = TRUE,
                                       missing_cells = TRUE)
  )
  designs$weighted$w[designs$weighted$id %% 17L == 0L] <- 0
  designs$weighted$w[
    designs$weighted$region == 2L & designs$weighted$time == max(designs$weighted$time)
  ] <- 500
  formulas <- list(numeric = ~x1 + x2 + female,
                   interaction = ~x1 * sector + female)
  for (design_name in names(designs)) {
    data <- designs[[design_name]]
    weights <- if (design_name == "weighted") "w" else NULL
    for (formula_name in names(formulas)) {
      xformla <- formulas[[formula_name]]
      for (control in c("never", "not_yet")) {
        fit <- fast_csdid(
          data, "y", "time", "cohort", "region", weights = weights,
          xformla = xformla, est_method = "reg", control_group = control,
          inference = "jackknife"
        )
        for (cluster in unique(data$region)) {
          deleted_data <- data[data$region != cluster, , drop = FALSE]
          expected <- reference_reg_fit(
            deleted_data, xformla, control, weights = weights,
            att_map = fit$att_map
          )
          expect_equal(
            fit$jackknife$att[as.character(cluster), ], expected$att,
            tolerance = 1e-10,
            info = paste(design_name, formula_name, control, "ATT", cluster)
          )
          expect_equal(
            fit$jackknife$aggregate[as.character(cluster), ], expected$aggregate,
            tolerance = 1e-10,
            info = paste(design_name, formula_name, control, "aggregate", cluster)
          )
        }
      }
    }
  }
})

test_that("unaffected regression ATT entries are exact copies", {
  data <- make_regression_fixture(unbalanced = TRUE, weighted = TRUE)
  never <- fast_csdid(
    data, "y", "time", "cohort", "region", weights = "w",
    xformla = ~x1 + x2, est_method = "reg", control_group = "never",
    inference = "jackknife"
  )
  diagnostics <- never$jackknife$nuisance_diagnostics
  unaffected <- diagnostics$att_influence_source == "unaffected"
  full <- matrix(never$att_gt, nrow = nrow(unaffected), ncol = ncol(unaffected),
                 byrow = TRUE)
  expect_true(any(unaffected))
  expect_identical(unname(never$jackknife$att[unaffected]), unname(full[unaffected]))
  expect_true(all(never$jackknife$att_status[unaffected] == 0L))

  not_yet <- fast_csdid(
    data, "y", "time", "cohort", "region", weights = "w",
    xformla = ~x1 + x2, est_method = "reg", control_group = "not_yet",
    inference = "jackknife"
  )
  diagnostics <- not_yet$jackknife$nuisance_diagnostics
  unaffected <- diagnostics$att_influence_source == "unaffected"
  full <- matrix(not_yet$att_gt, nrow = nrow(unaffected), ncol = ncol(unaffected),
                 byrow = TRUE)
  expect_true(any(unaffected))
  expect_identical(unname(not_yet$jackknife$att[unaffected]),
                   unname(full[unaffected]))
})

test_that("regression deletion reports cohort removal and exact rank loss", {
  data <- make_regression_fixture(unbalanced = FALSE)
  sole_cohort <- data[data$region != 2L, , drop = FALSE]
  cohort_fit <- fast_csdid(
    sole_cohort, "y", "time", "cohort", "region", xformla = ~x1,
    est_method = "reg", inference = "jackknife"
  )
  cohort_cells <- cohort_fit$att_map$group == 2006
  expect_true(all(is.na(cohort_fit$jackknife$att["1", cohort_cells])))
  expect_true(all(cohort_fit$jackknife$att_status_label["1", cohort_cells] ==
                    "cohort_removed"))

  rank_data <- make_test_data(FALSE)
  rank_data <- rank_data[rank_data$region %in% 1:6, , drop = FALSE]
  rank_data$x_control_cluster <- as.numeric(rank_data$region == 6L)
  rank_fit <- fast_csdid(
    rank_data, "y", "time", "cohort", "region",
    xformla = ~x_control_cluster, est_method = "reg",
    inference = "jackknife"
  )
  expect_true(all(is.na(rank_fit$jackknife$att["6", ])))
  expect_true(all(rank_fit$jackknife$att_status_label["6", ] == "rank_changed"))
  expect_true(all(rank_fit$jackknife$nuisance_diagnostics$att_rank_changed["6", ]))
  expect_true(all(grepl("block:",
                        rank_fit$jackknife$nuisance_diagnostics$att_failed_block["6", ])))

  table <- jackknife_cluster_table(rank_fit, names(rank_fit$att_gt)[1L])
  expect_true(table$rank_changed[table$cluster == "6"])
  expect_match(table$failed_nuisance_block[table$cluster == "6"], "block:")
  summary <- jackknife_influence_summary(rank_fit)
  att_summary <- summary[summary$estimand_type == "att", ]
  expect_true(all(att_summary$rank_changed))
  expect_true(all(grepl("block:", att_summary$failed_nuisance_blocks)))
})

test_that("regression influence diagnostics expose coefficient and source changes", {
  data <- make_regression_fixture(unbalanced = TRUE, weighted = TRUE)
  fit <- fast_csdid(
    data, "y", "time", "cohort", "region", weights = "w",
    xformla = ~x1 + x2 + female, est_method = "reg",
    control_group = "not_yet", inference = "jackknife"
  )
  label <- names(fit$att_gt)[1L]
  table <- jackknife_cluster_table(fit, label)
  expect_true(any(is.finite(table$max_coefficient_change)))
  expect_true(any(table$influence_source == "control_regression"))
  expect_true(any(grepl("treated_moments", table$influence_source)))
  summary <- jackknife_influence_summary(fit)
  row <- summary[summary$estimand == label, ]
  expect_true(is.finite(row$max_coefficient_change))
  expect_match(row$influence_sources, "control_regression")
  expect_equal(fit$diagnostics$jackknife_influence,
               jackknife_influence_summary(fit))
})
