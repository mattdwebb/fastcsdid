test_that("independent repeated-CS regression formula matches did 2.5.1", {
  skip_if_not_installed("did", minimum_version = "2.5.1")
  designs <- list(
    balanced = make_regression_fixture(),
    unbalanced = make_regression_fixture(unbalanced = TRUE),
    missing_cells = make_regression_fixture(unbalanced = TRUE, missing_cells = TRUE),
    weighted = make_regression_fixture(unbalanced = TRUE, weighted = TRUE)
  )
  formulas <- list(
    intercept = ~1,
    one_numeric = ~x1,
    multiple = ~x1 + x2 + female,
    factor_interaction = ~x1 * sector + female
  )

  for (design_name in names(designs)) {
    data <- designs[[design_name]]
    weighted <- identical(design_name, "weighted")
    for (formula_name in names(formulas)) {
      xformla <- formulas[[formula_name]]
      for (control in c("never", "not_yet")) {
        expected <- reference_reg_fit(
          data, xformla, control, weights = if (weighted) "w" else NULL
        )
        external <- suppressMessages(suppressWarnings(
          did_reg_reference(data, xformla, control, weighted)
        ))
        external_att <- data.frame(
          group = external$group, post = external$t, estimate = external$att
        )
        external_att <- external_att[external_att$post >= external_att$group, ]
        joined <- merge(
          data.frame(group = expected$att_map$group, post = expected$att_map$post,
                     reference = unname(expected$att)),
          external_att, by = c("group", "post"), all = TRUE
        )
        expect_false(anyNA(joined), info = paste(design_name, formula_name, control))
        expect_equal(joined$reference, joined$estimate, tolerance = 1e-10,
                     info = paste(design_name, formula_name, control, "ATT"))
      }
    }
  }
})

test_that("independent regression aggregations match did", {
  skip_if_not_installed("did", minimum_version = "2.5.1")
  data <- make_regression_fixture(unbalanced = TRUE, weighted = TRUE,
                                  missing_cells = TRUE)
  xformla <- ~x1 * sector + x2 + female
  for (control in c("never", "not_yet")) {
    expected <- reference_reg_fit(data, xformla, control, weights = "w")
    external <- suppressMessages(suppressWarnings(
      did_reg_reference(data, xformla, control, TRUE, compute_inffunc = TRUE)
    ))
    simple <- did::aggte(external, type = "simple", bstrap = FALSE, cband = FALSE)
    expect_equal(unname(expected$aggregate["simple"]), simple$overall.att,
                 tolerance = 1e-10)
    for (type in c("group", "calendar", "dynamic")) {
      result <- did::aggte(external, type = type, bstrap = FALSE, cband = FALSE)
      keep <- if (type == "dynamic") result$egt >= 0 else rep(TRUE, length(result$egt))
      prefix <- switch(type, group = "group:", calendar = "calendar:",
                       dynamic = "event:")
      labels <- paste0(prefix, result$egt[keep])
      expect_equal(unname(expected$aggregate[labels]), result$att.egt[keep],
                   tolerance = 1e-10, info = paste(control, type))
    }
  }
})
