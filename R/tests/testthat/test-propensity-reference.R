test_that("independent standardized IPW and locally efficient DR formulas match did", {
  skip_if_not_installed("did", minimum_version = "2.5.1")
  designs <- list(
    balanced = make_regression_fixture(),
    unbalanced_weighted = make_regression_fixture(
      unbalanced = TRUE, weighted = TRUE, missing_cells = TRUE
    )
  )
  formulas <- list(intercept = ~1,
                   numeric = ~x1 + x2 + female,
                   interaction = ~x1 * sector + female)
  for (design_name in names(designs)) {
    data <- designs[[design_name]]
    weighted <- design_name == "unbalanced_weighted"
    for (method in c("ipw", "dr")) {
      for (formula_name in names(formulas)) {
        if (weighted && method == "dr" && formula_name == "interaction") next
        for (control in c("never", "not_yet")) {
          expected <- reference_nonlinear_fit(
            data, formulas[[formula_name]], method, control,
            weights = if (weighted) "w" else NULL
          )
          external <- suppressMessages(suppressWarnings(
            did_nonlinear_reference(
              data, formulas[[formula_name]], method, control, weighted,
              compute_inffunc = TRUE
            )
          ))
          external_att <- data.frame(
            group = external$group, post = external$t, att = external$att
          )
          external_att <- external_att[
            external_att$post >= external_att$group, , drop = FALSE
          ]
          joined <- merge(
            data.frame(group = expected$att_map$group,
                       post = expected$att_map$post,
                       att = unname(expected$att)),
            external_att, by = c("group", "post"), suffixes = c("_ref", "_did")
          )
          expect_equal(
            joined$att_ref, joined$att_did, tolerance = 1e-9,
            info = paste(design_name, method, formula_name, control, "ATT")
          )
          expect_equal(
            expected$aggregate, reference_did_aggregates(external),
            tolerance = 1e-9,
            info = paste(design_name, method, formula_name, control, "aggregate")
          )
        }
      }
    }
  }
})

test_that("reference propensity coefficients and probabilities match DRDID fits", {
  skip_if_not_installed("DRDID", minimum_version = "1.3.0")
  data <- make_regression_fixture(unbalanced = TRUE, weighted = TRUE)
  prepared <- reference_nonlinear_prepare(data, ~x1 + x2 + female, "w")
  cell <- reference_reg_map(data)[1L, ]
  sample <- reference_att_sample(prepared, cell, "not_yet")
  reference <- reference_logit(sample$X, sample$D, sample$w)
  external <- suppressWarnings(fastglm::fastglmPure(
    x = sample$X, y = sample$D, family = stats::binomial(),
    weights = sample$w / mean(sample$w), method = 3L,
    tol = 1e-8, maxit = 100L
  ))
  expect_equal(reference$coefficients, unname(external$coefficients),
               tolerance = 1e-8)
  expect_equal(reference$probabilities, unname(external$fitted.values),
               tolerance = 1e-9)
})
