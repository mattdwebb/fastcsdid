test_that("complete-support covariance variants agree and are PSD", {
  deleted <- cbind(
    a = c(1.1, 0.9, 1.2, 0.8, 1.0),
    b = c(2.2, 1.8, 2.1, 1.9, 2.0),
    c = c(-0.8, -1.1, -0.9, -1.2, -1.0)
  )
  full <- c(1, 2, -1)
  covariance <- jackknife_covariance(deleted, full)
  expect_equal(covariance$cv3, covariance$common_cv3, tolerance = 0)
  # BLAS crossprod and scalar sum use different accumulation orders on some
  # platforms. Test agreement at double precision rather than bit identity.
  roundoff <- 64 * .Machine$double.eps
  expect_equal(covariance$cv3, covariance$pairwise_cv3, tolerance = roundoff)
  expect_equal(covariance$cv3j, covariance$common_cv3j, tolerance = 0)
  expect_equal(covariance$cv3j, covariance$pairwise_cv3j, tolerance = roundoff)
  expect_true(all(covariance$n_valid_cov == 5L))
  expect_true(all(covariance$valid_deletions))
  expect_true(all(covariance$positive_semidefinite))
})

test_that("pairwise covariance records each common deletion set", {
  deleted <- cbind(
    a = c(1.1, 0.9, 1.2, 0.8, 1.0),
    b = c(2.2, 1.8, 2.1, 1.9, 2.0),
    c = c(-0.8, -1.1, -0.9, -1.2, -1.0)
  )
  full <- c(1, 2, -1)
  deleted[1, "b"] <- NA_real_
  deleted[2, "c"] <- NA_real_
  covariance <- jackknife_covariance(deleted, full)

  expected_counts <- matrix(
    c(5, 4, 4, 4, 4, 3, 4, 3, 4), 3, 3,
    dimnames = list(colnames(deleted), colnames(deleted))
  )
  expect_equal(covariance$n_valid_cov, expected_counts)
  expect_equal(covariance$n_valid_by_stat, c(a = 5, b = 4, c = 4))
  expect_equal(covariance$common_repetitions, 3L)
  expect_true(is.na(covariance$cv3["b", "b"]))
  expect_true(is.finite(covariance$pairwise_cv3["b", "b"]))

  b_rows <- is.finite(deleted[, "b"])
  G <- sum(b_rows)
  expected_b <- (G - 1) / G * sum((deleted[b_rows, "b"] - full[2])^2)
  expect_equal(covariance$pairwise_cv3["b", "b"], expected_b,
               tolerance = 1e-15)

  bc_rows <- is.finite(deleted[, "b"]) & is.finite(deleted[, "c"])
  Gbc <- sum(bc_rows)
  expected_bc <- (Gbc - 1) / Gbc * sum(
    (deleted[bc_rows, "b"] - full[2]) *
      (deleted[bc_rows, "c"] - full[3])
  )
  expect_equal(covariance$pairwise_cv3["b", "c"], expected_bc,
               tolerance = 1e-15)
  expect_match(covariance$covariance_type["pairwise_cv3"], "pairwise")
  expect_type(covariance$positive_semidefinite, "logical")
})

test_that("one missing deletion retains its scalar pairwise variance", {
  deleted <- cbind(
    a = c(1.1, 0.9, 1.2, 0.8, 1.0),
    b = c(2.2, NA, 2.1, 1.9, 2.0)
  )
  covariance <- jackknife_covariance(deleted, c(1, 2))
  expect_equal(covariance$n_valid_cov["b", "b"], 4L)
  expect_equal(covariance$n_valid_cov["a", "b"], 4L)
  expect_true(is.finite(covariance$pairwise_se_cv3["b"]))
  expect_true(is.na(covariance$se_cv3["b"]))
  expect_equal(covariance$common_repetitions, 4L)
})
