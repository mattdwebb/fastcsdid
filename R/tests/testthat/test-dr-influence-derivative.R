test_that("DR influence equals empirical cluster-weight derivatives", {
  # DRDID 1.3.0's repeated-cross-section influence formula has the same old
  # pre-period OR sign that this regression fixes. Its point estimates remain
  # valid comparators in the point-estimation tests. Here the independent
  # reference is a symmetric numerical derivative of the entire point fit.
  data <- read.csv(test_path("fixtures", "dr_derivative.csv"))
  for (control in c("never", "not_yet")) for (weighted in c(FALSE, TRUE)) {
    d <- data
    if (!weighted) d$w <- 1
    estimate <- function(x) fast_csdid(
      x, "y", "time", "g", "cluster", weights = "w", xformla = ~x+z,
      est_method = "dr", control_group = control, inference = "none"
    )
    fit <- estimate(d)
    influence <- cbind(
      influence_function(fit, type = "att", cluster_sum = TRUE),
      influence_function(fit, type = "aggregate", cluster_sum = TRUE)
    )
    derivative <- matrix(NA_real_, nrow(influence), ncol(influence))
    step <- 1e-5
    clusters <- sort(unique(d$cluster))
    for (j in seq_along(clusters)) {
      selector <- as.numeric(d$cluster == clusters[j])
      upper <- lower <- d
      upper$w <- d$w * (1 + step * selector)
      lower$w <- d$w * (1 - step * selector)
      plus <- estimate(upper); minus <- estimate(lower)
      derivative[j, ] <- (c(plus$att_gt, plus$aggregation$estimates) -
                          c(minus$att_gt, minus$aggregation$estimates))/(2*step)
    }
    expect_lt(max(abs(influence/nrow(d)-derivative)), 2e-7)
    expect_lt(max(abs(fit$vcov$analytical-crossprod(derivative))), 2e-8)
    k <- length(fit$att_gt)+1L
    expect_equal(sqrt(fit$vcov$analytical[k,k]),
                 sqrt(sum(derivative[,k]^2)), tolerance = 2e-8)
  }
})
