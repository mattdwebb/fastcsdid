#' Collapse repeated-cross-section microdata
#'
#' Collapses individual observations to region-by-period sufficient statistics.
#' Treatment cohort must be constant within region. Never-treated regions are
#' represented internally by `Inf`.
#'
#' @param data A data frame.
#' @param outcome,time,cohort,cluster Column names supplied as strings.
#' @param weights Optional nonnegative sampling-weight column name. Weighted
#'   cell means are `sum(weights * outcome) / sum(weights)`. The same fixed
#'   weights determine treated-cohort aggregation mass, delete-one-cluster
#'   estimates, and analytical inference.
#' @param never_treated Value identifying never-treated regions.
#' @param xformla Optional one-sided covariate formula. When supplied, the model
#'   frame and model matrix are constructed once and cached for
#'   `est_method = "reg"`.
#' @return An object of class `fastrcsdid_moments`.
#' @export
collapse_rcs <- function(data, outcome, time, cohort, cluster,
                         weights = NULL, never_treated = 0, xformla = NULL) {
  if (!is.data.frame(data)) {
    stop("data must be a data frame.", call. = FALSE)
  }
  for (x in list(outcome = outcome, time = time, cohort = cohort, cluster = cluster)) {
    .assert_scalar_string(x, deparse(substitute(x)))
  }
  columns <- c(outcome, time, cohort, cluster)
  if (!is.null(weights)) {
    .assert_scalar_string(weights, "weights")
    columns <- c(columns, weights)
  }
  .require_columns(data, unique(columns))

  regression_data <- NULL
  if (!is.null(xformla)) {
    regression_data <- .prepare_regression_data(
      data, outcome, time, cohort, cluster, weights, xformla
    )
    data <- regression_data$data
  }

  y <- .as_numeric_design(data[[outcome]], "outcome")
  period <- .as_numeric_design(data[[time]], "time")
  cohort_value <- .as_numeric_design(data[[cohort]], "cohort")
  region_value <- data[[cluster]]
  w <- if (is.null(weights)) {
    rep.int(1, nrow(data))
  } else {
    .as_numeric_design(data[[weights]], "weights")
  }

  bad <- !is.finite(y) | !is.finite(period) | is.na(cohort_value) |
    is.na(region_value) | !is.finite(w) | w < 0
  if (any(bad)) {
    stop("Outcome, time, cohort, cluster, and weights must be observed and finite; weights must be nonnegative.", call. = FALSE)
  }
  if (!nrow(data)) {
    stop("data has no observations.", call. = FALSE)
  }

  regions <- unique(region_value)
  times <- sort(unique(period))
  region_index <- match(region_value, regions)
  time_index <- match(period, times)
  H <- length(regions)
  TT <- length(times)

  assignment_raw <- vapply(seq_len(H), function(h) {
    u <- unique(cohort_value[region_index == h])
    if (length(u) != 1L) {
      stop("cohort must be constant within cluster/region.", call. = FALSE)
    }
    u[[1L]]
  }, numeric(1))
  never_treated_numeric <- .as_scalar_numeric(never_treated, "never_treated")
  assignment <- .normalize_assignment(assignment_raw, never_treated_numeric)

  cell_index <- region_index + H * (time_index - 1L)
  n_obs <- matrix(tabulate(cell_index, nbins = H * TT), nrow = H, ncol = TT)
  sum_by_cell <- function(z) {
    out <- numeric(H * TT)
    collapsed <- rowsum(z, group = cell_index, reorder = FALSE)
    out[as.integer(rownames(collapsed))] <- collapsed[, 1L]
    matrix(out, nrow = H, ncol = TT)
  }
  weight_sum <- sum_by_cell(w)
  outcome_sum <- sum_by_cell(w * y)
  outcome_sq_sum <- sum_by_cell(w * y * y)

  dimnames(n_obs) <- dimnames(weight_sum) <- dimnames(outcome_sum) <-
    dimnames(outcome_sq_sum) <- list(as.character(regions), as.character(times))

  out <- list(
    regions = regions,
    times = times,
    assignment = assignment,
    assignment_raw = assignment_raw,
    never_treated = never_treated_numeric,
    n_obs = n_obs,
    weight_sum = weight_sum,
    outcome_sum = outcome_sum,
    outcome_sq_sum = outcome_sq_sum,
    region_total_weight = rowSums(weight_sum),
    N = nrow(data),
    total_weight = sum(w),
    weighting = list(
      weighted = !is.null(weights),
      column = weights,
      type = "sampling",
      cell_mean = "sum(weight * outcome) / sum(weight)",
      aggregation_mass = "sum of sampling weights in the treated cohort",
      zero_weight_observations = "retained in n_obs but excluded from weighted moments"
    ),
    call_columns = list(
      outcome = outcome,
      time = time,
      cohort = cohort,
      cluster = cluster,
      weights = weights
    )
  )
  class(out) <- "fastrcsdid_moments"
  if (!is.null(regression_data)) {
    out$regression_data <- .finalize_regression_data(
      regression_data, out, time, cluster
    )
  }
  out
}
