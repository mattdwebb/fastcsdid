.finite_model_frame_rows <- function(model_frame) {
  ok <- stats::complete.cases(model_frame)
  for (column in model_frame) {
    if (is.matrix(column) && is.numeric(column)) {
      ok <- ok & apply(is.finite(column), 1L, all)
    } else if (is.numeric(column) || is.integer(column)) {
      ok <- ok & is.finite(column)
    }
  }
  ok
}

.prepare_regression_data <- function(data, outcome, time, cohort, cluster,
                                     weights, xformla) {
  model_frame_start <- proc.time()[["elapsed"]]
  if (!inherits(xformla, "formula") || length(xformla) != 2L) {
    stop("xformla must be a one-sided R formula such as ~ age + education.",
         call. = FALSE)
  }
  formula_variables <- all.vars(xformla)
  .require_columns(data, formula_variables)

  y <- .as_numeric_design(data[[outcome]], "outcome")
  period <- .as_numeric_design(data[[time]], "time")
  cohort_value <- .as_numeric_design(data[[cohort]], "cohort")
  weight <- if (is.null(weights)) rep.int(1, nrow(data)) else
    .as_numeric_design(data[[weights]], "weights")
  core_ok <- is.finite(y) & is.finite(period) & !is.na(cohort_value) &
    !is.na(data[[cluster]]) & is.finite(weight)

  initial_model_frame <- stats::model.frame(
    xformla, data = data, na.action = stats::na.pass,
    drop.unused.levels = FALSE
  )
  model_ok <- .finite_model_frame_rows(initial_model_frame)
  keep <- core_ok & model_ok
  if (!any(keep)) {
    stop("No complete finite observations remain after applying xformla.",
         call. = FALSE)
  }
  if (!all(keep)) {
    warning("dropped ", sum(!keep),
            " rows due to missing or non-finite outcome, design, weight, or covariate values",
            call. = FALSE)
  }
  clean_data <- droplevels(data[keep, , drop = FALSE])
  model_frame <- droplevels(initial_model_frame[keep, , drop = FALSE])
  model_terms <- attr(initial_model_frame, "terms")
  model_frame_time <- proc.time()[["elapsed"]] - model_frame_start
  model_matrix_start <- proc.time()[["elapsed"]]
  model_matrix <- stats::model.matrix(model_terms, data = model_frame)
  model_matrix_time <- proc.time()[["elapsed"]] - model_matrix_start
  if (!ncol(model_matrix)) {
    stop("xformla must produce at least one model-matrix column.", call. = FALSE)
  }
  if (any(!is.finite(model_matrix))) {
    stop("xformla produced non-finite model-matrix values.", call. = FALSE)
  }

  list(
    data = clean_data,
    X = unname(model_matrix),
    y = .as_numeric_design(clean_data[[outcome]], "outcome"),
    weight = if (is.null(weights)) rep.int(1, nrow(clean_data)) else
      .as_numeric_design(clean_data[[weights]], "weights"),
    formula = xformla,
    terms = model_terms,
    column_names = colnames(model_matrix),
    contrasts = attr(model_matrix, "contrasts"),
    assign = attr(model_matrix, "assign"),
    original_rows = which(keep),
    dropped_rows = which(!keep),
    timing = c(model_frame = model_frame_time,
               model_matrix = model_matrix_time)
  )
}

.finalize_regression_data <- function(regression_data, moments, time, cluster) {
  regression_data$region_index <- match(
    regression_data$data[[cluster]], moments$regions
  )
  regression_data$time_index <- match(
    .as_numeric_design(regression_data$data[[time]], "time"), moments$times
  )
  regression_data$assignment <-
    moments$assignment[regression_data$region_index]
  regression_data$data <- NULL
  regression_data
}
