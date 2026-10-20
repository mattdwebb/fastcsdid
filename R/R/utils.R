.assert_scalar_string <- function(x, name) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    stop(name, " must be a single non-empty column name.", call. = FALSE)
  }
}

.require_columns <- function(data, columns) {
  missing <- setdiff(columns, names(data))
  if (length(missing)) {
    stop("Missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
}

.as_numeric_design <- function(x, name) {
  original <- x
  if (is.factor(x)) {
    x <- as.character(x)
  }
  out <- suppressWarnings(as.numeric(x))
  bad <- !is.na(original) & is.na(out)
  if (any(bad)) {
    stop(name, " must be numeric or contain numeric-like labels.", call. = FALSE)
  }
  out
}

.as_scalar_numeric <- function(x, name) {
  out <- .as_numeric_design(x, name)
  if (length(out) != 1L || is.na(out)) {
    stop(name, " must be one nonmissing numeric value.", call. = FALSE)
  }
  out
}

.safe_ratio <- function(numerator, denominator) {
  out <- rep(NA_real_, length(numerator))
  ok <- is.finite(denominator) & denominator > 0
  out[ok] <- numerator[ok] / denominator[ok]
  out
}

.att_names <- function(att_map) {
  paste0("ATT(", att_map$group, ",", att_map$post, ")")
}

.region_constant <- function(values, region_index, regions, variable_name) {
  pieces <- split(values, region_index)
  out <- vapply(pieces, function(z) {
    u <- unique(z)
    u <- u[!is.na(u)]
    if (length(u) != 1L) {
      stop(variable_name, " must be constant within region.", call. = FALSE)
    }
    u[[1L]]
  }, FUN.VALUE = values[[1L]])
  unname(out[as.character(seq_along(regions))])
}

.match_aggregation_types <- function(types) {
  choices <- c("simple", "group", "calendar", "dynamic")
  unique(match.arg(types, choices, several.ok = TRUE))
}

.normalize_assignment <- function(assignment, never_treated) {
  assignment <- .as_numeric_design(assignment, "assignment")
  never_treated <- .as_scalar_numeric(never_treated, "never_treated")
  assignment[assignment == never_treated] <- Inf
  assignment
}

.restore_assignment <- function(assignment, never_treated) {
  out <- assignment
  out[is.infinite(out)] <- never_treated
  out
}

.matrix_failures <- function(status, reason = NULL, unsupported_code = 2L,
                              row_name = "cluster", column_name = "statistic",
                              status_label = NULL) {
  index <- which(status == unsupported_code, arr.ind = TRUE)
  if (!nrow(index)) {
    out <- data.frame(character(), character(), character(), character(),
                      stringsAsFactors = FALSE)
    names(out) <- c(row_name, column_name, "status", "reason")
    return(out)
  }
  failure_reason <- if (is.null(reason)) {
    rep("unsupported", nrow(index))
  } else {
    reason[index]
  }
  failure_status <- if (is.null(status_label)) {
    rep("unsupported", nrow(index))
  } else {
    status_label[index]
  }
  out <- data.frame(
    row = rownames(status)[index[, "row"]],
    column = colnames(status)[index[, "col"]],
    status = failure_status,
    reason = failure_reason,
    stringsAsFactors = FALSE
  )
  names(out)[1:2] <- c(row_name, column_name)
  out
}
