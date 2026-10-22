#' Compile the post-treatment ATT map
#'
#' For each treated cohort, uses the latest observed period before treatment as
#' the base period and includes every observed post-treatment period.
#'
#' @param moments A `fastrcsdid_moments` object.
#' @param anticipation Anticipation length. Only zero is supported.
#' @return A data frame with one row per valid group-time effect.
#' @export
compile_att_map <- function(moments, anticipation = 0L) {
  if (!inherits(moments, "fastrcsdid_moments")) {
    stop("moments must be returned by collapse_rcs().", call. = FALSE)
  }
  anticipation <- as.integer(anticipation)
  if (length(anticipation) != 1L || is.na(anticipation) || anticipation != 0L) {
    stop("Only anticipation = 0 is supported.", call. = FALSE)
  }

  treated <- sort(unique(moments$assignment[is.finite(moments$assignment)]))
  omitted <- vapply(treated, function(g) {
    pre <- moments$times[moments$times < g]
    post <- moments$times[moments$times >= g]
    !length(pre) || !length(post)
  }, logical(1))
  rows <- lapply(treated, function(g) {
    pre <- moments$times[moments$times < g]
    post <- moments$times[moments$times >= g]
    if (!length(pre) || !length(post)) {
      return(NULL)
    }
    base <- max(pre)
    data.frame(
      group = rep(g, length(post)),
      base = rep(base, length(post)),
      post = post,
      event = post - g,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  if (is.null(out) || !nrow(out)) {
    stop("No valid post-treatment ATT(g,t) comparisons were found.", call. = FALSE)
  }
  rownames(out) <- NULL
  out$att_id <- seq_len(nrow(out))
  out <- out[, c("att_id", "group", "base", "post", "event")]
  attr(out, "omitted_groups") <- if (any(omitted)) {
    data.frame(
      group = treated[omitted],
      status = vapply(treated[omitted], function(g) {
        if (!any(moments$times < g)) "missing_base_period" else "missing_post_period"
      }, character(1)),
      reason = "no observed pre-treatment or post-treatment period",
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(group = numeric(), status = character(), reason = character(),
               stringsAsFactors = FALSE)
  }
  out
}
