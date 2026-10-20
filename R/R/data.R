#' Example repeated-cross-section data
#'
#' A small simulated repeated-cross-section data set for examples and smoke
#' tests. It contains 12 regions observed over 8 periods, with first treatment
#' in periods 4, 5, or 6 for treated regions and `0` for never-treated regions.
#'
#' @format A data frame with 960 rows and 8 variables:
#' \describe{
#'   \item{region}{Region or cluster identifier.}
#'   \item{time}{Time period.}
#'   \item{cohort}{First treated period, with `0` denoting never treated.}
#'   \item{unit}{Within-cell observation identifier.}
#'   \item{y}{Outcome.}
#'   \item{x1}{Continuous covariate.}
#'   \item{x2}{Binary covariate.}
#'   \item{weight}{Positive sampling weight.}
#' }
#' @examples
#' data(fastrcsdid_example)
#' fit <- fast_csdid(
#'   fastrcsdid_example,
#'   outcome = "y", time = "time", cohort = "cohort", cluster = "region",
#'   weights = "weight", inference = "none"
#' )
#' fit$aggregation$estimates
"fastrcsdid_example"
