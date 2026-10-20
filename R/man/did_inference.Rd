\name{did_design}
\alias{did_design}\alias{did_point}\alias{did_jackknife}
\title{No-covariate state-level inference adapter}
\description{Caches repeated-cross-section moments for native simple CSDID and mean-centered delete-one-state jackknife t inference.}
\usage{
did_design(data, outcome, cluster, time, timing)
did_point(design, assignment = NULL, engine = c("fast", "reference"))
did_jackknife(design, assignment = NULL, null = 0, engine = c("fast", "reference"))
}
\arguments{
\item{data}{Repeated-cross-section data frame with complete state-period cells.}
\item{outcome, cluster, time, timing}{Column names. Never treated timing is zero.}
\item{design}{Object returned by did_design.}
\item{assignment}{One cohort per state in the design's region order. Defaults to observed cohorts.}
\item{engine}{Cached calculation or raw-data refit oracle.}
\item{null}{Scalar null for the jackknife t statistic.}
}
\details{Jackknife SE uses mean-centered CV3J with factor (G-1)/G, and p-values use Student t with G-1 degrees of freedom. This adapter requires fixed group-time support; use fast_csdid for changing-support inference. Undefined deletions make inference unavailable.}
\value{A design list, scalar point estimate, or inference list with estimates, p-values and diagnostics.}
