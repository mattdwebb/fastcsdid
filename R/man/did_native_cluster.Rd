\name{did_native_cluster}
\alias{did_native_cluster}
\title{Native analytical state-clustered CSDID inference}
\usage{did_native_cluster(design, null = 0)}
\arguments{\item{design}{A no-covariate did_design object.}\item{null}{Scalar null value.}}
\description{Computes the simple-aggregate influence-function clustered covariance, including estimated aggregation weights. Uses normal p-values and no additional finite-sample multiplier.}
\value{A list with estimate, native cluster method, SE, statistic, df (NA), p-value, distribution and status.}
