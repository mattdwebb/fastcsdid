test_that("native analytical inference exposes the existing clustered covariance", {
 d <- make_test_data(FALSE)
 design <- did_design(d,"y","region","time","cohort")
 z <- did_native_cluster(design)
 f <- fast_csdid(d,"y","time","cohort","region",inference="none",aggregation="simple")
 expect_equal(z$native_cluster_se,sqrt(f$vcov$analytical["simple","simple"]))
 expect_equal(z$native_cluster_p,2*pnorm(-abs(z$native_cluster_stat)))
 expect_equal(z$estimate,did_point(design))
})
