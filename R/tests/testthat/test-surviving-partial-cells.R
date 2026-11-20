test_that("changing support retains identified cells in each deletion", {
 d <- read.csv(test_path("fixtures","partial_deletion.csv"))
 for(method in c("none","reg","ipw","dr")) {
  f <- fast_csdid(d,"y","time","g","cluster",weights="w",est_method=method,
   xformla=if(method=="none") NULL else ~x+z,inference="jackknife",jackknife_method="refit")
  deleted <- f$jackknife$aggregate[,"simple"]
  expect_length(deleted,13)
  expect_true(all(is.finite(deleted)),info=method)
  expect_true(f$jackknife$support_changed["13","simple"])
  broken <- fast_csdid(d[d$cluster!=13,],"y","time","g","cluster",weights="w",est_method=method,
   xformla=if(method=="none") NULL else ~x+z,inference="none")
  expect_false(is.finite(broken$aggregation$estimates["simple"]))
 }
})

test_that("legacy no-covariate helpers retain surviving cohorts", {
 d <- expand.grid(unit=1:4,time=1:6,cluster=1:8)
 d$g <- c(3,5,rep(0,6))[d$cluster]
 d$y <- .2*d$cluster+.1*d$time+sin(d$cluster*d$time+d$unit)+
  (d$g>0 & d$time>=d$g)*(1+.1*d$cluster)
 design <- did_design(d,"y","cluster","time","g")
 a <- did_jackknife(design,engine="fast")
 b <- did_jackknife(design,engine="reference")
 expect_true(all(is.finite(a$delete_estimates)))
 expect_equal(a$delete_estimates,b$delete_estimates,tolerance=1e-10)
 expect_equal(a$jackknife_df,7)
})
