test_that("not-yet DRIPW matches DRDID 1.4.1 and official Stata csdid 2.0.3", {
  d <- read.csv(test_path("fixtures","notyet_upstream_data.csv"))
  reference <- read.csv(test_path("fixtures","notyet_upstream_reference.csv"))
  z <- fast_csdid(d,"y","time","g","cluster",xformla=~x+z,est_method="dr",control_group="not_yet",inference="none")
  expect_lt(max(abs(z$aggregation$estimates["simple"]-reference$estimate)),1e-8)
  expect_lt(max(abs(sqrt(z$vcov$analytical["simple","simple"])-reference$native_se)),1e-8)
})
