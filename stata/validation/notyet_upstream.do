version 16
clear all
set more off
adopath ++ "."
insheet using "validation/frozen/notyet_upstream_reference.csv", comma names clear double
scalar reference_att=estimate[1]
scalar reference_se=native_se[1]
insheet using "validation/frozen/notyet_upstream_data.csv", comma names clear double
fastcsdid y, time(time) gvar(g) cluster(cluster) covariates(x z) method(dripw) control(not_yet)
assert abs(_b[ATT]-reference_att)<1e-8
assert abs(e(native_se)-reference_se)<1e-8
display "NOTYET_UPSTREAM_REGRESSION_PASS"
