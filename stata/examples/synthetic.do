version 16
clear
set type double
insheet using "validation/frozen/fixture.csv", comma names clear
adopath ++ "."
fastcsdid y, time(time) gvar(g) cluster(cluster) covariates(x z) weights(w) jackknife
matrix list e(attgt)
matrix list e(aggregation)
