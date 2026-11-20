version 16
clear all
set more off
set type double
capture set processors 1
local official : env CSDID_OFFICIAL_ADO
if "`official'"=="" {
    display as error "Set CSDID_OFFICIAL_ADO to the frozen official csdid 2.0.3 directory."
    exit 198
}
adopath ++ "`official'"
adopath ++ "../../stata"
csdid version
insheet using "data.csv", comma names clear double
tempname output
postfile `output' str20 engine double estimate native_se using "stata_results.dta", replace
csdid y x z, time(time) gvar(g) cluster(cluster) rcs notyet base_period(universal) method(dr) analytical pointwise
estat simple, post
scalar reference_att=_b[ATT]
scalar reference_se=_se[ATT]
post `output' ("official_Stata") (reference_att) (reference_se)
fastcsdid y, time(time) gvar(g) cluster(cluster) covariates(x z) method(dripw) control(not_yet)
assert abs(_b[ATT]-reference_att)<1e-8
assert abs(e(native_se)-reference_se)<1e-8
post `output' ("fast_Stata") (_b[ATT]) (e(native_se))
postclose `output'
use "stata_results.dta", clear
format estimate native_se %24.17g
outsheet using "Stata_results.csv", comma replace
display "NOTYET_OFFICIAL_STATA_203_PASS"
