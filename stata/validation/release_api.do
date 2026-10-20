version 16
clear all
set more off
adopath ++ ".."
import delimited using "frozen/fixture.csv", clear
foreach opt in rijack assignments(B) strata(cluster) {
    capture fastcsdid y, time(time) gvar(g) cluster(cluster) `opt'
    assert _rc == 198
}
capture which fastcsdid_assignments
* Check removal from this distribution via the repository source audit;
* another installed package may still have the old helper on adopath.
fastcsdid y, time(time) gvar(g) cluster(cluster) method(none) jackknife
assert e(jk_se) > 0 & e(jk_se) < .
preserve
drop if inlist(cluster,2,3,4)
fastcsdid y, time(time) gvar(g) cluster(cluster) method(none) jackknife
assert rowsof(e(deleted)) == 9
capture fastcsdid y, time(time) gvar(g) cluster(cluster) method(none) jackknife jksupport(fixed)
assert _rc == 498
restore
display "FASTCSDID_RELEASE_API_PASS"
exit, clear
