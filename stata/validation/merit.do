version 16
clear all
set more off
set type double
adopath ++ ".."

import delimited using "frozen/merit.csv", clear
fastcsdid outcome, time(year) gvar(cohort) cluster(state) weights(weight) ///
    method(none) jackknife

assert abs(_b[ATT] - 0.0477819967004955) < 1e-12
assert abs(e(jk_se) - 0.0146168254236831) < 1e-12
assert rowsof(e(attgt)) == 33
assert rowsof(e(deleted)) == 51
assert e(jackknife_support) == "changing"
display "FASTCSDID_STATA_MERIT_COMPLETE"
exit, clear
