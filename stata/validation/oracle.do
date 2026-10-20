version 16
clear all
set more off
set type double
capture set processors 1
adopath ++ ".."
insheet using "frozen/fixture.csv", comma names clear
tempname results components M A
postfile `results' str8 method str8 control byte weighted double estimate double se double jk using actual, replace
postfile `components' str8 method str8 control byte weighted str3 kind int index double estimate using components, replace
foreach method in none reg ipw dr {
    local cov
    if "`method'"!="none" local cov covariates(x z)
    foreach control in never not_yet {
        forvalues weighted=0/1 {
            local weight
            if `weighted' local weight weights(w)
            fastcsdid y, time(time) gvar(g) cluster(cluster) method(`method') control(`control') `cov' `weight' jackknife
            post `results' ("`method'") ("`control'") (`weighted') (_b[ATT]) (e(native_se)) (e(jk_se))
            matrix `M' = e(attgt)
            matrix `A' = e(aggregation)
            forvalues k=1/`=rowsof(`M')' {
                post `components' ("`method'") ("`control'") (`weighted') ("att") (`k') (`M'[`k',5])
            }
            forvalues k=1/`=rowsof(`A')' {
                post `components' ("`method'") ("`control'") (`weighted') ("agg") (`k') (`A'[`k',3])
            }
        }
    }
}
postclose `results'
postclose `components'
use actual, clear
format estimate se jk %24.17g
outsheet using actual.csv, comma replace
use components, clear
format estimate %24.17g
outsheet using components.csv, comma replace
display "CSDID_STATA_ORACLE_COMPLETE"
exit, clear
