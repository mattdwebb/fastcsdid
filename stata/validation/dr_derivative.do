version 16
clear all
set more off
set type double
capture set processors 1
adopath ++ ".."
run "../fastcsdid.ado"
insheet using "frozen/fixture.csv", comma names clear
quietly fastcsdid y, time(time) gvar(g) cluster(cluster) covariates(x z) weights(w) method(dr)
mata:
Z0=st_data(.,tokens("y time g cluster w x z"))
ids=uniqrows(sort(Z0[,4],1)); step=1e-5
M=A=IF=CV=J(0,0,.)
for (weighted=0;weighted<=1;weighted++) {
  Z=Z0
  if(!weighted) Z[,5]=J(rows(Z),1,1)
  for (ci=1;ci<=2;ci++) {
    control=ci==1 ? "never" : "not_yet"
    fcs_point(Z,"dr",control,0,.,M,A,IF,CV)
    reference=CV
    D=J(rows(ids),rows(M)+rows(A),.)
    for(h=1;h<=rows(ids);h++) {
      saved=Z[,5]
      Z[,5]=saved:*(1:+step:*(Z[,4]:==ids[h]))
      fcs_point(Z,"dr",control,0,.,M,A,IF,CV)
      upper=M[,5]\A[,3]
      Z[,5]=saved:*(1:-step:*(Z[,4]:==ids[h]))
      fcs_point(Z,"dr",control,0,.,M,A,IF,CV)
      lower=M[,5]\A[,3]
      Z[,5]=saved
      D[h,]=(upper-lower)'/(2*step)
    }
    assert(max(abs(vec(reference-quadcross(D,D))))<2e-8)
  }
}
end
display "CSDID_DR_DERIVATIVE_PASS"
exit, clear
