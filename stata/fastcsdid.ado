*! fastcsdid 1.0.0 20oct2026 MIT; Matt Webb
program define fastcsdid, eclass
    version 16.0
    syntax varname(numeric) [if] [in], Time(varname numeric) Gvar(varname numeric) ///
        Cluster(varname) [Covariates(varlist numeric) Weights(varname numeric) ///
        Method(string) Control(string) EVENTWindow(numlist min=2 max=2) JACKknife ///
        JKSUPport(string)]
    if "`method'" == "" local method dripw
    if "`method'" == "dripw" local method dr
    if "`method'" == "ra" local method reg
    if "`method'" == "stdipw" local method ipw
    if !inlist("`method'","none","reg","ipw","dr") error 198
    if "`control'" == "" local control never
    if !inlist("`control'","never","not_yet") error 198
    if "`jksupport'" == "" local jksupport changing
    if !inlist("`jksupport'","changing","fixed") {
        di as error "jksupport() must be changing or fixed"
        exit 198
    }
    if "`method'" == "none" & "`covariates'" != "" {
        di as error "method(none) cannot ignore covariates"
        exit 198
    }
    local lo 0
    local hi .
    if "`eventwindow'" != "" {
        tokenize `eventwindow'
        local lo `1'
        local hi `2'
        if `lo' < 0 | `lo' > `hi' error 198
    }
    marksample use, novarlist
    tempvar cid w bad
    quietly egen long `cid' = group(`cluster') if `use'
    quietly gen double `w' = 1 if `use'
    if "`weights'" != "" quietly replace `w' = `weights' if `use'
    quietly egen byte `bad' = rowmiss(`varlist' `time' `gvar' `cid' `w' `covariates') if `use'
    quietly count if `use' & (`bad' > 0 | `w' < 0 | `gvar' < 0)
    if r(N) {
        di as error "Required data must be observed; weights and cohort nonnegative"
        exit 459
    }
    tempname b V M A J NV NP JP JS NC
    mata: fcs_run("`varlist' `time' `gvar' `cid' `w' `covariates'", "`use'", ///
        "`method'", "`control'", `lo', `hi', "`jackknife'", "`jksupport'", ///
        "`b'", "`V'", "`M'", "`A'", "`J'", "`NV'", "`NP'", "`JS'", "`JP'", "`NC'")
    matrix colnames `b' = ATT
    matrix colnames `V' = ATT
    matrix rownames `V' = ATT
    ereturn post `b' `V', esample(`use')
    ereturn matrix attgt = `M'
    ereturn matrix aggregation = `A'
    ereturn scalar native_se = scalar(`NV')
    ereturn scalar native_p = scalar(`NP')
    ereturn scalar N_clust = scalar(`NC')
    quietly count if e(sample)
    ereturn scalar N = r(N)
    if "`jackknife'" != "" {
        ereturn matrix deleted = `J'
        ereturn scalar jk_se = scalar(`JS')
        ereturn scalar jk_p = scalar(`JP')
        ereturn scalar jk_df = scalar(`NC')-1
    }
    ereturn local method "`method'"
    ereturn local jackknife_support "`jksupport'"
    ereturn local data_mode "rcs"
    ereturn local target "cohort-mass weighted group-time ATT"
    ereturn local status "pass"
    ereturn local cmd "fastcsdid"
    di as text "RCS aggregate ATT = " as result %12.8f _b[ATT]
    di as text "Native SE (normal reference) = " as result %12.8f e(native_se)
end

mata:
real scalar fcs_ratio(real colvector w, real colvector y)
{
    if (sum(w)<=0) return(.)
    return(sum(w:*y)/sum(w))
}
real colvector fcs_rif(real colvector w, real colvector y)
{
    return(w:*(y:-fcs_ratio(w,y))/mean(w))
}
real scalar fcs_wls(real matrix X, real colvector y, real colvector w,
    real colvector m, real matrix IF)
{
    real matrix H
    real colvector b
    H = quadcross(X,w,X)
    if (rank(H)<cols(X)) return(0)
    b = qrsolve(H,quadcross(X,w,y))
    m = X*b
    IF = qrsolve(H/rows(X),(X:*(w:*(y-m)))')'
    return(1)
}
real scalar fcs_logit(real matrix X, real colvector D, real colvector w,
    real colvector p, real matrix IF)
{
    real colvector b, step, eta, bn
    real matrix H
    real scalar iter, scale, loss, newer
    b = J(cols(X),1,0)
    for(iter=1;iter<=100;iter++) {
        p = invlogit(X*b)
        H = quadcross(X,w:*p:*(1:-p),X)
        if(rank(H)<cols(X)) return(0)
        step = qrsolve(H,quadcross(X,w,D-p))
        eta = X*b
        loss = sum(w:*(ln(1:+exp(eta)):-D:*eta))
        scale = 1
        while(scale>2^-25) {
            bn = b+scale*step
            eta = X*bn
            newer = sum(w:*(ln(1:+exp(eta)):-D:*eta))
            // The summed objective grows with N. Allow only floating-point
            // roundoff in its comparison; keep the Newton-step convergence
            // threshold unchanged.
            if(newer<=loss+64*2^-52*max((1,abs(loss),abs(newer)))) break
            scale = scale/2
        }
        b = bn
        if(max(abs(step))<1e-10) break
        if(max(abs(b))>40) return(0)
    }
    if(iter>100) return(0)
    p = rowmin((invlogit(X*b),J(rows(X),1,1-1e-6)))
    H = quadcross(X,w:*p:*(1:-p),X)
    if(rank(H)<cols(X)) return(0)
    IF = qrsolve(H/rows(X),(X:*(w:*(D-p)))')'
    return(1)
}
real colvector fcs_cell(real colvector y, real colvector P, real colvector D,
    real matrix X, real colvector w, string scalar method)
{
    real colvector w10,w11,w00,w01,wd,p,trim,z,m00,m01,m10,m11,inf,mp,mb
    real matrix o00,o01,o10,o11,pi
    real scalar att,r10,r11,r00,r01
    if(rows(y)==0 | sum(w)<=0) return(J(rows(y)+1,1,.))
    w=w/mean(w)
    w10=w:*D:*(1:-P); w11=w:*D:*P
    w00=w:*(1:-D):*(1:-P); w01=w:*(1:-D):*P
    if(min((sum(w10),sum(w11),sum(w00),sum(w01)))<=0) return(J(rows(y)+1,1,.))
    if(method=="none") {
        att=fcs_ratio(w11,y)-fcs_ratio(w10,y)-fcs_ratio(w01,y)+fcs_ratio(w00,y)
        inf=fcs_rif(w11,y)-fcs_rif(w10,y)-fcs_rif(w01,y)+fcs_rif(w00,y)
        return(att\inf)
    }
    m00=m01=m10=m11=J(rows(y),1,.)
    o00=o01=o10=o11=pi=J(rows(y),cols(X),.)
    if(method=="reg" | method=="dr") {
        if(!fcs_wls(X,y,w00,m00,o00) | !fcs_wls(X,y,w01,m01,o01)) return(J(rows(y)+1,1,.))
    }
    if(method=="reg") {
        wd=w:*D
        att=fcs_ratio(w11,y)-fcs_ratio(w10,y)-fcs_ratio(wd,m01-m00)
        inf=fcs_rif(w11,y)-fcs_rif(w10,y)-fcs_rif(wd,m01-m00)-(o01-o00)*(X'*wd/sum(wd))
        return(att\inf)
    }
    p=J(rows(y),1,.)
    if(!fcs_logit(X,D,w,p,pi)) return(J(rows(y)+1,1,.))
    trim=(D:==1):|(p:<.995)
    w10=trim:*w10; w11=trim:*w11
    w00=trim:*p:/(1:-p):*w00; w01=trim:*p:/(1:-p):*w01
    if(min((sum(w00),sum(w01)))<=0) return(J(rows(y)+1,1,.))
    z=y
    if(method=="dr") z=y-P:*m01-(1:-P):*m00
    r10=fcs_ratio(w10,z); r11=fcs_ratio(w11,z)
    r00=fcs_ratio(w00,z); r01=fcs_ratio(w01,z)
    att=r11-r10-r01+r00
    inf=fcs_rif(w11,z)-fcs_rif(w10,z)-fcs_rif(w01,z)+fcs_rif(w00,z) -
        pi*(X'*(w01:*(z:-r01))/sum(w01)-X'*(w00:*(z:-r00))/sum(w00))
    if(method=="dr") {
        if(!fcs_wls(X,y,w:*D:*(1:-P),m10,o10) | !fcs_wls(X,y,w:*D:*P,m11,o11)) return(J(rows(y)+1,1,.))
        wd=trim:*w:*D
        att=att+fcs_ratio(wd,m11-m01)-fcs_ratio(w11,m11-m01)-fcs_ratio(wd,m10-m00)+fcs_ratio(w10,m10-m00)
        mp=X'*(wd/sum(wd)-w11/sum(w11))
        mb=X'*(wd/sum(wd)-w10/sum(w10))
        inf=inf+o01*(X'*(w01/sum(w01)-w11/sum(w11)))+o00*(X'*(w00/sum(w00)-w10/sum(w10))) +
            fcs_rif(wd,m11-m01)-fcs_rif(w11,m11-m01)-fcs_rif(wd,m10-m00)+fcs_rif(w10,m10-m00) +
            (o11-o01)*mp-(o10-o00)*mb
    }
    return(att\inf)
}
void fcs_point(real matrix Z, string scalar method, string scalar control,
    real scalar lo, real scalar hi, real matrix M, real matrix A, real matrix IF, real matrix CV)
{
    real colvector y,t,g,c,w,groups,times,pre,post,idx,v,inf,mass,sel,raw,W,active,ids,levels
    real matrix X,I,JIF,CI
    real scalar N,k,j,h,b,period,group,event,a,typ,lev,theta,den
    N=rows(Z); y=Z[,1]; t=Z[,2]; g=Z[,3]; c=Z[,4]; w=Z[,5]
    X=J(N,1,1)
    if(cols(Z)>5) X=(X,Z[,6..cols(Z)])
    groups=uniqrows(sort(select(g,g:>0),1)); times=uniqrows(sort(t,1))
    M=J(0,5,.); I=J(N,0,.)
    for(j=1;j<=rows(groups);j++) {
        group=groups[j]
        pre=select(times,times:<group); post=select(times,times:>=group)
        if(!rows(pre) | !rows(post)) continue
        b=max(pre)
        for(k=1;k<=rows(post);k++) {
            period=post[k]; event=period-group
            if(event<lo | event>hi) continue
            sel=(g:==group):|(g:==0)
            if(control=="not_yet") sel=sel:|(g:>period)
            idx=selectindex(sel:&((t:==b):|(t:==period)))
            v=fcs_cell(y[idx],t[idx]:==period,g[idx]:==group,X[idx,],w[idx],method)
            inf=J(N,1,.)
            if(v[1]<.) {
                inf=J(N,1,0)
                inf[idx]=N/rows(idx)*v[2..rows(v)]
            }
            M=M\(group,b,period,event,v[1]); I=(I,inf)
        }
    }
    if(!rows(M)) _error(2000)
    mass=J(rows(M),1,0)
    for(k=1;k<=rows(M);k++) mass[k]=sum(select(w,g:==M[k,1]))
    A=J(0,3,.); IF=J(N,0,.)
    for(typ=1;typ<=4;typ++) {
        levels=.
        if(typ==2) levels=uniqrows(sort(M[,1],1))
        if(typ==3) levels=uniqrows(sort(M[,3],1))
        if(typ==4) levels=uniqrows(sort(M[,4],1))
        for(lev=1;lev<=rows(levels);lev++) {
            sel=J(rows(M),1,1)
            if(typ==2) sel=M[,1]:==levels[lev]
            if(typ==3) sel=M[,3]:==levels[lev]
            if(typ==4) sel=M[,4]:==levels[lev]
            raw=sel:*mass
            if(typ==2) raw=sel
            active=selectindex(raw:>0)
            den=sum(raw); theta=.
            inf=J(N,1,.)
            if(rows(active)>0 & den>0) {
                if(max(M[active,5])<.) {
                    W=raw/den
                    theta=W[active]'*M[active,5]
                    inf=I[,active]*W[active]
                    if(typ!=2) {
                        for(j=1;j<=rows(groups);j++) {
                            idx=selectindex((M[,1]:==groups[j]):&sel)
                            if(rows(idx)) inf=inf+N*w:*(g:==groups[j])*sum(M[idx,5]:-theta)/den
                        }
                    }
                }
            }
            A=A\(typ,levels[lev],theta); IF=(IF,inf)
        }
    }
    JIF=(I,IF); ids=uniqrows(sort(c,1)); CI=J(rows(ids),cols(JIF),.)
    for(h=1;h<=rows(ids);h++) CI[h,]=colsum(select(JIF,c:==ids[h]))
    CV=quadcross(CI,CI)/N^2
}
void fcs_run(string scalar vars, string scalar use, string scalar method,
    string scalar control, real scalar lo, real scalar hi, string scalar jk,
    string scalar support,
    string scalar bn,string scalar vn,string scalar mn,string scalar an,string scalar jn,
    string scalar sn,string scalar pn,string scalar jsn,string scalar jpn,string scalar cn)
{
    real matrix Z,M,A,IF,CV,DM,DA,DI,DV
    real colvector ids,deleted,gg
    real scalar h,G,v,se,theta,js
    Z=st_data(.,tokens(vars),use)
    ids=uniqrows(sort(Z[,4],1)); G=rows(ids)
    for(h=1;h<=G;h++) {
        gg=select(Z[,3],Z[,4]:==ids[h])
        if(min(gg)!=max(gg)) _error(459)
    }
    fcs_point(Z,method,control,lo,hi,M,A,IF,CV)
    theta=A[1,3]; v=CV[rows(M)+1,rows(M)+1]; se=sqrt(v)
    if(theta>=. | se>=.) _error(498)
    st_matrix(bn,theta); st_matrix(vn,v); st_matrix(mn,M); st_matrix(an,A)
    st_numscalar(sn,se); st_numscalar(pn,se>0 ? 2*normal(-abs(theta/se)) : .)
    st_numscalar(cn,G)
    if(jk!="") {
        if(G<2) _error(498)
        deleted=J(G,1,.)
        for(h=1;h<=G;h++) {
            fcs_point(select(Z,Z[,4]:!=ids[h]),method,control,lo,hi,DM,DA,DI,DV)
            if(support=="fixed") {
                if(rows(DM)!=rows(M)) _error(498)
                if(any(DM[,1..4]:!=M[,1..4])) _error(498)
            }
            deleted[h]=DA[1,3]
        }
        if(max(deleted)>=.) _error(498)
        js=sqrt((G-1)/G*sum((deleted:-mean(deleted)):^2))
        st_matrix(jn,deleted); st_numscalar(jsn,js)
        st_numscalar(jpn,js>0 ? 2*ttail(G-1,abs(theta/js)) : .)
    }
}
real colvector fcs_delete(real matrix Z,string scalar method,string scalar control,
    real scalar lo,real scalar hi,real matrix M,string scalar support)
{
    real colvector ids,deleted
    real matrix DM,DA,DI,DV
    real scalar h
    ids=uniqrows(sort(Z[,4],1)); deleted=J(rows(ids),1,.)
    for(h=1;h<=rows(ids);h++) {
        fcs_point(select(Z,Z[,4]:!=ids[h]),method,control,lo,hi,DM,DA,DI,DV)
        if(support=="fixed") {
            if(rows(DM)!=rows(M)) _error(498)
            if(any(DM[,1..4]:!=M[,1..4])) _error(498)
        }
        deleted[h]=DA[1,3]
    }
    if(any(deleted:>=.)) _error(498)
    return(deleted)
}
end
