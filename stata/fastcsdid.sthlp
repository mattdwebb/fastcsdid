{smcl}
{* *! version 1.0.0 20oct2026}{...}
{title:fastcsdid — Group-time DiD for stationary repeated cross sections}

{title:Syntax}
{p 8 12 2}
{cmd:fastcsdid} {it:outcome} {ifin},
{cmd:time(}{it:timevar}{cmd:)} {cmd:gvar(}{it:cohortvar}{cmd:)}
{cmd:cluster(}{it:clustervar}{cmd:)}
[{cmd:covariates(}{it:varlist}{cmd:)} {cmd:weights(}{it:weightvar}{cmd:)}
{cmd:method(}{it:method}{cmd:)} {cmd:control(}{it:group}{cmd:)}
{cmd:eventwindow(}{it:min max}{cmd:)} {cmd:jackknife}
{cmd:jksupport(}{it:policy}{cmd:)}]

{title:Description}
{pstd}
Estimates post-treatment ATT(g,t) and simple, cohort, calendar and event-time
aggregations. The command always uses repeated-cross-section equations.
Each cluster must have one absorbing treatment cohort; zero means never treated.
Covariates must be numeric and already encoded.

{title:Options}
{pstd}
{cmd:method()} accepts {cmd:none}, {cmd:reg} (alias {cmd:ra}),
{cmd:ipw} (alias {cmd:stdipw}), and {cmd:dr} (alias {cmd:dripw}, default).
{cmd:none} rejects covariates. {cmd:control()} is {cmd:never} (default) or
{cmd:not_yet}. {cmd:weights()} supplies nonnegative sampling weights.
Required missing data are rejected.

{pstd}
{cmd:eventwindow()} restricts the nonnegative post-treatment event-time window.
{cmd:jackknife} computes complete delete-cluster refits, mean-centered CV3J
standard errors and Student t(G-1) p-values.
{cmd:jksupport(changing)} (default) recomputes cohort masses and allows a cohort
to disappear if the remaining aggregate is identified.
{cmd:jksupport(fixed)} additionally requires unchanged group-time support.
An undefined deletion causes an error; failed deletions are not dropped.

{title:Stored results}
{pstd}
{cmd:e(b)}: simple ATT. {cmd:e(V)}: analytical cluster influence variance.
{cmd:e(native_se)}, {cmd:e(native_p)}: analytical standard error and normal-reference p-value.
{cmd:e(N)}, {cmd:e(N_clust)}: observations and clusters.
{cmd:e(attgt)} columns: group, base, post, event, estimate.
{cmd:e(aggregation)} columns: family (1=simple, 2=group, 3=calendar, 4=event),
level, estimate.
{cmd:e(deleted)}, {cmd:e(jk_se)}, {cmd:e(jk_p)}, {cmd:e(jk_df)}:
optional jackknife results. {cmd:e(V)} remains analytical.

{title:Example}
{phang2}{cmd:fastcsdid y, time(year) gvar(first_treat) cluster(state) method(reg) covariates(x z) jackknife}{p_end}

{title:Author and source}
{pstd}Matthew D. Webb. MIT license.
{browse "https://github.com/mattdwebb/fastcsdid":GitHub repository}.
