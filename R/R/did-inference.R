# No-covariate cross-package adapter. Existing RA/IPW/DR APIs are unchanged.
.did_abort <- function(reason, detail="") {
  stop(structure(list(message=paste(reason,detail),call=NULL,reason=reason),class=c("did_inference_error","error","condition")))
}

did_design <- function(data, outcome, cluster, time, timing) {
  moments <- collapse_rcs(data,outcome,time,timing,cluster)
  if(any(moments$weight_sum<=0)) .did_abort("undefined_estimand","complete state-period cells required")
  list(data=data,outcome=outcome,cluster=cluster,time=time,timing=timing,moments=moments)
}
.did_assignment <- function(design,assignment=NULL) {
  a <- if(is.null(assignment)) design$moments$assignment else as.numeric(assignment)
  a[a==0] <- Inf
  if(length(a)!=length(design$moments$regions)||anyNA(a)||any(a<0)) .did_abort("invalid_timing_assignment")
  for(g in unique(a[is.finite(a)])) {
    if(!any(design$moments$times<g)||!any(design$moments$times>=g)) .did_abort("invalid_timing_assignment")
  }
  if(!any(is.finite(a))) .did_abort("constant_assignment")
  a
}

did_point <- function(design, assignment=NULL, engine=c("fast","reference")) {
  engine <- match.arg(engine);a <- .did_assignment(design,assignment)
  if(engine=="reference") {
    d <- design$data; d[[design$timing]] <- a[match(d[[design$cluster]],design$moments$regions)]
    fit <- fast_csdid(d,design$outcome,design$time,design$timing,design$cluster,inference="none",aggregation="simple")
    value <- fit$aggregation$estimates[["simple"]]
  } else {
    m <- design$moments;m$assignment <- a;map <- compile_att_map(m)
    v <- evaluate_assignment(m,att_map=map)
    value <- aggregate_att(v$att,map,v$group_size,"simple")$estimates[["simple"]]
  }
  if(!is.finite(value)) .did_abort("undefined_estimand")
  as.numeric(value)
}

did_jackknife <- function(design,assignment=NULL,null=0,engine=c("fast","reference")) {
  engine <- match.arg(engine);a <- .did_assignment(design,assignment);G <- length(a)
  if(G<3) .did_abort("insufficient_clusters")
  estimate <- did_point(design,a,engine);m <- design$moments;m$assignment <- a
  map <- compile_att_map(m);v <- evaluate_assignment(m,att_map=map)
  if(engine=="fast") {
    deleted <- .delete_att_matrix(m,map,v$att)$estimates
    # Native simple ATT weights are cohort total sample weights, re-estimated after deletion.
    weights <- matrix(v$group_size[as.character(map$group)],G,nrow(map),byrow=TRUE)
    for(h in seq_len(G)) weights[h,map$group==a[h]] <- weights[h,map$group==a[h]]-m$region_total_weight[h]
    theta <- rowSums(deleted*weights)/rowSums(weights)
  } else {
    theta <- vapply(seq_len(G),function(h) {
      d <- design$data[design$data[[design$cluster]]!=m$regions[h],,drop=FALSE]
      d[[design$timing]] <- a[match(d[[design$cluster]],m$regions)]
      if(is.finite(a[h])&&sum(a==a[h])==1) return(NA_real_)
      if(!any(is.infinite(a[-h]))) return(NA_real_)
      dd <- did_design(d,design$outcome,design$cluster,design$time,design$timing)
      did_point(dd,engine="reference")
    },numeric(1))
  }
  names(theta) <- as.character(m$regions)
  failed <- !is.finite(theta)
  se <- if(any(failed)) NA_real_ else sqrt((G-1)/G*sum((theta-mean(theta))^2))
  stat <- if(is.finite(se)&&se>0) (estimate-null)/se else NA_real_
  p <- if(is.finite(stat)) 2*stats::pt(-abs(stat),G-1) else NA_real_
  status <- if(any(failed)) "singular_jackknife_replication" else if(!is.finite(se)||se==0) "zero_standard_error" else "ok"
  list(estimate=estimate,standard_error=se,statistic=stat,p_value=p,method="cluster_jackknife",cluster_count=G,
       jackknife_standard_error=se,jackknife_t=stat,jackknife_df=G-1,jackknife_p_value=p,
       delete_estimates=theta,diagnostics=list(status=status,failed_clusters=m$regions[failed],centering="mean of delete-one estimates",engine=engine))
}
