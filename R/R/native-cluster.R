#' Native no-covariate CSDID clustered analytical inference
#' @param design A did_design object.
#' @param null Scalar null value.
#' @export
did_native_cluster <- function(design, null=0) {
  fit <- fast_csdid(design$data,design$outcome,design$time,design$timing,design$cluster,
                    inference='none',aggregation='simple')
  theta <- unname(fit$aggregation$estimates['simple'])
  se <- sqrt(fit$vcov$analytical['simple','simple'])
  stat <- (theta-null)/se
  list(estimate=theta,native_cluster_method='CSDID_cluster_influence_analytical',
       native_cluster_se=se,native_cluster_stat=stat,native_cluster_df=NA_real_,
       native_cluster_p=2*stats::pnorm(-abs(stat)),native_cluster_distribution='normal',
       status=if(is.finite(se)&&se>0) 'ok' else 'zero_standard_error')
}
