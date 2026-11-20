library(fastrcsdid)
stopifnot(packageVersion("DRDID")>=package_version("1.4.0"))
d <- read.csv("data.csv")
fit <- fast_csdid(d,"y","time","g","cluster",xformla=~x+z,est_method="dr",control_group="not_yet",inference="none")
N <- nrow(d); map <- fit$att_map
values <- numeric(nrow(map)); influence <- matrix(0,N,nrow(map))
for(k in seq_len(nrow(map))) {
  m <- map[k,]
  rows <- which(d$time %in% c(m$base,m$post) & (d$g==m$group | d$g==0 | d$g>m$post))
  x <- d[rows,]
  z <- DRDID::drdid_rc(x$y,as.numeric(x$time==m$post),as.numeric(x$g==m$group),model.matrix(~x+z,x),inffunc=TRUE)
  values[k] <- z$ATT
  influence[rows,k] <- N/length(rows)*z$att.inf.func
}
cell <- influence_function(fit,type="att",cluster_sum=TRUE)
reference <- rowsum(influence,d$cluster)
stopifnot(max(abs(values-fit$att_gt))<1e-9,max(abs(cell-reference))<1e-7)
mass <- vapply(map$group,function(g) sum(d$g==g),numeric(1))
weights <- mass/sum(mass); theta <- sum(weights*values)
score <- as.vector(influence%*%weights)
for(k in seq_along(values)) score <- score+N*(d$g==map$group[k])*(values[k]-theta)/sum(mass)
se <- sqrt(sum(rowsum(score,d$cluster)^2))/N
native <- sqrt(fit$vcov$analytical["simple","simple"])
stopifnot(abs(theta-fit$aggregation$estimates["simple"])<1e-9,abs(se-native)<1e-9)
write.csv(data.frame(engine=c("fast_R","DRDID"),version=c(as.character(packageVersion("fastrcsdid")),as.character(packageVersion("DRDID"))),estimate=c(unname(fit$aggregation$estimates["simple"]),theta),native_se=c(native,se)),"R_results.csv",row.names=FALSE)
write.csv(data.frame(map,fast_att=unname(fit$att_gt),upstream_att=values),"R_cells.csv",row.names=FALSE)
cat("NOTYET_DRDID_141_CELL_IF_AND_AGGREGATE_PASS\n")
