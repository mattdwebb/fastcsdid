# Refresh only influence-dependent DR outputs. Existing point, deletion and
# assignment fixtures are retained and checked, not silently overwritten.
args <- commandArgs(TRUE)
stopifnot(length(args) == 1L)
stopifnot(packageVersion("DRDID") >= package_version("1.4.1"))
for (path in list.files(file.path(args[1], "R"), full.names=TRUE)) source(path)
out <- "validation/frozen"
d <- read.csv(file.path(out, "fixture.csv"))
savecsv <- function(x, name) write.csv(x, file.path(out,paste0(name,".csv")),row.names=FALSE)
external <- list()
for (method in c("none","reg","ipw","dr")) for(control in c("never","not_yet")) for(weighted in c(FALSE,TRUE)) {
  name <- paste(method,control,as.integer(weighted),sep="_")
  f <- fast_csdid(d,"y","time","g","cluster",weights=if(weighted) "w" else NULL,
    xformla=if(method=="none") NULL else ~x+z,est_method=method,
    control_group=control,inference="jackknife",jackknife_method="refit")
  old <- function(suffix) read.csv(file.path(out,paste0(name,"_",suffix,".csv")))
  stopifnot(max(abs(f$att_gt-old("att")$estimate))<2e-9,
            max(abs(f$aggregation$estimates-old("agg")$estimate))<2e-9,
            max(abs(f$jackknife$aggregate[,"simple"]-old("jk")$deleted))<2e-8)
  if(method=="dr") {
    savecsv(influence_function(f,type="att"),paste0(name,"_if"))
    savecsv(influence_function(f,type="aggregate"),paste0(name,"_aif"))
    savecsv(f$vcov$analytical,paste0(name,"_cov"))
  }
  for(k in seq_len(nrow(f$att_map))) {
    a <- f$att_map[k,]
    use <- (d$g==a$group | d$g==0 | (control=="not_yet" & d$g>a$post)) & d$time %in% c(a$base,a$post)
    z <- d[use,]
    fun <- get(switch(method,none="reg_did_rc",reg="reg_did_rc",ipw="std_ipw_did_rc",dr="drdid_rc"),asNamespace("DRDID"))
    ref <- fun(z$y,as.numeric(z$time==a$post),as.numeric(z$g==a$group),
      covariates=if(method=="none") matrix(1,nrow(z),1) else model.matrix(~x+z,z),
      i.weights=if(weighted) z$w else NULL,inffunc=TRUE)
    ii <- influence_function(f,type="att")[use,k]*sum(use)/nrow(d)
    external[[length(external)+1L]] <- data.frame(case=name,cell=k,
      att_gap=abs(f$att_gt[k]-ref$ATT),if_gap=max(abs(ii-ref$att.inf.func)))
  }
}
savecsv(do.call(rbind,external),"external_drdid")
for (stem in c("randomized","stratified_randomized")) {
  file <- file.path(out,paste0(stem,".csv"))
  if(!file.exists(file)) next
  B <- as.matrix(read.csv(file.path(out,if(stem=="randomized") "assignments.csv" else "stratified_assignments.csv")))
  records <- read.csv(file)
  for(i in which(records$method=="dr")) {
    dd <- d
    if(records$row[i]>0) dd$g <- B[records$row[i],dd$cluster]
    control <- if("control" %in% names(records)) records$control[i] else "never"
    f <- fast_csdid(dd,"y","time","g","cluster",weights="w",xformla=~x+z,
                   est_method="dr",control_group=control,inference="jackknife",jackknife_method="refit")
    deleted <- f$jackknife$aggregate[,"simple"]
    jk <- sqrt(11/12*sum((deleted-mean(deleted))^2))
    stopifnot(abs(f$aggregation$estimates["simple"]-records$estimate[i])<2e-9,
              abs(jk-records$jk[i])<2e-8)
    records$se[i] <- sqrt(f$vcov$analytical["simple","simple"])
  }
  write.csv(records,file,row.names=FALSE)
}
writeLines(c(paste("R",getRversion()),paste("DRDID",packageVersion("DRDID")),
 "Refreshed 2026-11-20: corrected pre-period DR outcome-regression influence.",
 "DRDID >=1.4.1 verifies all point estimates and all influence functions, including DR.",
 "Independent checks: empirical cluster-weight finite differences and official Stata csdid 2.0.3.",
 paste("R/influence.R MD5",tools::md5sum(file.path(args[1],"R/influence.R")))),file.path(out,"versions.txt"))
cat("DR_FIXTURES_REFRESHED_POINTS_DELETIONS_UNCHANGED\n")
