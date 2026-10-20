args <- commandArgs(TRUE)
source_file <- args[1L]
method <- args[2L]
repetitions <- as.integer(args[3L])
output <- args[4L]
# Use the checked-in source; installation and package startup are outside timing.
for (file in list.files("R/R", full.names = TRUE, pattern = "\\.R$")) source(file)
d <- read.csv(source_file)
results <- list()
for (rep in 0:repetitions) {
  gc()
  start <- proc.time()[["elapsed"]]
  fit <- fast_csdid(d, "y", "time", "g", "cluster",
                   xformla = if (method == "none") NULL else ~x + z,
                   est_method = method, control_group = "never", inference = "none")
  elapsed <- proc.time()[["elapsed"]] - start
  results[[length(results) + 1L]] <- data.frame(
    engine = "fast_R", n = nrow(d), method = method, rep = rep, seconds = elapsed,
    estimate = unname(fit$aggregation$estimates[["simple"]]),
    native_se = sqrt(fit$vcov$analytical["simple", "simple"]))
  write.csv(do.call(rbind, results), output, row.names = FALSE)
}
