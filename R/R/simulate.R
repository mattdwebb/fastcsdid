# Internal deterministic fixture used by examples and tests.
.make_fixture <- function(unbalanced = FALSE) {
  regions <- 1:8
  times <- 1:6
  cohort <- c(3, 3, 5, 5, 0, 0, 0, 0)
  rows <- list()
  idx <- 0L
  for (h in regions) {
    for (t in times) {
      n <- if (unbalanced) 3L + ((h + 2L * t) %% 5L) else 5L
      for (j in seq_len(n)) {
        idx <- idx + 1L
        treated <- cohort[h] > 0 && t >= cohort[h]
        y <- 0.4 * h + 0.2 * t + 1.5 * treated + 0.03 * j
        rows[[idx]] <- data.frame(id = idx, region = h, time = t,
                                  cohort = cohort[h], y = y)
      }
    }
  }
  do.call(rbind, rows)
}
