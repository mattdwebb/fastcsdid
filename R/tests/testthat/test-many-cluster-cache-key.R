test_that("key cache handles long keys and exact collisions", {
  env <- new.env(parent = emptyenv())
  long <- paste(seq_len(5000), collapse = ",")
  expect_gt(nchar(long, type = "bytes"), 10000)
  .key_cache_set(env, long, 1L)
  expect_true(.key_cache_has(env, long))
  expect_identical(.key_cache_get(env, long), 1L)
  expect_false(.key_cache_has(env, paste0(long, ",5001")))
  # Same multiset of bytes in a different order must not alias.
  a <- "12|1,2"; b <- "21|1,2"
  .key_cache_set(env, a, "a"); .key_cache_set(env, b, "b")
  expect_identical(.key_cache_get(env, a), "a")
  expect_identical(.key_cache_get(env, b), "b")
})

test_that("one cluster per row works with not_yet + dr (long cache keys)", {
  set.seed(20261016)
  n_per <- 1200L
  times <- 1:4
  d <- do.call(rbind, lapply(times, function(t) {
    cohort <- sample(c(0, 3, 4), n_per, replace = TRUE)
    x1 <- rnorm(n_per)
    data.frame(time = t, cohort = cohort, x1 = x1,
               y = 0.5 * x1 + 0.2 * t + 1 * (cohort > 0 & t >= cohort) +
                 rnorm(n_per, sd = 0.5))
  }))
  d$region <- seq_len(nrow(d))
  fit <- fast_csdid(
    d, "y", "time", "cohort", "region", xformla = ~x1,
    est_method = "dr", control_group = "not_yet", inference = "none"
  )
  expect_s3_class(fit, "fastcsdid")
  att <- fit[["att_gt"]]
  expect_true(length(att) > 0)
  expect_true(all(is.finite(unlist(att))))
  reg <- fast_csdid(
    d, "y", "time", "cohort", "region", xformla = ~x1,
    est_method = "reg", control_group = "not_yet", inference = "none"
  )
  expect_s3_class(reg, "fastcsdid")
})
