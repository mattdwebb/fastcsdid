test_that("cache buckets distinguish real collisions and long UTF-8 keys", {
  cache <- new.env(parent = emptyenv())
  expect_identical(.key_cache_bucket("abc"), .key_cache_bucket("b`d"))
  .key_cache_set(cache, "abc", 1L)
  expect_false(.key_cache_has(cache, "b`d"))
  .key_cache_set(cache, "b`d", 2L)
  .key_cache_set(cache, "abc", 3L)
  expect_identical(.key_cache_get(cache, "abc"), 3L)
  expect_identical(.key_cache_get(cache, "b`d"), 2L)
  expect_length(get(.key_cache_bucket("abc"), cache)$keys, 2L)
  keys <- c("", "\u00e9\u4e2d", strrep("\u00e9", 6000L),
            paste0(strrep("\u00e9", 6000L), "-different"))
  expect_gt(nchar(keys[3L], type = "bytes"), 10000L)
  for (i in seq_along(keys)) .key_cache_set(cache, keys[i], i)
  for (i in seq_along(keys)) expect_identical(.key_cache_get(cache, keys[i]), i)
  expect_error(.key_cache_get(cache, "missing"), "cache key not found")
  expect_true(all(nchar(ls(cache), type = "bytes") < 100L))
})

audit_key_data <- function() {
  set.seed(914)
  d <- expand.grid(region = 1:12, time = 1:4, person = 1:8)
  d$cohort <- rep(c(0, 3, 4), each = 4)[d$region]
  d$x <- rnorm(nrow(d))
  d$y <- d$region/5 + d$time/3 + d$x/2 +
    2*(d$cohort > 0 & d$time >= d$cohort) + rnorm(nrow(d))/4
  d
}

audit_key_fit <- function(d, method) {
  fit <- fast_csdid(d, "y", "time", "cohort", "region", xformla = ~x,
                    est_method = method, control_group = "not_yet",
                    inference = "jackknife")
  # Strip only the region-label attributes that intentionally change.
  jackknife <- fit$jackknife$aggregate_covariance
  jackknife$valid_deletions <- unname(jackknife$valid_deletions)
  jackknife$common_support_rows <- unname(jackknife$common_support_rows)
  # Compare statistical outputs, not labels or execution timings.
  list(att = fit$att_gt, aggregate = fit$aggregation$estimates,
       covariance = fit$vcov$analytical,
       jackknife = jackknife)
}

test_that("long region labels preserve estimates and inference", {
  d <- audit_key_data()
  labelled <- d
  labelled$region <- paste0(sprintf("r%04d-", d$region), strrep("\u00e9", 6000L))
  expect_gt(nchar(labelled$region[1L], type = "bytes"), 10000L)
  for (method in c("reg", "ipw", "dr")) {
    expect_identical(audit_key_fit(d, method), audit_key_fit(labelled, method))
  }
})

test_that("forced cache collisions preserve estimates and inference", {
  d <- audit_key_data()
  methods <- c("reg", "ipw", "dr")
  expected <- lapply(methods, function(method) audit_key_fit(d, method))
  local_mocked_bindings(.key_cache_bucket = function(key) "all-keys-collide",
                        .package = "fastrcsdid")
  actual <- lapply(methods, function(method) audit_key_fit(d, method))
  expect_identical(actual, expected)
})

test_that("hashed caches preserve the pre-fix full-key results on short keys", {
  d <- audit_key_data()
  methods <- c("reg", "ipw", "dr")
  expected <- lapply(methods, function(method) audit_key_fit(d, method))
  local_mocked_bindings(
    .key_cache_has = function(env, key) exists(key, env, inherits = FALSE),
    .key_cache_get = function(env, key) get(key, env, inherits = FALSE),
    .key_cache_set = function(env, key, value) {
      assign(key, value, env)
      invisible(value)
    }, .package = "fastrcsdid")
  actual <- lapply(methods, function(method) audit_key_fit(d, method))
  expect_identical(actual, expected)
})
