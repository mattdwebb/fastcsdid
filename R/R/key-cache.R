# Environment-backed caches keyed by long strings.
#
# R environment names are symbols, which are limited to 10000 bytes. Cache
# keys built from comma-joined cluster/region ids exceed that limit once there
# are a few thousand clusters (e.g. one cluster per row), so keys are hashed to
# a short bucket name. Each bucket stores the full keys and values, and lookups
# compare the full key, so hash collisions can never alias two entries.

.key_cache_bucket <- function(key) {
  bytes <- as.numeric(utf8ToInt(enc2utf8(key)))
  n <- length(bytes)
  if (n == 0L) return("h0:0")
  idx <- seq_len(n)
  h1 <- sum((bytes * ((idx %% 65521) + 1)) %% 2147483647) %% 2147483647
  h2 <- sum((bytes * ((rev(idx) %% 65519) + 7)) %% 2147483629) %% 2147483629
  sprintf("h%.0f:%.0f:%d", h1, h2, n)
}

.key_cache_find <- function(env, key, bucket = .key_cache_bucket(key)) {
  if (!exists(bucket, envir = env, inherits = FALSE)) return(0L)
  entry <- get(bucket, envir = env, inherits = FALSE)
  pos <- match(key, entry$keys)
  if (is.na(pos)) 0L else pos
}

.key_cache_has <- function(env, key) {
  .key_cache_find(env, key) > 0L
}

.key_cache_set <- function(env, key, value) {
  bucket <- .key_cache_bucket(key)
  entry <- if (exists(bucket, envir = env, inherits = FALSE)) {
    get(bucket, envir = env, inherits = FALSE)
  } else {
    list(keys = character(), values = list())
  }
  pos <- match(key, entry$keys)
  if (is.na(pos)) {
    pos <- length(entry$keys) + 1L
    entry$keys[pos] <- key
  }
  entry$values[[pos]] <- value
  assign(bucket, entry, envir = env)
  invisible(value)
}

.key_cache_get <- function(env, key) {
  bucket <- .key_cache_bucket(key)
  pos <- .key_cache_find(env, key, bucket)
  if (pos == 0L) stop("internal error: cache key not found", call. = FALSE)
  get(bucket, envir = env, inherits = FALSE)$values[[pos]]
}
