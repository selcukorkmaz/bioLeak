# Time-series block permutation utilities --------------------------------------

#' Circular block permutation indices
#'
#' Generates a permutation of time indices by concatenating random-length blocks
#' sampled circularly from the ordered sequence. Used for creating
#' block-permuted surrogates that preserve short-range temporal structure.
#'
#' @param idx Integer vector of ordered indices.
#' @param block_len Positive integer block length (>= 1).
#'
#' @return Integer vector of permuted indices of the same length as `idx`.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' bioLeak:::.circular_block_permute(1:10, block_len = 3)
#'
#' }
#' @keywords internal
.circular_block_permute <- function(idx, block_len) {
  stopifnot(length(idx) > 0L, is.numeric(block_len), block_len > 0L)
  n <- length(idx)
  block_len <- max(1L, as.integer(round(block_len)))

  out <- integer(0L)
  while (length(out) < n) {
    start <- sample.int(n, 1L)
    block_seq <- ((start - 1L) + seq_len(block_len)) %% n + 1L
    out <- c(out, idx[block_seq])
  }

  out <- out[seq_len(n)]
  stopifnot(length(out) == n)
  out
}


#' Stationary bootstrap indices
#'
#' Implements the stationary bootstrap of Politis & Romano (1994),
#' which resamples contiguous blocks of variable length to preserve
#' weak temporal dependence while maintaining ergodicity.
#'
#' @param idx Integer vector of ordered indices.
#' @param mean_block Positive numeric, expected block length.
#'
#' @return Integer vector of permuted indices of the same length as `idx`.
#'
#' @references
#' Politis, D. N., & Romano, J. P. (1994).
#' *The stationary bootstrap.* Journal of the American Statistical Association, 89(428), 1303-1313.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' bioLeak:::.stationary_bootstrap(1:10, mean_block = 3)
#'
#' }
#' @keywords internal
.stationary_bootstrap <- function(idx, mean_block) {
  stopifnot(length(idx) > 0L, is.numeric(mean_block), mean_block > 0)
  n <- length(idx)
  mean_block <- max(1, mean_block)
  p <- 1 / mean_block

  out <- integer(n)
  pos <- sample.int(n, 1L)
  out[1L] <- idx[pos]

  for (i in 2:n) {
    if (runif(1L) < p) {
      pos <- sample.int(n, 1L)
    } else {
      pos <- pos + 1L
      if (pos > n) pos <- 1L
    }
    out[i] <- idx[pos]
  }

  stopifnot(length(out) == n)
  out
}
