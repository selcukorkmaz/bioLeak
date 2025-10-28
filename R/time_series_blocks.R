# Time-series block permutation utilities

#' Circular block permutation indices
#'
#' @param idx integer vector of ordered indices
#' @param block_len positive integer block length
#' @return integer vector of permuted indices of same length as `idx`
.circular_block_permute <- function(idx, block_len) {
  stopifnot(length(idx) > 0L, block_len > 0L)
  n <- length(idx)
  block_len <- max(1L, as.integer(block_len))
  out <- integer(0L)
  while (length(out) < n) {
    start <- sample.int(n, 1L)
    block_seq <- ((start - 1L) + seq_len(block_len)) %% n + 1L
    out <- c(out, idx[block_seq])
  }
  out[seq_len(n)]
}

#' Stationary bootstrap indices
#'
#' @param idx integer vector of ordered indices
#' @param mean_block positive numeric, expected block length
#' @return integer vector of permuted indices
.stationary_bootstrap <- function(idx, mean_block) {
  stopifnot(length(idx) > 0L, mean_block > 0)
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
  out
}
