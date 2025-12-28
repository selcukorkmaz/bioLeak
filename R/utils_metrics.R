# Metric utilities -----------------------------------------------------------

.cindex_pairwise <- function(pred, truth) {
  pred <- as.numeric(pred)
  truth <- as.numeric(truth)
  ok <- is.finite(pred) & is.finite(truth)
  pred <- pred[ok]
  truth <- truth[ok]
  n <- length(pred)
  if (n < 2L) return(NA_real_)
  if (length(unique(truth)) < 2L) return(NA_real_)

  conc <- 0L
  ties <- 0L
  total <- 0L

  for (i in seq_len(n - 1L)) {
    yi <- truth[i]
    pi <- pred[i]
    dy <- yi - truth[(i + 1L):n]
    valid <- dy != 0
    if (!any(valid)) next
    pj <- pred[(i + 1L):n][valid]
    dy <- dy[valid]
    dp <- pi - pj
    total <- total + length(dp)
    prod <- dp * dy
    conc <- conc + sum(prod > 0)
    ties <- ties + sum(dp == 0)
  }

  if (!total) return(NA_real_)
  (conc + 0.5 * ties) / total
}
