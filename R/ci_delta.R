# Gap statistic Δ = observed - mean(perm) standard errors ----------------------

.se_ci_delta <- function(delta, se_obs, perm_values, level = 0.95) {
  perm_values <- perm_values[is.finite(perm_values)]
  if (!length(perm_values)) {
    return(list(se = NA_real_, ci = c(NA_real_, NA_real_), z = NA_real_))
  }
  var_perm <- stats::var(perm_values)
  se_perm <- sqrt(var_perm / length(perm_values))
  se_obs <- ifelse(is.finite(se_obs), se_obs, 0)
  se_delta <- sqrt(se_obs^2 + se_perm^2)
  z <- if (se_delta > 0) delta / se_delta else NA_real_
  alpha <- (1 - level) / 2
  q <- stats::qnorm(c(alpha, 1 - alpha))
  ci <- delta + q * se_delta
  list(se = se_delta, ci = ci, z = z)
}
