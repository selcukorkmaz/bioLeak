# Train-only and no-silent-no-op invariants for guard_fit() preprocessing.
# Run from the bioLeak package root:
#   Rscript .claude/skills/verify-change/probes/guard_invariants.R

source(".claude/skills/verify-change/probes/_probe_helpers.R")
load_dev(".")

set.seed(3)
counts <- function(n, p = 6) as.data.frame(matrix(rpois(n * p, 50), n, p,
                                                  dimnames = list(NULL, paste0("g", seq_len(p)))))
Xtr <- counts(40); Xte <- counts(10)
ytr <- factor(rep(0:1, 20))

# --- Row independence: transform(row | other test rows) == transform(row alone).
for (m in c("zscore", "robust", "quantile", "log_cpm", "clr")) {
  r <- capture({
    gf <- guard_fit(Xtr, ytr, steps = list(normalize = list(method = m)), task = "binomial")
    all_rows <- gf$transform(Xte)
    one_row  <- gf$transform(Xte[1, , drop = FALSE])
    isTRUE(all.equal(unname(unlist(all_rows[1, ])), unname(unlist(one_row[1, ]))))
  })
  check("INV-2", is.null(r$error) && isTRUE(r$value),
        sprintf("normalize='%s': test row transformed alone == with other test rows%s", m,
                if (!is.null(r$error)) paste0(" (error: ", r$error, ")") else ""))
}

# --- State hash must change when configuration changes.
h1 <- guard_fit(Xtr, ytr, steps = list(normalize = list(method = "zscore")), task = "binomial")$state$hash
h2 <- guard_fit(Xtr, ytr, steps = list(normalize = list(method = "quantile")), task = "binomial")$state$hash
check("INV-3", !is.null(h1) && !identical(h1, h2), "state$hash differs between zscore and quantile configs")

# --- BL-3: unseen batch at transform time must warn or error, never pass silently.
set.seed(2)
Xb  <- data.frame(a = c(rnorm(20, 0), rnorm(20, 3)), b = rnorm(40))
btr <- rep(c("A", "B"), each = 20)
gfb <- guard_fit(Xb, factor(rep(0:1, 20)),
                 steps = list(batch_adjust = list(method = "combat_frozen", batch = btr),
                              normalize = list(method = "none")),
                 task = "binomial")
Xnew <- data.frame(a = rnorm(10, 10), b = rnorm(10))
r <- capture(gfb$transform(Xnew, batch = rep("C", 10)))
check("BL-3", length(r$warnings) > 0L || !is.null(r$error),
      sprintf("unseen batch 'C' at transform: %s",
              if (!is.null(r$error)) paste("error:", r$error)
              else if (length(r$warnings)) paste("warning:", r$warnings[1])
              else "no condition raised (silent no-op)"))

# --- BL-4: batch confounded with outcome, no true batch effect.
# Adjustment must keep the biological effect, or warn about the confounding.
set.seed(3)
y  <- rep(c(0, 1), each = 40)
bt <- ifelse(runif(80) < ifelse(y == 1, 0.85, 0.15), "B1", "B2")
Xc <- data.frame(sig = y * 2 + rnorm(80))
r <- capture({
  gfc <- guard_fit(Xc, factor(y),
                   steps = list(batch_adjust = list(method = "location_scale", batch = bt),
                                normalize = list(method = "none")),
                   task = "binomial")
  gfc$transform(Xc, batch = bt)
})
before <- diff(tapply(Xc$sig, y, mean))
after  <- if (is.null(r$error)) diff(tapply(r$value$sig, y, mean)) else NA_real_
check("BL-4", (!is.na(after) && after / before >= 0.9) || length(r$warnings) > 0L,
      sprintf("confounded batch: effect %.2f -> %.2f (retained %.0f%%), warnings: %d",
              before, after, 100 * after / before, length(r$warnings)))

probe_finish()
