pkgname <- "bioLeak"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('bioLeak')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("as_rsample")
### * as_rsample

flush(stderr()); flush(stdout())

### Name: as_rsample
### Title: Convert LeakSplits to an rsample resample set
### Aliases: as_rsample

### ** Examples

if (requireNamespace("rsample", quietly = TRUE)) {
  df <- data.frame(
    subject = rep(1:10, each = 2),
    outcome = rbinom(20, 1, 0.5),
    x1 = rnorm(20),
    x2 = rnorm(20)
  )
  splits <- make_split_plan(df, outcome = "outcome",
                        mode = "subject_grouped", group = "subject", v = 5)
  rset <- as_rsample(splits, data = df)
}



cleanEx()
nameEx("audit_leakage")
### * audit_leakage

flush(stderr()); flush(stdout())

### Name: audit_leakage
### Title: Audit leakage and confounding
### Aliases: audit_leakage

### ** Examples

set.seed(1)
df <- data.frame(
  subject = rep(1:6, each = 2),
  outcome = rbinom(12, 1, 0.5),
  x1 = rnorm(12),
  x2 = rnorm(12)
)

splits <- make_split_plan(df, outcome = "outcome",
                      mode = "subject_grouped", group = "subject", v = 3,
                      progress = FALSE)

custom <- list(
  glm = list(
    fit = function(x, y, task, weights, ...) {
      stats::glm(y ~ ., data = as.data.frame(x),
                 family = stats::binomial(), weights = weights)
    },
    predict = function(object, newdata, task, ...) {
      as.numeric(stats::predict(object,
                                newdata = as.data.frame(newdata),
                                type = "response"))
    }
  )
)

fit <- fit_resample(df, outcome = "outcome", splits = splits,
                    learner = "glm", custom_learners = custom,
                    metrics = "auc", refit = FALSE, seed = 1)

audit <- audit_leakage(fit, metric = "auc", B = 10,
                       X_ref = df[, c("x1", "x2")])




cleanEx()
nameEx("audit_leakage_by_learner")
### * audit_leakage_by_learner

flush(stderr()); flush(stdout())

### Name: audit_leakage_by_learner
### Title: Audit leakage per learner
### Aliases: audit_leakage_by_learner

### ** Examples

set.seed(1)
df <- data.frame(
  subject = rep(1:6, each = 2),
  outcome = factor(rep(c(0, 1), 6)),
  x1 = rnorm(12),
  x2 = rnorm(12)
)
splits <- make_split_plan(df, outcome = "outcome",
                      mode = "subject_grouped", group = "subject",
                      v = 3, progress = FALSE)
custom <- list(
  glm = list(
    fit = function(x, y, task, weights, ...) {
      stats::glm(y ~ ., data = data.frame(y = y, x),
                 family = stats::binomial(), weights = weights)
    },
    predict = function(object, newdata, task, ...) {
      as.numeric(stats::predict(object,
                                newdata = as.data.frame(newdata),
                                type = "response"))
    }
  )
)
custom$glm2 <- custom$glm
fit <- fit_resample(df, outcome = "outcome", splits = splits,
                    learner = c("glm", "glm2"), custom_learners = custom,
                    metrics = "auc", refit = FALSE, seed = 1)
audits <- audit_leakage_by_learner(fit, metric = "auc", B = 10,
                                   perm_stratify = FALSE)
names(audits)




cleanEx()
nameEx("audit_report")
### * audit_report

flush(stderr()); flush(stdout())

### Name: audit_report
### Title: Render an HTML audit report
### Aliases: audit_report

### ** Examples

set.seed(1)
df <- data.frame(
  subject = rep(1:6, each = 2),
  outcome = factor(rep(c(0, 1), 6)),
  x1 = rnorm(12),
  x2 = rnorm(12)
)

splits <- make_split_plan(df, outcome = "outcome",
                      mode = "subject_grouped", group = "subject",
                      v = 3, progress = FALSE)

custom <- list(
  glm = list(
    fit = function(x, y, task, weights, ...) {
      stats::glm(y ~ ., data = data.frame(y = y, x),
                 family = stats::binomial(), weights = weights)
    },
    predict = function(object, newdata, task, ...) {
      as.numeric(stats::predict(object,
                                newdata = as.data.frame(newdata),
                                type = "response"))
    }
  )
)

fit <- fit_resample(df, outcome = "outcome", splits = splits,
                    learner = "glm", custom_learners = custom,
                    metrics = "auc", refit = FALSE, seed = 1)

audit <- audit_leakage(fit, metric = "auc", B = 5, perm_stratify = FALSE)

if (requireNamespace("rmarkdown", quietly = TRUE) &&
    requireNamespace("ggplot2", quietly = TRUE) &&
    isTRUE(rmarkdown::pandoc_available("1.12.3"))) {
  out_file <- audit_report(audit, output_dir = tempdir(), quiet = TRUE)
  out_file
}



cleanEx()
nameEx("calibration_summary")
### * calibration_summary

flush(stderr()); flush(stdout())

### Name: calibration_summary
### Title: Calibration diagnostics for binomial predictions
### Aliases: calibration_summary

### ** Examples

set.seed(42)
df <- data.frame(
  subject = rep(1:15, each = 2),
  outcome = factor(rep(c(0, 1), 15)),
  x1 = rnorm(30),
  x2 = rnorm(30)
)
splits <- make_split_plan(df, outcome = "outcome",
                          mode = "subject_grouped", group = "subject",
                          v = 3, progress = FALSE)
custom <- list(
  glm = list(
    fit = function(x, y, task, weights, ...) {
      stats::glm(y ~ ., data = as.data.frame(x),
                 family = stats::binomial(), weights = weights)
    },
    predict = function(object, newdata, task, ...) {
      as.numeric(stats::predict(object, newdata = as.data.frame(newdata),
                                type = "response"))
    }
  )
)
fit <- fit_resample(df, outcome = "outcome", splits = splits,
                    learner = "glm", custom_learners = custom,
                    metrics = "auc", refit = FALSE, seed = 1)
cal <- calibration_summary(fit, bins = 5)
cal$metrics




cleanEx()
nameEx("confounder_sensitivity")
### * confounder_sensitivity

flush(stderr()); flush(stdout())

### Name: confounder_sensitivity
### Title: Confounder sensitivity summaries
### Aliases: confounder_sensitivity

### ** Examples

set.seed(42)
df <- data.frame(
  subject = rep(1:15, each = 2),
  outcome = factor(rep(c(0, 1), 15)),
  batch = factor(rep(c("A", "B", "C"), 10)),
  x1 = rnorm(30),
  x2 = rnorm(30)
)
splits <- make_split_plan(df, outcome = "outcome",
                          mode = "subject_grouped", group = "subject",
                          v = 3, progress = FALSE)
custom <- list(
  glm = list(
    fit = function(x, y, task, weights, ...) {
      stats::glm(y ~ ., data = as.data.frame(x),
                 family = stats::binomial(), weights = weights)
    },
    predict = function(object, newdata, task, ...) {
      as.numeric(stats::predict(object, newdata = as.data.frame(newdata),
                                type = "response"))
    }
  )
)
fit <- fit_resample(df, outcome = "outcome", splits = splits,
                    learner = "glm", custom_learners = custom,
                    metrics = "auc", refit = FALSE, seed = 1)
confounder_sensitivity(fit, confounders = "batch", coldata = df)




cleanEx()
nameEx("dot-guard_fit")
### * dot-guard_fit

flush(stderr()); flush(stdout())

### Name: .guard_fit
### Title: Fit leakage-safe preprocessing pipeline
### Aliases: .guard_fit

### ** Examples

x <- data.frame(a = c(1, 2, NA), b = c(3, 4, 5))
fit <- .guard_fit(x, y = c(1, 2, 3),
                  steps = list(impute = list(method = "median")),
                  task = "gaussian")
fit$transform(x)



cleanEx()
nameEx("fit_resample")
### * fit_resample

flush(stderr()); flush(stdout())

### Name: fit_resample
### Title: Fit and evaluate with leakage guards over predefined splits
### Aliases: fit_resample

### ** Examples

set.seed(1)
df <- data.frame(
  subject = rep(1:10, each = 2),
  outcome = rbinom(20, 1, 0.5),
  x1 = rnorm(20),
  x2 = rnorm(20)
)
splits <- make_split_plan(df, outcome = "outcome",
                      mode = "subject_grouped", group = "subject", v = 5)

# glmnet learner (requires glmnet package)
fit <- fit_resample(df, outcome = "outcome", splits = splits,
                      learner = "glmnet", metrics = "auc")
summary(fit)

# Custom learner (logistic regression) - no extra packages needed
custom <- list(
  glm = list(
    fit = function(x, y, task, weights, ...) {
      stats::glm(y ~ ., data = as.data.frame(x),
                 family = stats::binomial(), weights = weights)
    },
    predict = function(object, newdata, task, ...) {
      as.numeric(stats::predict(object, newdata = as.data.frame(newdata), type = "response"))
    }
  )
)
fit2 <- fit_resample(df, outcome = "outcome", splits = splits,
                     learner = "glm", custom_learners = custom,
                     metrics = "accuracy")

summary(fit2)



cleanEx()
nameEx("impute_guarded")
### * impute_guarded

flush(stderr()); flush(stdout())

### Name: impute_guarded
### Title: Leakage-safe data imputation via guarded preprocessing
### Aliases: impute_guarded

### ** Examples

train <- data.frame(x = c(1, 2, NA, 4), y = c(NA, 1, 1, 0))
test <- data.frame(x = c(NA, 5), y = c(1, NA))
imp <- impute_guarded(train, test, method = "median", winsor = FALSE)
imp$train
imp$test



cleanEx()
nameEx("make_split_plan")
### * make_split_plan

flush(stderr()); flush(stdout())

### Name: make_split_plan
### Title: Create leakage-resistant splits
### Aliases: make_split_plan

### ** Examples

set.seed(1)
df <- data.frame(
  subject = rep(1:10, each = 2),
  outcome = rbinom(20, 1, 0.5),
  x1 = rnorm(20),
  x2 = rnorm(20)
)
splits <- make_split_plan(df, outcome = "outcome",
                      mode = "subject_grouped", group = "subject", v = 5)



cleanEx()
nameEx("plot_calibration")
### * plot_calibration

flush(stderr()); flush(stdout())

### Name: plot_calibration
### Title: Plot calibration curve for binomial predictions
### Aliases: plot_calibration

### ** Examples

if (requireNamespace("ggplot2", quietly = TRUE)) {
  set.seed(42)
  df <- data.frame(
    subject = rep(1:15, each = 2),
    outcome = factor(rep(c(0, 1), 15)),
    x1 = rnorm(30),
    x2 = rnorm(30)
  )
  splits <- make_split_plan(df, outcome = "outcome",
                            mode = "subject_grouped", group = "subject",
                            v = 3, progress = FALSE)
  custom <- list(
    glm = list(
      fit = function(x, y, task, weights, ...) {
        stats::glm(y ~ ., data = as.data.frame(x),
                   family = stats::binomial(), weights = weights)
      },
      predict = function(object, newdata, task, ...) {
        as.numeric(stats::predict(object, newdata = as.data.frame(newdata),
                                  type = "response"))
      }
    )
  )
  fit <- fit_resample(df, outcome = "outcome", splits = splits,
                      learner = "glm", custom_learners = custom,
                      metrics = "auc", refit = FALSE, seed = 1)
  plot_calibration(fit, bins = 5)
}




cleanEx()
nameEx("plot_confounder_sensitivity")
### * plot_confounder_sensitivity

flush(stderr()); flush(stdout())

### Name: plot_confounder_sensitivity
### Title: Plot confounder sensitivity
### Aliases: plot_confounder_sensitivity

### ** Examples

if (requireNamespace("ggplot2", quietly = TRUE)) {
  set.seed(42)
  df <- data.frame(
    subject = rep(1:15, each = 2),
    outcome = factor(rep(c(0, 1), 15)),
    batch = factor(rep(c("A", "B", "C"), 10)),
    x1 = rnorm(30),
    x2 = rnorm(30)
  )
  splits <- make_split_plan(df, outcome = "outcome",
                            mode = "subject_grouped", group = "subject",
                            v = 3, progress = FALSE)
  custom <- list(
    glm = list(
      fit = function(x, y, task, weights, ...) {
        stats::glm(y ~ ., data = as.data.frame(x),
                   family = stats::binomial(), weights = weights)
      },
      predict = function(object, newdata, task, ...) {
        as.numeric(stats::predict(object, newdata = as.data.frame(newdata),
                                  type = "response"))
      }
    )
  )
  fit <- fit_resample(df, outcome = "outcome", splits = splits,
                      learner = "glm", custom_learners = custom,
                      metrics = "auc", refit = FALSE, seed = 1)
  plot_confounder_sensitivity(fit, confounders = "batch", coldata = df)
}




cleanEx()
nameEx("plot_fold_balance")
### * plot_fold_balance

flush(stderr()); flush(stdout())

### Name: plot_fold_balance
### Title: Plot fold balance of class counts per fold
### Aliases: plot_fold_balance

### ** Examples

if (requireNamespace("ggplot2", quietly = TRUE)) {
  set.seed(42)
  df <- data.frame(
    subject = rep(1:15, each = 2),
    outcome = factor(rep(c(0, 1), 15)),
    x1 = rnorm(30),
    x2 = rnorm(30)
  )
  splits <- make_split_plan(df, outcome = "outcome",
                            mode = "subject_grouped", group = "subject",
                            v = 3, progress = FALSE)
  custom <- list(
    glm = list(
      fit = function(x, y, task, weights, ...) {
        stats::glm(y ~ ., data = as.data.frame(x),
                   family = stats::binomial(), weights = weights)
      },
      predict = function(object, newdata, task, ...) {
        as.numeric(stats::predict(object, newdata = as.data.frame(newdata),
                                  type = "response"))
      }
    )
  )
  fit <- fit_resample(df, outcome = "outcome", splits = splits,
                      learner = "glm", custom_learners = custom,
                      metrics = "auc", refit = FALSE, seed = 1)
  plot_fold_balance(fit)
}




cleanEx()
nameEx("plot_overlap_checks")
### * plot_overlap_checks

flush(stderr()); flush(stdout())

### Name: plot_overlap_checks
### Title: Plot overlap diagnostics between train/test groups
### Aliases: plot_overlap_checks

### ** Examples

set.seed(1)
df <- data.frame(
  subject = rep(1:6, each = 2),
  outcome = rbinom(12, 1, 0.5),
  x1 = rnorm(12),
  x2 = rnorm(12)
)
splits <- make_split_plan(df, outcome = "outcome",
                      mode = "subject_grouped", group = "subject", v = 3)
custom <- list(
  glm = list(
    fit = function(x, y, task, weights, ...) {
      stats::glm(y ~ ., data = as.data.frame(x),
                 family = stats::binomial(), weights = weights)
    },
    predict = function(object, newdata, task, ...) {
      as.numeric(stats::predict(object, newdata = as.data.frame(newdata),
                                type = "response"))
    }
  )
)
fit <- fit_resample(df, outcome = "outcome", splits = splits,
                    learner = "glm", custom_learners = custom,
                    metrics = "accuracy", refit = FALSE)
if (requireNamespace("ggplot2", quietly = TRUE)) {
  out <- plot_overlap_checks(fit, column = "subject")
  out$overlap_counts
}




cleanEx()
nameEx("plot_perm_distribution")
### * plot_perm_distribution

flush(stderr()); flush(stdout())

### Name: plot_perm_distribution
### Title: Plot permutation distribution for a LeakAudit object
### Aliases: plot_perm_distribution

### ** Examples

if (requireNamespace("ggplot2", quietly = TRUE)) {
  set.seed(42)
  df <- data.frame(
    subject = rep(1:15, each = 2),
    outcome = factor(rep(c(0, 1), 15)),
    x1 = rnorm(30),
    x2 = rnorm(30)
  )
  splits <- make_split_plan(df, outcome = "outcome",
                            mode = "subject_grouped", group = "subject",
                            v = 3, progress = FALSE)
  custom <- list(
    glm = list(
      fit = function(x, y, task, weights, ...) {
        stats::glm(y ~ ., data = as.data.frame(x),
                   family = stats::binomial(), weights = weights)
      },
      predict = function(object, newdata, task, ...) {
        as.numeric(stats::predict(object, newdata = as.data.frame(newdata),
                                  type = "response"))
      }
    )
  )
  fit <- fit_resample(df, outcome = "outcome", splits = splits,
                      learner = "glm", custom_learners = custom,
                      metrics = "auc", refit = FALSE, seed = 1)
  audit <- audit_leakage(fit, metric = "auc", B = 20)
  plot_perm_distribution(audit)
}




cleanEx()
nameEx("plot_time_acf")
### * plot_time_acf

flush(stderr()); flush(stdout())

### Name: plot_time_acf
### Title: Plot ACF of test predictions for time-series leakage checks
### Aliases: plot_time_acf

### ** Examples

if (requireNamespace("ggplot2", quietly = TRUE)) {
  set.seed(42)
  df <- data.frame(
    id = 1:30,
    time = seq.Date(as.Date("2020-01-01"), by = "day", length.out = 30),
    y = rnorm(30),
    x1 = rnorm(30),
    x2 = rnorm(30)
  )
  splits <- make_split_plan(df, outcome = "y", mode = "time_series",
                            time = "time", v = 3, progress = FALSE)
  custom <- list(
    lm = list(
      fit = function(x, y, task, weights, ...) {
        stats::lm(y ~ ., data = data.frame(y = y, x))
      },
      predict = function(object, newdata, task, ...) {
        as.numeric(stats::predict(object, newdata = as.data.frame(newdata)))
      }
    )
  )
  fit <- fit_resample(df, outcome = "y", splits = splits,
                      learner = "lm", custom_learners = custom,
                      metrics = "rmse", refit = FALSE, seed = 1)
  plot_time_acf(fit, lag.max = 10)
}




cleanEx()
nameEx("predict_guard")
### * predict_guard

flush(stderr()); flush(stdout())

### Name: predict_guard
### Title: Apply a fitted GuardFit transformer to new data
### Aliases: predict_guard

### ** Examples

x_train <- data.frame(a = c(1, 2, NA, 4), b = c(10, 11, 12, 13))
fit <- .guard_fit(
  x_train,
  y = c(0.1, 0.2, 0.3, 0.4),
  steps = list(impute = list(method = "median")),
  task = "gaussian"
)
x_new <- data.frame(a = c(NA, 5), b = c(9, 14))
out <- predict_guard(fit, x_new)
out



cleanEx()
nameEx("show-LeakSplits-method")
### * show-LeakSplits-method

flush(stderr()); flush(stdout())

### Name: show,LeakSplits-method
### Title: Display summary for LeakSplits objects
### Aliases: show,LeakSplits-method

### ** Examples

df <- data.frame(
  subject = rep(1:10, each = 2),
  outcome = rbinom(20, 1, 0.5),
  x1 = rnorm(20),
  x2 = rnorm(20)
)
splits <- make_split_plan(df, outcome = "outcome",
                      mode = "subject_grouped", group = "subject", v = 5)
show(splits)



cleanEx()
nameEx("simulate_leakage_suite")
### * simulate_leakage_suite

flush(stderr()); flush(stdout())

### Name: simulate_leakage_suite
### Title: Simulate leakage scenarios and audit results
### Aliases: simulate_leakage_suite

### ** Examples




cleanEx()
nameEx("summary.LeakAudit")
### * summary.LeakAudit

flush(stderr()); flush(stdout())

### Name: summary.LeakAudit
### Title: Summarize a leakage audit
### Aliases: summary.LeakAudit

### ** Examples

set.seed(1)
df <- data.frame(
  subject = rep(1:6, each = 2),
  outcome = rbinom(12, 1, 0.5),
  x1 = rnorm(12),
  x2 = rnorm(12)
)
splits <- make_split_plan(df, outcome = "outcome",
                      mode = "subject_grouped", group = "subject", v = 3)
custom <- list(
  glm = list(
    fit = function(x, y, task, weights, ...) {
      stats::glm(y ~ ., data = as.data.frame(x),
                 family = stats::binomial(), weights = weights)
    },
    predict = function(object, newdata, task, ...) {
      as.numeric(stats::predict(object, newdata = as.data.frame(newdata),
                                type = "response"))
    }
  )
)
fit <- fit_resample(df, outcome = "outcome", splits = splits,
                    learner = "glm", custom_learners = custom,
                    metrics = "auc", refit = FALSE, seed = 1)
audit <- audit_leakage(fit, metric = "auc", B = 5,
                       X_ref = df[, c("x1", "x2")], seed = 1)
summary(audit) # prints the audit report and returns `audit` invisibly




cleanEx()
nameEx("summary.LeakFit")
### * summary.LeakFit

flush(stderr()); flush(stdout())

### Name: summary.LeakFit
### Title: Summarize a LeakFit object
### Aliases: summary.LeakFit

### ** Examples

set.seed(1)
df <- data.frame(
  subject = rep(1:6, each = 2),
  outcome = factor(rep(c(0, 1), each = 6)),
  x1 = rnorm(12),
  x2 = rnorm(12)
)
splits <- make_split_plan(
  df,
  outcome = "outcome",
  mode = "subject_grouped",
  group = "subject",
  v = 3,
  stratify = TRUE,
  progress = FALSE
)
custom <- list(
  glm = list(
    fit = function(x, y, task, weights, ...) {
      stats::glm(y ~ ., data = data.frame(y = y, x),
                 family = stats::binomial(), weights = weights)
    },
    predict = function(object, newdata, task, ...) {
      as.numeric(stats::predict(object,
                                newdata = as.data.frame(newdata),
                                type = "response"))
    }
  )
)
fit <- fit_resample(df, outcome = "outcome", splits = splits,
                    learner = "glm", custom_learners = custom,
                    metrics = "auc", seed = 1)
summary_df <- summary(fit)
summary_df




cleanEx()
nameEx("tune_resample")
### * tune_resample

flush(stderr()); flush(stdout())

### Name: tune_resample
### Title: Leakage-aware nested tuning with tidymodels
### Aliases: tune_resample

### ** Examples





### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
