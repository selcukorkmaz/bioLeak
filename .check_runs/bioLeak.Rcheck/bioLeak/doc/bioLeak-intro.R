## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  message = FALSE,
  warning = FALSE,
  eval = TRUE
)

## ----example-data-------------------------------------------------------------
library(bioLeak)

set.seed(123)
n <- 160
subject <- rep(seq_len(40), each = 4)
batch <- sample(paste0("B", 1:6), n, replace = TRUE)
study <- sample(paste0("S", 1:4), n, replace = TRUE)
time <- seq_len(n)

x1 <- rnorm(n)
x2 <- rnorm(n)
x3 <- rnorm(n)
linpred <- 0.7 * x1 - 0.4 * x2 + 0.2 * x3 + rnorm(n, sd = 0.5)
p <- stats::plogis(linpred)
outcome <- factor(ifelse(runif(n) < p, "case", "control"),
                  levels = c("control", "case"))

df <- data.frame(
  subject = subject,
  batch = batch,
  study = study,
  time = time,
  outcome = outcome,
  x1 = x1,
  x2 = x2,
  x3 = x3
)

df_leaky <- within(df, {
  leak_subject <- ave(as.numeric(outcome == "case"), subject, FUN = mean)
  leak_batch <- ave(as.numeric(outcome == "case"), batch, FUN = mean)
  leak_global <- mean(as.numeric(outcome == "case"))
})

df_time <- df
df_time$leak_future <- c(as.numeric(df_time$outcome == "case")[-1], 0)
predictors <- c("x1", "x2", "x3")


# Example data (first 6 rows)
head(df)

# Outcome class counts
as.data.frame(table(df$outcome))

## ----leaky-splits-------------------------------------------------------------
leaky_splits <- make_split_plan(
  df,
  outcome = "outcome",
  mode = "subject_grouped",
  group = "row_id",
  v = 5,
  repeats = 1,
  stratify = TRUE
)

cat("Leaky splits summary (sample-wise CV):\n")
leaky_splits

## ----safe-splits--------------------------------------------------------------
safe_splits <- make_split_plan(
  df,
  outcome = "outcome",
  mode = "subject_grouped",
  group = "subject",
  v = 5,
  repeats = 1,
  stratify = TRUE,
  seed = 10
)

cat("Leakage-safe splits summary (subject-grouped CV):\n")
safe_splits

## ----splits-other-modes-------------------------------------------------------
batch_splits <- make_split_plan(
  df,
  outcome = "outcome",
  mode = "batch_blocked",
  batch = "batch",
  v = 4,
  stratify = TRUE
)

cat("Batch-blocked splits summary:\n")
batch_splits

study_splits <- make_split_plan(
  df,
  outcome = "outcome",
  mode = "study_loocv",
  study = "study"
)

cat("Study leave-one-out splits summary:\n")
study_splits

time_splits <- make_split_plan(
  df,
  outcome = "outcome",
  mode = "time_series",
  time = "time",
  v = 4,
  horizon = 2
)

cat("Time-series splits summary:\n")
time_splits

nested_splits <- make_split_plan(
  df,
  outcome = "outcome",
  mode = "subject_grouped",
  group = "subject",
  v = 3,
  nested = TRUE,
  stratify = TRUE
)

cat("Nested CV splits summary:\n")
nested_splits

## ----combined-splits----------------------------------------------------------
# For combined mode, constraint-axis levels should not span all folds.
# Here each subject belongs to exactly one site, so site exclusion only
# removes training samples from the same site as the test subjects.
df_comb <- df
df_comb$site <- paste0("site", rep(1:8, each = 5)[df_comb$subject])

combined_splits <- make_split_plan(
  df_comb,
  outcome = "outcome",
  mode = "combined",
  constraints = list(
    list(type = "subject", col = "subject"),
    list(type = "batch",   col = "site")
  ),
  v = 4,
  stratify = TRUE,
  seed = 42
)

cat("Combined (N-axis) splits summary:\n")
combined_splits

## ----combined-overlap-check---------------------------------------------------
# Verify zero overlap on all constraint axes
overlap_result <- check_split_overlap(combined_splits)
overlap_result

## ----check-split-overlap------------------------------------------------------
# Run on the subject-grouped splits from earlier
overlap_safe <- check_split_overlap(safe_splits)
overlap_safe

## ----compact-splits-----------------------------------------------------------
# Efficient storage for large N
big_splits <- make_split_plan(
  df,
  outcome = "outcome",
  mode = "subject_grouped",
  group = "subject",
  v = 5,
  compact = TRUE  # <--- Saves memory
)

cat("Compact-mode splits summary:\n")
big_splits
cat(sprintf("Compact storage enabled: %s\n", big_splits@info$compact))

## ----strict-mode-demo, error=TRUE---------------------------------------------
try({
# Enable strict mode temporarily
withr::with_options(list(bioLeak.strict = TRUE), {
  strict_splits <- make_split_plan(
    df,
    outcome = "outcome",
    mode = "subject_grouped",
    group = "subject",
    v = 5,
    stratify = TRUE,
    seed = 42
  )
  cat("Strict mode splits completed without error.\n")
})
})

## ----strict-mode-trained-recipe, error=TRUE-----------------------------------
try({
# Strict mode catches a pre-trained recipe
if (requireNamespace("recipes", quietly = TRUE)) {
  rec <- recipes::recipe(outcome ~ ., data = df[, c("outcome", predictors)]) |>
    recipes::step_normalize(recipes::all_numeric_predictors()) |>
    recipes::prep(training = df[, c("outcome", predictors)])

  withr::with_options(list(bioLeak.strict = TRUE), {
    tryCatch(
      fit_resample(
        df, outcome = "outcome",
        splits = safe_splits,
        learner = "glmnet",
        preprocess = rec
      ),
      error = function(e) cat("Strict mode error:", conditionMessage(e), "\n")
    )
  })
}
})

## ----leaky-scaling------------------------------------------------------------
df_leaky_scaled <- df
df_leaky_scaled[predictors] <- scale(df_leaky_scaled[predictors])
scaled_summary <- data.frame(
  feature = predictors,
  mean = colMeans(df_leaky_scaled[predictors]),
  sd = apply(df_leaky_scaled[predictors], 2, stats::sd)
)
scaled_summary$mean <- round(scaled_summary$mean, 3)
scaled_summary$sd <- round(scaled_summary$sd, 3)

# Leaky global scaling: means ~0 and SDs ~1 computed on all samples
scaled_summary

## ----guarded-preprocess-------------------------------------------------------
fold1 <- safe_splits@indices[[1]]
train_x <- df[fold1$train, predictors]
test_x <- df[fold1$test, predictors]

guard <- .guard_fit(
  X = train_x,
  y = df$outcome[fold1$train],
  steps = list(
    impute = list(method = "median", winsor = TRUE),
    normalize = list(method = "zscore"),
    filter = list(var_thresh = 0, iqr_thresh = 0),
    fs = list(method = "none")
  ),
  task = "binomial"
)

train_x_guarded <- predict_guard(guard, train_x)
test_x_guarded <- predict_guard(guard, test_x)

cat("GuardFit object:\n")
guard
cat("\nGuardFit summary:\n")
summary(guard)

# Guarded training data (first 6 rows)
head(train_x_guarded)

# Guarded test data (first 6 rows)
head(test_x_guarded)

## ----guard-ensure-levels------------------------------------------------------
raw_levels <- data.frame(
  site = c("A", "B", "B"),
  status = c("yes", "no", "yes"),
  stringsAsFactors = FALSE
)

level_state <- .guard_ensure_levels(raw_levels)

# Aligned factor data with consistent levels
level_state$data

# Levels map
level_state$levels

## ----leaky-impute-------------------------------------------------------------
train <- data.frame(a = c(1, 2, NA, 4), b = c(NA, 1, 1, 0))
test <- data.frame(a = c(NA, 5), b = c(1, NA))

all_median <- vapply(rbind(train, test),
                     function(col) median(col, na.rm = TRUE),
                     numeric(1))
train_leaky <- as.data.frame(Map(function(col, m) { col[is.na(col)] <- m; col },
                                 train, all_median))
test_leaky <- as.data.frame(Map(function(col, m) { col[is.na(col)] <- m; col },
                                test, all_median))

# Leaky medians computed on train + test
data.frame(feature = names(all_median), median = all_median)

# Leaky-imputed training data
train_leaky

# Leaky-imputed test data
test_leaky

## ----safe-impute--------------------------------------------------------------
imp <- impute_guarded(
  train = train,
  test = test,
  method = "median",
  winsor = FALSE
)

# Guarded-imputed training data
imp$train

# Guarded-imputed test data
imp$test

## ----guard-to-recipe----------------------------------------------------------
if (requireNamespace("recipes", quietly = TRUE)) {
  guard_steps <- list(
    impute    = list(method = "median"),
    normalize = list(method = "zscore")
  )

  rec <- guard_to_recipe(
    steps = guard_steps,
    formula = outcome ~ .,
    training_data = df[, c("outcome", predictors)]
  )

  cat("Converted recipe:\n")
  rec
}

## ----parsnip-spec-------------------------------------------------------------
spec <- parsnip::logistic_reg(mode = "classification") |>
  parsnip::set_engine("glm")

## ----fit-leaky----------------------------------------------------------------
fit_leaky <- fit_resample(
  df_leaky,
  outcome = "outcome",
  splits = leaky_splits,
  learner = spec,
  metrics = c("auc", "accuracy"),
  preprocess = list(
    impute = list(method = "median"),
    normalize = list(method = "zscore"),
    filter = list(var_thresh = 0),
    fs = list(method = "none")
  )
)

cat("Leaky fit summary:\n")
summary(fit_leaky)
metrics_leaky <- as.data.frame(fit_leaky@metric_summary)
num_cols <- vapply(metrics_leaky, is.numeric, logical(1))
metrics_leaky[num_cols] <- lapply(metrics_leaky[num_cols], round, digits = 3)

# Leaky fit: mean and SD of metrics across folds
metrics_leaky

## ----fit-safe-----------------------------------------------------------------
fit_safe <- fit_resample(
  df,
  outcome = "outcome",
  splits = safe_splits,
  learner = spec,
  metrics = c("auc", "accuracy"),
  preprocess = list(
    impute = list(method = "median"),
    normalize = list(method = "zscore"),
    filter = list(var_thresh = 0),
    fs = list(method = "none")
  ),
  positive_class = "case",
  class_weights = c(control = 1, case = 1),
  refit = TRUE
)

cat("Leakage-safe fit summary:\n")
summary(fit_safe)
metrics_safe <- as.data.frame(fit_safe@metric_summary)
num_cols <- vapply(metrics_safe, is.numeric, logical(1))
metrics_safe[num_cols] <- lapply(metrics_safe[num_cols], round, digits = 3)

# Leakage-safe fit: mean and SD of metrics across folds
metrics_safe

# Per-fold metrics (first 6 rows)
head(fit_safe@metrics)

## ----fit-fold-status----------------------------------------------------------
# Fold-level diagnostics for reproducible troubleshooting
head(fit_safe@info$fold_status, 5)

## ----multiclass-example-------------------------------------------------------
if (requireNamespace("ranger", quietly = TRUE)) {
  set.seed(11)
  df_multi <- df
  df_multi$outcome3 <- factor(sample(c("A", "B", "C"),
                                     nrow(df_multi), replace = TRUE))

  multi_splits <- make_split_plan(
    df_multi,
    outcome = "outcome3",
    mode = "subject_grouped",
    group = "subject",
    v = 4,
    stratify = TRUE,
    seed = 11
  )

  fit_multi <- fit_resample(
    df_multi,
    outcome = "outcome3",
    splits = multi_splits,
    learner = "ranger",
    metrics = c("accuracy", "macro_f1", "log_loss"),
    refit = FALSE
  )

  cat("Multiclass fit summary:\n")
  summary(fit_multi)
} else {
  cat("ranger not installed; skipping multiclass example.\n")
}

## ----survival-example---------------------------------------------------------
if (requireNamespace("survival", quietly = TRUE) &&
    requireNamespace("glmnet", quietly = TRUE)) {
  set.seed(12)
  df_surv <- df
  df_surv$time_to_event <- rexp(nrow(df_surv), rate = 0.1)
  df_surv$event <- rbinom(nrow(df_surv), 1, 0.7)

  surv_splits <- make_split_plan(
    df_surv,
    outcome = c("time_to_event", "event"),
    mode = "subject_grouped",
    group = "subject",
    v = 4,
    stratify = FALSE,
    seed = 12
  )

  fit_surv <- fit_resample(
    df_surv,
    outcome = c("time_to_event", "event"),
    splits = surv_splits,
    learner = "glmnet",
    metrics = "cindex",
    refit = FALSE
  )

  cat("Survival fit summary:\n")
  summary(fit_surv)
} else {
  cat("survival or glmnet not installed; skipping survival example.\n")
}

## ----learner-args-------------------------------------------------------------
if (requireNamespace("glmnet", quietly = TRUE)) {
  fit_glmnet <- fit_resample(
    df,
    outcome = "outcome",
    splits = safe_splits,
    learner = "glmnet",
    metrics = "auc",
    learner_args = list(glmnet = list(alpha = 0.5)),
    preprocess = list(
      impute = list(method = "median"),
      normalize = list(method = "zscore"),
      filter = list(var_thresh = 0),
      fs = list(method = "none")
    )
  )
  cat("GLMNET summary with learner-specific arguments:\n")
  summary(fit_glmnet)
} else {
  cat("glmnet not installed; skipping learner_args example.\n")
}

## ----se-example---------------------------------------------------------------
if (requireNamespace("SummarizedExperiment", quietly = TRUE)) {
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = t(as.matrix(df[, predictors]))),
    colData = df[, c("subject", "batch", "study", "time", "outcome"), drop = FALSE]
  )

  se_splits <- make_split_plan(
    se,
    outcome = "outcome",
    mode = "subject_grouped",
    group = "subject",
    v = 3
  )

  se_fit <- fit_resample(
    se,
    outcome = "outcome",
    splits = se_splits,
    learner = spec,
    metrics = "auc"
  )
  cat("SummarizedExperiment fit summary:\n")
  summary(se_fit)
} else {
  cat("SummarizedExperiment not installed; skipping SE example.\n")
}

## ----tidymodels-interop-------------------------------------------------------
library(bioLeak)
library(parsnip)
library(recipes)
library(yardstick)

set.seed(123)
N <- 60
df <- data.frame(
  subject = factor(rep(paste0("S", 1:20), length.out = N)), 
  outcome = factor(rep(c("ClassA", "ClassB"), length.out = N)),
  x1 = rnorm(N),
  x2 = rnorm(N),
  x3 = rnorm(N)
)

spec <- logistic_reg() |> set_engine("glm")

# Use bioLeak's native split planner to avoid conversion errors
set.seed(13)

# Use make_split_plan instead of rsample::group_vfold_cv
# This creates a subject-grouped CV directly compatible with fit_resample
splits <- make_split_plan(
  df, 
  outcome = "outcome", 
  mode = "subject_grouped", 
  group = "subject", 
  v = 3
)

rec <- recipes::recipe(outcome ~ x1 + x2 + x3, data = df) |>
  recipes::step_impute_median(recipes::all_numeric_predictors()) |>
  recipes::step_normalize(recipes::all_numeric_predictors())

metrics_set <- yardstick::metric_set(yardstick::roc_auc, yardstick::accuracy)

fit_rs <- fit_resample(
  df,
  outcome = "outcome",
  splits = splits,
  learner = spec,
  preprocess = rec,
  metrics = metrics_set,
  refit = FALSE
)

if (exists("as_rsample", where = asNamespace("bioLeak"), mode = "function")) {
    rs_export <- as_rsample(fit_rs@splits, data = df)
    print(rs_export)
}

## ----workflow-example---------------------------------------------------------
if (requireNamespace("workflows", quietly = TRUE)) {
  wf <- workflows::workflow() |>
    workflows::add_model(spec) |>
    workflows::add_formula(outcome ~ x1 + x2 + x3)

  fit_wf <- fit_resample(
    df,
    outcome = "outcome",
    splits = splits,
    learner = wf,
    metrics = "auc",
    refit = FALSE
  )

  cat("Workflow fit summary:\n")
  summary(fit_wf)
} else {
  cat("workflows not installed; skipping workflow example.\n")
}

## ----tune-example-------------------------------------------------------------
if (requireNamespace("parsnip", quietly = TRUE) &&
    requireNamespace("recipes", quietly = TRUE) &&
    requireNamespace("tune", quietly = TRUE) &&
    requireNamespace("glmnet", quietly = TRUE)) {
  # --- 1. Create Data ---
  set.seed(123)
  N <- 60
  df <- data.frame(
    subject = factor(rep(paste0("S", 1:20), length.out = N)),
    outcome = factor(sample(c("ClassA", "ClassB"), N, replace = TRUE)),
    x1 = rnorm(N),
    x2 = rnorm(N),
    x3 = rnorm(N)
  )

  # --- 2. Generate Nested Splits ---
  set.seed(1)
  nested_splits <- make_split_plan(
    df,
    outcome = "outcome",
    mode = "subject_grouped",
    group = "subject",
    v = 3,
    nested = TRUE
  )

  # --- 3. Define Recipe & Model ---
  rec <- recipes::recipe(outcome ~ x1 + x2 + x3, data = df) |>
    recipes::step_impute_median(recipes::all_numeric_predictors()) |>
    recipes::step_normalize(recipes::all_numeric_predictors())

  spec_tune <- parsnip::logistic_reg(
    penalty = tune::tune(),
    mixture = 1,
    mode = "classification"
  ) |>
    parsnip::set_engine("glmnet")

  # --- 4. Run Tuning ---
  tuned <- tune_resample(
    df,
    outcome = "outcome",
    splits = nested_splits,
    learner = spec_tune,
    preprocess = rec,
    inner_v = 2,
    grid = 3,
    metrics = c("auc", "accuracy"),
    selection = "one_std_err",
    refit = TRUE,
    seed = 14
  )

  summary(tuned)
} else {
  cat("parsnip/recipes/tune/glmnet not installed; skipping nested tuning example.\n")
}

## ----tune-diagnostics---------------------------------------------------------
if (exists("tuned")) {
  # Fold-level status from the outer loop
  head(tuned$fold_status, 5)

  # Final tuned model is available when refit = TRUE
  is.null(tuned$final_model)
} else {
  cat("Nested tuning object not available (dependencies missing).\n")
}

## ----fit-xgboost--------------------------------------------------------------
if (requireNamespace("parsnip", quietly = TRUE) &&
    requireNamespace("xgboost", quietly = TRUE) &&
    requireNamespace("recipes", quietly = TRUE)) {
  
  # 1. Define the model spec
  xgb_spec <- parsnip::boost_tree(
    mode = "classification",
    trees = 100,
    tree_depth = 6,
    learn_rate = 0.01
  ) |>
    parsnip::set_engine("xgboost")
  
  # 2. Define a recipe on data without 'subject' because split metadata
  # columns are excluded from predictors by fit_resample().
  df_for_rec <- df[, !names(df) %in% "subject"]
  
  rec_xgb <- recipes::recipe(outcome ~ ., data = df_for_rec) |>
    recipes::step_dummy(recipes::all_nominal_predictors()) |>
    recipes::step_impute_median(recipes::all_numeric_predictors())
  # Note: No need for step_rm(subject) because it's already gone!
  
  # 3. Fit
  fit_xgb <- fit_resample(
    df,
    outcome = "outcome",
    splits = splits,
    learner = xgb_spec,
    metrics = "auc",
    preprocess = rec_xgb 
  )
  
  cat("XGBoost parsnip fit summary:\n")
  print(summary(fit_xgb))
} else {
  cat("parsnip/xgboost/recipes not installed.\n")
}

## ----custom-learners----------------------------------------------------------
custom_learners <- list(
  glm = list(
    fit = function(x, y, task, weights, ...) {
      df_fit <- data.frame(y = y, x, check.names = FALSE)
      stats::glm(y ~ ., data = df_fit,
                 family = stats::binomial(), weights = weights)
    },
    predict = function(object, newdata, task, ...) {
      as.numeric(stats::predict(object, newdata = as.data.frame(newdata),
                                type = "response"))
    }
  )
)

cat("Custom learner names:\n")
names(custom_learners)
cat("Custom learner components (fit/predict):\n")
lapply(custom_learners, names)

## ----plot-fold-balance--------------------------------------------------------
if (requireNamespace("ggplot2", quietly = TRUE)) {
  plot_fold_balance(fit_safe)
} else {
  cat("ggplot2 not installed; skipping fold balance plot.\n")
}

## ----plot-calibration---------------------------------------------------------
if (requireNamespace("ggplot2", quietly = TRUE)) {
  plot_calibration(fit_safe, bins = 10)
} else {
  cat("ggplot2 not installed; skipping calibration plot.\n")
}

## ----plot-confounder-sensitivity----------------------------------------------
if (requireNamespace("ggplot2", quietly = TRUE)) {
  plot_confounder_sensitivity(fit_safe, confounders = c("batch", "study"),
                              metric = "auc", min_n = 3)
} else {
  cat("ggplot2 not installed; skipping confounder sensitivity plot.\n")
}

## ----diagnostics-tables-------------------------------------------------------
cal <- calibration_summary(fit_safe, bins = 10, min_bin_n = 5)
conf_tbl <- confounder_sensitivity(fit_safe,
                                   confounders = c("batch", "study"),
                                   metric = "auc",
                                   min_n = 3)

# Calibration metrics
cal$metrics

# Confounder sensitivity table (first 6 rows)
head(conf_tbl)

## ----plot-overlap-------------------------------------------------------------
if (requireNamespace("ggplot2", quietly = TRUE)) {
  plot_overlap_checks(fit_leaky, column = "subject")
  plot_overlap_checks(fit_safe, column = "subject")
} else {
  cat("ggplot2 not installed; skipping overlap plots.\n")
}

## ----audit-leakage------------------------------------------------------------
# Use df_leaky (the original 160-row data) so X_ref aligns with fit_safe's predictions.
# df may have been redefined to a smaller dataset in the tidymodels section above.
X_ref <- df_leaky[, predictors]
X_ref[c(1, 5), ] <- X_ref[1, ]

audit <- audit_leakage(
  fit_safe,
  metric = "auc",
  B = 20,
  perm_stratify = TRUE,
  batch_cols = c("batch", "study"),
  X_ref = X_ref,
  sim_method = "cosine",
  sim_threshold = 0.995,
  return_perm = TRUE
)

cat("Leakage audit summary:\n")
summary(audit)
if (!is.null(audit@permutation_gap) && nrow(audit@permutation_gap) > 0) {
  # Permutation significance results
  audit@permutation_gap
}
if (!is.null(audit@batch_assoc) && nrow(audit@batch_assoc) > 0) {
  # Batch/study association with folds (Cramer's V)
  audit@batch_assoc
} else {
  cat("No batch or study associations detected.\n")
}
if (!is.null(audit@target_assoc) && nrow(audit@target_assoc) > 0) {
  # Top features by target association score
  head(audit@target_assoc)
} else {
  cat("No target leakage scan results available.\n")
}
if (!is.null(audit@duplicates) && nrow(audit@duplicates) > 0) {
  # Top duplicate/near-duplicate pairs by similarity
  head(audit@duplicates)
} else {
  cat("No near-duplicates detected.\n")
}

## ----mechanism-summary--------------------------------------------------------
mech <- audit@info$mechanism_summary
if (is.data.frame(mech) && nrow(mech) > 0) {
  mech
} else {
  cat("No mechanism summary available.\n")
}

## ----plot-perm----------------------------------------------------------------
if (requireNamespace("ggplot2", quietly = TRUE)) {
  plot_perm_distribution(audit)
} else {
  cat("ggplot2 not installed; skipping permutation plot.\n")
}

## ----audit-by-learner---------------------------------------------------------
if (requireNamespace("ranger", quietly = TRUE) && 
    requireNamespace("parsnip", quietly = TRUE)) {
    
    # 1. Define specs
    # Standard GLM (no tuning)
    spec_glm <- parsnip::logistic_reg() |> 
        parsnip::set_engine("glm")
    
    # Random Forest
    spec_rf <- parsnip::rand_forest(
        mode = "classification",
        trees = 100
    ) |>
        parsnip::set_engine("ranger")
    
    # 2. Fit using the current split object
    fit_multi <- fit_resample(
        df,
        outcome = "outcome",
        splits = splits,
        learner = list(glm = spec_glm, rf = spec_rf),
        metrics = "auc"
    )
    
    # 3. Run the audit
    audits <- audit_leakage_by_learner(fit_multi, metric = "auc", B = 20)
    cat("Per-learner audit summary:\n")
    print(audits)
    
} else {
    cat("ranger/parsnip not installed.\n")
}

## ----audit-report, eval = FALSE-----------------------------------------------
# if (requireNamespace("rmarkdown", quietly = TRUE) && rmarkdown::pandoc_available()) {
#   report_path <- audit_report(audit, output_dir = ".")
#   cat("HTML report written to:\n", report_path, "\n")
# } else {
#   cat("rmarkdown or pandoc not available; skipping audit report rendering.\n")
# }

## ----time-series-fit----------------------------------------------------------
time_splits <- make_split_plan(
  df_time,
  outcome = "outcome",
  mode = "time_series",
  time = "time",
  v = 4,
  horizon = 1
)

cat("Time-series splits summary:\n")
time_splits

fit_time_leaky <- fit_resample(
  df_time,
  outcome = "outcome",
  splits = time_splits,
  learner = spec,
  metrics = "auc"
)

cat("Time-series leaky fit summary:\n")
summary(fit_time_leaky)

## ----time-series-safe---------------------------------------------------------
df_time_safe <- df_time
df_time_safe$leak_future <- NULL

fit_time_safe <- fit_resample(
  df_time_safe,
  outcome = "outcome",
  splits = time_splits,
  learner = spec,
  metrics = "auc"
)

cat("Time-series safe fit summary:\n")
summary(fit_time_safe)

audit_time <- audit_leakage(
  fit_time_safe,
  metric = "auc",
  B = 20,
  time_block = "stationary",
  block_len = 5
)

cat("Time-series leakage audit summary:\n")
summary(audit_time)
if (!is.null(audit_time@permutation_gap) && nrow(audit_time@permutation_gap) > 0) {
  # Time-series permutation significance results
  audit_time@permutation_gap
}

if (requireNamespace("ggplot2", quietly = TRUE)) {
  plot_time_acf(fit_time_safe, lag.max = 20)
} else {
  cat("ggplot2 not installed; skipping ACF plot.\n")
}

## ----cv-ci--------------------------------------------------------------------
# Standard 95% CI from the leakage-safe fit
ci_std <- cv_ci(fit_safe@metrics, level = 0.95, method = "normal")
cat("Standard 95% CI:\n")
print(ci_std)

# Nadeau-Bengio corrected CI
# For K-fold CV: n_train ≈ (K-1)/K × n, n_test ≈ n/K
K_folds  <- length(unique(fit_safe@metrics$fold))
n_total  <- nrow(fit_safe@splits@info$coldata)

ci_nb <- cv_ci(
  fit_safe@metrics,
  level   = 0.95,
  method  = "nadeau_bengio",
  n_train = round(n_total * (K_folds - 1L) / K_folds),
  n_test  = round(n_total / K_folds)
)
cat("Nadeau-Bengio corrected 95% CI:\n")
print(ci_nb)

## ----delta-lsi-data-----------------------------------------------------------
set.seed(100)
n_d    <- 120L
subj_d <- rep(paste0("P", seq_len(30L)), each = 4L)   # 30 subjects, 4 obs each

df_dlsi_guarded <- data.frame(
  subject = subj_d,
  outcome = factor(sample(c("case", "control"), n_d, replace = TRUE)),
  x1      = rnorm(n_d),
  x2      = rnorm(n_d)
)

# Leaky feature: subject-level mean of outcome
# With subject-grouped splits, test subjects' leak value encodes their outcome
df_dlsi_naive <- df_dlsi_guarded
df_dlsi_naive$leak <- ave(
  as.numeric(df_dlsi_guarded$outcome == "case"),
  subj_d, FUN = mean
)

# Shared splits: 5-fold × 5 repeats → R_eff = 5 (tier C: point + p-value)
splits_dlsi <- make_split_plan(
  df_dlsi_guarded,
  outcome = "outcome",
  mode    = "subject_grouped",
  group   = "subject",
  v       = 5L,
  repeats = 5L
)

## ----delta-lsi-fits-----------------------------------------------------------
# Reuse the parsnip spec (base R glm) defined earlier
# Naive fit: uses leaky feature (leak encodes test-subject outcome label)
fit_dlsi_naive <- fit_resample(
  df_dlsi_naive,
  outcome = "outcome",
  splits  = splits_dlsi,
  learner = spec,
  metrics = "auc",
  seed    = 1L
)

# Guarded fit: same splits, clean features only
fit_dlsi_guarded <- fit_resample(
  df_dlsi_guarded,
  outcome = "outcome",
  splits  = splits_dlsi,
  learner = spec,
  metrics = "auc",
  seed    = 1L
)

cat("Naive   AUC:", round(mean(fit_dlsi_naive@metrics$auc,   na.rm = TRUE), 3), "\n")
cat("Guarded AUC:", round(mean(fit_dlsi_guarded@metrics$auc, na.rm = TRUE), 3), "\n")

## ----delta-lsi-run------------------------------------------------------------
result_dlsi <- delta_lsi(
  fit_leaky   = fit_dlsi_naive,
  fit_guarded = fit_dlsi_guarded,
  metric      = "auc",
  M_boot      = 300L,
  M_flip      = 300L,
  seed        = 1L
)

summary(result_dlsi)

## ----delta-lsi-tier-----------------------------------------------------------
cat("Tier:          ", result_dlsi@tier, "\n")
cat("R_eff:         ", result_dlsi@R_eff, "\n")
cat("inference_ok:  ", result_dlsi@inference_ok, "\n")

## ----delta-lsi-components-----------------------------------------------------
# delta_metric: arithmetic mean of {Δ_r}
cat("delta_metric:  ", round(result_dlsi@delta_metric, 4), "\n")
# delta_lsi: Huber M-estimate of {Δ_r} (k=1.345, fixed MAD scale)
cat("delta_lsi:     ", round(result_dlsi@delta_lsi,    4), "\n")
# delta_lsi_ci: 95% BCa interval for delta_lsi (NA below tier B)
cat("95% BCa CI:    [",
    paste(round(result_dlsi@delta_lsi_ci, 4), collapse = ", "), "]\n")
# p_value: sign-flip randomization p-value (NA below tier C or if unpaired)
cat("p-value:       ",
    if (is.na(result_dlsi@p_value)) "NA (tier D or unpaired)"
    else format(result_dlsi@p_value, digits = 3), "\n")

## ----delta-lsi-unpaired, warning=TRUE-----------------------------------------
# Unpaired: naive uses sample-wise CV, guarded uses subject-grouped splits
splits_rand <- make_split_plan(
  df_dlsi_naive,
  outcome = "outcome",
  mode    = "subject_grouped",
  group   = "row_id",
  v       = 5L,
  repeats = 5L
)

fit_naive_rand <- fit_resample(
  df_dlsi_naive,
  outcome = "outcome",
  splits  = splits_rand,
  learner = spec,
  metrics = "auc",
  seed    = 1L
)

r_unpaired <- suppressWarnings(
  delta_lsi(fit_naive_rand, fit_dlsi_guarded, metric = "auc", seed = 1L)
)

cat(sprintf("Paired:        %s\n", r_unpaired@info$paired))
cat(sprintf("R_eff:         %d  (0 = unpaired, inference disabled)\n", r_unpaired@R_eff))
cat(sprintf("delta_metric:  %.4f  (naive raw mean - guarded raw mean)\n",
            r_unpaired@delta_metric))
cat(sprintf("p_value:       %s\n",
            if (is.na(r_unpaired@p_value)) "NA" else r_unpaired@p_value))

## ----delta-lsi-direction------------------------------------------------------
r_hib_true  <- suppressWarnings(delta_lsi(
  fit_dlsi_naive, fit_dlsi_guarded,
  metric           = "auc",
  higher_is_better = TRUE,
  seed             = 1L
))
r_hib_false <- suppressWarnings(delta_lsi(
  fit_dlsi_naive, fit_dlsi_guarded,
  metric           = "auc",
  higher_is_better = FALSE,
  seed             = 1L
))

cat(sprintf("hib = TRUE:   delta_metric = %+.4f  (positive = naive inflated)\n",
            r_hib_true@delta_metric))
cat(sprintf("hib = FALSE:  delta_metric = %+.4f  (sign flipped)\n",
            r_hib_false@delta_metric))

## ----delta-lsi-hib-slot-------------------------------------------------------
cat("higher_is_better used:", result_dlsi@info$higher_is_better, "\n")

## ----delta-lsi-slots----------------------------------------------------------
# Per-repeat summaries for custom plotting
head(result_dlsi@repeats_naive,   3)
head(result_dlsi@repeats_guarded, 3)

# Per-fold diagnostics
head(result_dlsi@folds_naive,   3)
head(result_dlsi@folds_guarded, 3)

# Metadata
str(result_dlsi@info)

## ----delta-lsi-audit-connect--------------------------------------------------
# Step 1: detect and characterise leakage
audit_naive_dlsi <- audit_leakage(fit_dlsi_naive, metric = "auc", B = 20L)
cat("Naive audit summary:\n")
summary(audit_naive_dlsi)

# Step 2: confirmed guards reduced inflation; check residual DLSI
cat("\ndelta_metric (naive vs guarded):", round(result_dlsi@delta_metric, 4), "\n")
if (!is.na(result_dlsi@p_value)) {
  cat("Sign-flip p-value:              ", format(result_dlsi@p_value, digits = 3), "\n")
}

## ----parallel-setup, eval=FALSE-----------------------------------------------
# library(future)
# 
# # Use multiple cores (works on all OS)
# plan(multisession, workers = 4)
# 
# # Run a heavy simulation
# sim <- simulate_leakage_suite(..., parallel = TRUE)
# 
# # Parallel folds or audits
# # fit_resample(..., parallel = TRUE)
# # audit_leakage(..., parallel = TRUE)
# # audit_leakage_by_learner(..., parallel_learners = TRUE)
# 
# # Return to sequential processing
# plan(sequential)

## ----simulate-suite-----------------------------------------------------------
if (requireNamespace("glmnet", quietly = TRUE)) {
  sim <- simulate_leakage_suite(
    n = 80,
    p = 6,
    mode = "subject_grouped",
    learner = "glmnet",
    leakage = "subject_overlap",
    seeds = 1:2,
    B = 20,
    parallel = FALSE,
    signal_strength = 1
  )
  
  # Simulation results (first 6 rows)
  head(sim)
} else {
  cat("glmnet not installed; skipping simulation suite example.\n")
}

