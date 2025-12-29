# One-time installs (if needed):
# install.packages("BiocManager")
# BiocManager::install("GEOquery")
# install.packages("remotes")
# remotes::install_github("selcukorkmaz/bioLeak")

library(GEOquery)
library(Biobase)
library(bioLeak)

gse <- getGEO("GSE2034", GSEMatrix = TRUE)
eset <- gse[[1]]

meta <- pData(eset)

# Bone relapse status is encoded in characteristics_ch1 for this series.
outcome_num <- as.integer(sub(".*: ", "", meta$characteristics_ch1))
outcome <- factor(outcome_num, levels = c(0, 1), labels = c("no_relapse", "relapse"))

X <- t(exprs(eset))

df <- data.frame(
  outcome = outcome,
  subject_id = meta$geo_accession,
  X,
  check.names = FALSE
)

splits <- make_splits(
  df,
  outcome = "outcome",
  mode = "subject_grouped",
  group = "subject_id",
  stratify = TRUE,
  v = 5
)

splits

fit <- fit_resample(
  df,
  outcome = "outcome",
  splits = splits,
  preprocess = list(
    impute = list(method = "median"),
    normalize = list(method = "zscore"),
    filter = list(var_thresh = 0.01, iqr_thresh = 0),
    fs = list(method = "ttest")
  ),
  learner = "glmnet",
  metrics = c("auc")
)

summary(fit)

audit <- audit_leakage(
  fit,
  metric = "auc",
  B = 100,
  X_ref = X,
  sim_method = "pearson",
  sim_threshold = 0.99
)

summary(audit)

audit_report(audit, output_dir = ".")
