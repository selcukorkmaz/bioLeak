# Suppress NSE notes for ggplot2 aesthetics in R CMD check.
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    "acf",
    "count",
    "delta",
    "fold",
    "lag",
    "metric",
    "mid",
    "n",
    "obs_rate",
    "pred_mean",
    ".pred",
    ".pred_class",
    "prop_scaled",
    "repeat_idx",
    "series",
    "truth",
    "type",
    "level",
    "value"
  ))
}
