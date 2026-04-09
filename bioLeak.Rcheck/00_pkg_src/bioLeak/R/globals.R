# Suppress NSE notes for ggplot2 aesthetics in R CMD check.
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    "acf",
    "count",
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
    "series",
    "truth",
    "type",
    "level",
    "value"
  ))
}
