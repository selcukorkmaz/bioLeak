# Suppress NSE notes for ggplot2 aesthetics in R CMD check.
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    "acf",
    "count",
    "fold",
    "lag",
    "metric",
    "mid",
    "prop_scaled",
    "series",
    "type",
    "value"
  ))
}
