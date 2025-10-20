#' Render an HTML leakage audit report
#' @param audit LeakAudit
#' @param file output file path ending with .html
#' @param theme R Markdown theme ("flatly", "cosmo", "default")
#' @param self_contained logical; embed dependencies (default TRUE)
#' @param logo optional logo path or URL
#' @return invisible file path
#' @export
render_report <- function(audit,
                          file = "bioLeak_report.html",
                          theme = "flatly",
                          self_contained = TRUE,
                          logo = NULL) {
  if (!inherits(audit, "LeakAudit"))
    stop("'audit' must be a 'LeakAudit' object.", call. = FALSE)
  if (!requireNamespace("rmarkdown", quietly = TRUE))
    stop("Package 'rmarkdown' is required.", call. = FALSE)
  if (!rmarkdown::pandoc_available())
    stop("Pandoc is not available. Install RStudio or Pandoc manually.", call. = FALSE)

  # Locate or build template
  template <- system.file("rmd", "bioLeak_report.Rmd", package = "bioLeak")
  if (template == "" || is.null(template)) {
    template <- tempfile(fileext = ".Rmd")
    writeLines(.bio_report_template_inline(theme = theme, logo = logo), con = template, useBytes = TRUE)
  }

  # Output directory
  outdir <- normalizePath(dirname(file), mustWork = FALSE)
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  outfile <- normalizePath(file, mustWork = FALSE)

  opts <- list(
    input = template,
    output_file = outfile,
    params = list(audit = audit),
    quiet = TRUE,
    envir = new.env(parent = globalenv())
  )

  res <- tryCatch({
    rmarkdown::render(
      input = opts$input,
      output_file = opts$output_file,
      params = opts$params,
      quiet = opts$quiet,
      envir = opts$envir,
      output_options = list(
        self_contained = self_contained,
        theme = theme,
        toc = TRUE,
        toc_depth = 2
      )
    )
  }, error = function(e) {
    warning("HTML rendering failed. Saving text summary instead: ", conditionMessage(e))
    txtfile <- sub("\\.html$", ".txt", outfile)
    capture.output(str(audit), file = txtfile)
    return(txtfile)
  })

  invisible(res %||% outfile)
}

.bio_report_template_inline <- function(theme = "flatly", logo = NULL) {
  paste0(
    "---\n",
    "title: \"bioLeak Leakage Audit Report\"\n",
    "output:\n",
    "  html_document:\n",
    "    theme: ", theme, "\n",
    "    toc: true\n",
    "    toc_depth: 2\n",
    "params:\n",
    "  audit: NULL\n",
    "---\n\n",
    if (!is.null(logo)) paste0("![](", logo, "){width=120px}\n\n") else "",
    "```{r setup, include=FALSE}\n",
    "knitr::opts_chunk$set(echo=FALSE, message=FALSE, warning=FALSE)\n",
    "aud <- params$audit\n",
    "fit <- aud@fit\n",
    "```\n\n",
    "## Overview\n",
    "- Split mode: `r fit@splits@mode`  \n",
    "- Outcome: `r fit@outcome`  \n",
    "- Task: `r fit@task`  \n",
    "- Hash: `r fit@info$hash`  \n\n",
    "## Cross-validated Metrics\n",
    "```{r}\n",
    "if (nrow(fit@metrics) > 0) knitr::kable(fit@metrics) else cat('No metrics available.')\n",
    "```\n\n",
    "## Permutation Test\n",
    "```{r, fig.height=4}\n",
    "if (!is.null(aud@perm_distribution)) {\n",
    "  if (requireNamespace('ggplot2', quietly = TRUE)) {\n",
    "    df <- data.frame(metric = aud@perm_distribution)\n",
    "    ggplot2::ggplot(df, ggplot2::aes(metric)) +\n",
    "      ggplot2::geom_histogram(bins=30, fill='grey70') +\n",
    "      ggplot2::geom_vline(xintercept = aud@permutation_gap$metric_obs, color='red', linetype='dashed') +\n",
    "      ggplot2::labs(title='Permutation Test Distribution')\n",
    "  } else print(aud@permutation_gap)\n",
    "} else cat('No permutation results.')\n",
    "```\n\n",
    "## Batch / Study Association\n",
    "```{r}\n",
    "if (!is.null(aud@batch_assoc) && nrow(aud@batch_assoc) > 0) knitr::kable(aud@batch_assoc) else cat('No batch association results.')\n",
    "```\n\n",
    "## Near-Duplicates\n",
    "```{r}\n",
    "if (!is.null(aud@duplicates) && nrow(aud@duplicates) > 0) utils::head(aud@duplicates, 15) else cat('No duplicates found above threshold.')\n",
    "```\n\n",
    "## Reproducibility\n",
    "Generated on: `r Sys.time()`  \n",
    "R version: `r getRversion()`  \n",
    "bioLeak version: `r utils::packageVersion('bioLeak')`\n"
  )
}
