#!/usr/bin/env Rscript
## =================================================================
## Build rendered manuscript by substituting placeholder macros
## =================================================================
## Reads placeholder_values.rds (named list: "[NAME]" -> "value")
## and replaces every \placeholder{NAME} in manuscript_jss.tex
## with the corresponding value. Writes manuscript_jss_rendered.tex.
## =================================================================

cat("=== Building Manuscript ===\n\n")

base_dir <- "paper"

## 1. Load placeholder values
ph_file <- file.path(base_dir, "placeholder_values.rds")
if (!file.exists(ph_file)) {
  stop("placeholder_values.rds not found. Run compile_results.R first.")
}
placeholders <- readRDS(ph_file)
cat(sprintf("Loaded %d placeholder values\n", length(placeholders)))

## 2. Read manuscript source
tex_file <- file.path(base_dir, "manuscript_jss.tex")
if (!file.exists(tex_file)) {
  stop("manuscript_jss.tex not found.")
}
tex <- readLines(tex_file, warn = FALSE)
tex_text <- paste(tex, collapse = "\n")

## 3. Replace \placeholder{NAME} with the corresponding value
n_replaced <- 0L
for (key in names(placeholders)) {
  ## key format: "[NAME]" — extract the inner name
  name <- gsub("^\\[|\\]$", "", key)
  ## Escape special regex characters in the name
  name_escaped <- gsub("([\\\\\\[\\]{}()+*?^$.|])", "\\\\\\1", name)
  pattern <- paste0("\\\\placeholder\\{", name_escaped, "\\}")
  value <- placeholders[[key]]

  count <- length(gregexpr(pattern, tex_text)[[1]])
  if (count > 0 && gregexpr(pattern, tex_text)[[1]][1] != -1) {
    tex_text <- gsub(pattern, value, tex_text)
    n_replaced <- n_replaced + count
    cat(sprintf("  %s -> %s (%d occurrences)\n", key, value, count))
  }
}

## 4. Check for unresolved placeholders
unresolved <- gregexpr("\\\\placeholder\\{[^}]+\\}", tex_text)[[1]]
n_unresolved <- if (unresolved[1] == -1) 0L else length(unresolved)

## 5. Write rendered output
out_file <- file.path(base_dir, "manuscript_jss_rendered.tex")
writeLines(strsplit(tex_text, "\n")[[1]], out_file)

cat(sprintf("\nTotal replacements: %d\n", n_replaced))
if (n_unresolved > 0) {
  ## Extract unresolved placeholder names for reporting
  matches <- regmatches(tex_text, gregexpr("\\\\placeholder\\{[^}]+\\}", tex_text))[[1]]
  cat(sprintf("WARNING: %d unresolved placeholder(s) remain:\n", n_unresolved))
  for (m in unique(matches)) cat(sprintf("  %s\n", m))
} else {
  cat("All placeholders resolved.\n")
}

cat(sprintf("Rendered manuscript written to %s\n", out_file))
cat("=== Done ===\n")
