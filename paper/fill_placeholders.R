#!/usr/bin/env Rscript
## =================================================================
## Auto-substitute \placeholder{KEY} tokens in manuscript_jss.tex
## with computed values from placeholder_values.rds
## =================================================================
## Usage:  Rscript paper/fill_placeholders.R
## Effect: Reads manuscript_jss.tex, replaces every \placeholder{KEY}
##         with the corresponding value, writes the result in-place
##         (a .bak backup is created first).
## =================================================================

base_dir <- "paper"

## 1. Load placeholder values (named list: "[KEY]" -> value)
pv_file <- file.path(base_dir, "placeholder_values.rds")
if (!file.exists(pv_file)) {
  stop("placeholder_values.rds not found. Run compile_results.R first.")
}
placeholders <- readRDS(pv_file)
cat(sprintf("Loaded %d placeholder entries.\n", length(placeholders)))

## 2. Read manuscript
tex_file <- file.path(base_dir, "manuscript_jss.tex")
if (!file.exists(tex_file)) {
  stop("manuscript_jss.tex not found.")
}
tex <- readLines(tex_file, warn = FALSE)
tex_orig <- tex  # keep for comparison

## 3. Create backup
bak_file <- paste0(tex_file, ".bak")
writeLines(tex_orig, bak_file)
cat(sprintf("Backup written to %s\n", bak_file))

## 4. Replace each \placeholder{KEY} with its value
n_replaced <- 0L
for (ph_key in names(placeholders)) {
  ## ph_key is like "[KEY]" — extract the bare key

  bare_key <- gsub("^\\[|\\]$", "", ph_key)
  ## The LaTeX macro is \placeholder{KEY}
  pattern <- paste0("\\placeholder{", bare_key, "}")
  value   <- as.character(placeholders[[ph_key]])

  ## Use fixed-string replacement (no regex)
  hits <- grep(pattern, tex, fixed = TRUE)
  if (length(hits) > 0) {
    tex <- gsub(pattern, value, tex, fixed = TRUE)
    n_replaced <- n_replaced + length(hits)
    cat(sprintf("  %s -> %s  (%d lines)\n", pattern, value, length(hits)))
  }
}

## 5. Write result
writeLines(tex, tex_file)
cat(sprintf("\nDone: %d substitutions across %d placeholder keys.\n",
            n_replaced, length(placeholders)))

## 6. Check for any remaining \placeholder{} tokens
remaining <- grep("\\\\placeholder\\{", tex)
if (length(remaining) > 0) {
  cat(sprintf("WARNING: %d lines still contain \\placeholder{} tokens:\n",
              length(remaining)))
  for (ln in remaining) {
    cat(sprintf("  Line %d: %s\n", ln, trimws(tex[ln])))
  }
} else {
  cat("All placeholder tokens have been replaced.\n")
}
