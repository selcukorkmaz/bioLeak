## ---------------------------------------------------------------------
## test-vignette-suggests.R
## ---------------------------------------------------------------------
## Tests for the 0.3.7 cleanup of Suggests-package handling in the
## vignette. The editor's Comment 9 noted that the vignette mixed
## `requireNamespace("recipes", quietly = TRUE)` (defensive) with
## `library(recipes)` (will hard-error if recipes is not installed),
## defeating the defensive check. After the fix, every code chunk
## that calls `library(<suggests_pkg>)` in eval = TRUE mode must be
## gated by a chunk option `eval = requireNamespace(...)` or by an
## inline `if (requireNamespace(...))` guard.
## ---------------------------------------------------------------------

vignette_lines <- function() {
  paths <- c(
    file.path("..", "..", "vignettes", "bioLeak-intro.Rmd"),
    file.path("vignettes", "bioLeak-intro.Rmd"),
    system.file("doc", "bioLeak-intro.Rmd", package = "bioLeak"),
    system.file("doc", "bioLeak-intro.R", package = "bioLeak")
  )
  paths <- paths[nzchar(paths) & file.exists(paths)]
  if (!length(paths)) skip("vignette source not found in this build")
  readLines(paths[1])
}

## Packages listed in DESCRIPTION's Suggests field. These are the
## packages whose `library()` calls must be guarded.
suggests_pkgs <- function() {
  desc_path <- c(
    file.path("..", "..", "DESCRIPTION"),
    "DESCRIPTION",
    system.file("DESCRIPTION", package = "bioLeak")
  )
  desc_path <- desc_path[nzchar(desc_path) & file.exists(desc_path)]
  if (!length(desc_path)) skip("DESCRIPTION not found in this build")
  d <- read.dcf(desc_path[1])
  if (!"Suggests" %in% colnames(d)) return(character(0))
  raw <- d[1, "Suggests"]
  pkgs <- trimws(strsplit(raw, ",")[[1]])
  ## Strip version constraints like " (>= 3.0.0)"
  pkgs <- sub(" *\\(.*\\)$", "", pkgs)
  pkgs[nzchar(pkgs)]
}

test_that("Vignette has no unguarded `library(<suggests_pkg>)` calls", {
  lines <- vignette_lines()
  pkgs  <- suggests_pkgs()
  if (!length(pkgs)) skip("no Suggests packages declared")

  ## Walk the vignette, tracking which chunk each line is in and
  ## whether that chunk is gated by an `eval = ...` option that
  ## references requireNamespace() or by an explicit eval = FALSE.
  in_chunk <- FALSE
  chunk_header <- ""
  offending <- character(0)

  for (i in seq_along(lines)) {
    L <- lines[i]
    if (grepl("^```\\{r", L)) {
      in_chunk <- TRUE
      chunk_header <- L
      next
    }
    if (grepl("^```\\s*$", L)) {
      in_chunk <- FALSE
      chunk_header <- ""
      next
    }
    if (!in_chunk) next

    ## Match `library(<pkg>)` (top-level call, not commented).
    m <- regmatches(L, regexec("^\\s*library\\(([A-Za-z0-9._]+)\\)", L))[[1]]
    if (length(m) == 2L) {
      pkg <- m[2]
      if (!pkg %in% pkgs) next   # Imports/Depends → fine
      ## Suggests pkg → check if the chunk header contains a guard
      header_safe <-
        grepl("eval\\s*=\\s*FALSE", chunk_header) ||
        grepl(sprintf("requireNamespace\\(\\s*[\"']%s[\"']", pkg),
              chunk_header)
      if (!header_safe) {
        offending <- c(offending,
                       sprintf("line %d: %s   (chunk: %s)", i, L, chunk_header))
      }
    }
  }

  expect_length(offending, 0L)
})

test_that("Suggests-using chunks at the documented locations are properly gated", {
  lines <- vignette_lines()

  ## tidymodels-interop chunk header should carry an
  ## `eval = requireNamespace(...)` gate.
  idx <- grep("^```\\{r tidymodels-interop", lines)
  expect_length(idx, 1L)
  hdr <- lines[idx]
  expect_match(hdr, "requireNamespace\\(\\s*['\"]recipes['\"]")
  expect_match(hdr, "requireNamespace\\(\\s*['\"]yardstick['\"]")

  ## parallel-setup chunk should remain eval = FALSE (documentation-only).
  idx <- grep("^```\\{r parallel-setup", lines)
  expect_length(idx, 1L)
  expect_match(lines[idx], "eval\\s*=\\s*FALSE")
})

test_that("Defensive requireNamespace guards remain in place for all gated chunks", {
  lines <- vignette_lines()
  ## Confirm no chunks regressed to bare library() inside an unguarded
  ## chunk (the editor's specific complaint pattern).
  ## Count: at least one inline `requireNamespace(<pkg>)` guard for every
  ## Suggests package that the vignette uses.
  used_with_guard <- regmatches(
    lines,
    gregexpr("requireNamespace\\(\\s*['\"]([^'\"]+)['\"]", lines)
  )
  flat <- unique(unlist(used_with_guard))
  pkg_names <- gsub("requireNamespace\\(\\s*['\"]([^'\"]+)['\"]", "\\1", flat)
  pkg_names <- pkg_names[nzchar(pkg_names)]

  ## At minimum, recipes/yardstick guards must exist somewhere in the
  ## vignette (either in chunk options or inline `if`).
  expect_true("recipes"   %in% pkg_names)
  expect_true("yardstick" %in% pkg_names)
})
