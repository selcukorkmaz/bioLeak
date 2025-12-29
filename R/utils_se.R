# Utilities to work with SummarizedExperiment (SE) or plain matrices

.bio_is_se <- function(x) inherits(x, "SummarizedExperiment")

.bio_get_x <- function(x, assay_name = NULL) {
  if (.bio_is_se(x)) {
    if (is.null(assay_name)) {
      an <- SummarizedExperiment::assayNames(x)
      if (length(an) == 0) stop("No assays found in SummarizedExperiment.")
      assay_name <- an[1L]
    }
    out <- t(SummarizedExperiment::assay(x, assay_name))  # samples x features
    out <- as.data.frame(out, stringsAsFactors = FALSE)
  } else if (is.data.frame(x)) {
    out <- x  # keep original structure
  } else if (is.matrix(x)) {
    # convert to numeric data.frame if possible
    out <- as.data.frame(x, stringsAsFactors = FALSE)
    # attempt numeric conversion for columns that look numeric
    out[] <- lapply(out, function(col) {
      suppressWarnings(num <- as.numeric(col))
      if (is.numeric(col)) {
        col
      } else if (!anyNA(num)) {
        num
      } else {
        col
      }
    })
  } else {
    stop("Unsupported x type.")
  }
  return(out)
}


.bio_get_y <- function(x, outcome) {
  if (.bio_is_se(x)) {
    cd <- SummarizedExperiment::colData(x)
    if (!(outcome %in% colnames(cd))) stop("Outcome column not in colData.")
    cd[[outcome]]
  } else {
    if (is.null(outcome))
      stop("Provide outcome column name when x is not a SummarizedExperiment.")
    if (is.character(outcome) && outcome %in% colnames(x)) {
      if (is.data.frame(x)) {
        x[[outcome]]
      } else if (is.matrix(x)) {
        x[, outcome]
      } else {
        x[[outcome]]
      }
    } else {
      stop("Outcome not found in data frame or matrix.")
    }
  }
}

.bio_get_meta <- function(x, cols) {
  if (.bio_is_se(x)) {
    cd <- SummarizedExperiment::colData(x)
    out <- list()
    for (nm in cols) {
      out[[nm]] <- if (nm %in% colnames(cd)) cd[[nm]] else NULL
    }
    out
  } else {
    as.list(rep(list(NULL), length(cols))) |> stats::setNames(cols)
  }
}

.bio_hash_indices <- function(idx) {
  if (requireNamespace("digest", quietly = TRUE))
    return(paste0("h", substr(digest::digest(idx), 1, 8)))
  paste0("h", sprintf("%08X", as.integer(sum(unlist(idx)) %% .Machine$integer.max)))
}

.bio_is_classification <- function(y) is.factor(y) || (is.numeric(y) && length(unique(y)) <= 10)
.bio_is_regression <- function(y) is.numeric(y) && !.bio_is_binomial(y)
.bio_is_binomial <- function(y) {
  if (is.factor(y)) return(nlevels(y) == 2)
  if (is.numeric(y)) {
    u <- sort(unique(y))
    return(length(u) == 2 && all(u %in% c(0,1)))
  }
  FALSE
}
