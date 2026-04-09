## =================================================================
## Shared helper functions for bioLeak paper scripts
## =================================================================

#' Load and prepare the curatedOvarianData case study dataset
#'
#' Loads all eligible ExpressionSets from curatedOvarianData, creates a
#' 3-year overall survival binary endpoint, finds common genes across
#' studies, and returns a combined data.frame ready for modelling.
#'
#' @param max_genes Integer. Maximum number of genes to retain (most
#'   variable by variance in the first study). Default 2000.
#' @param min_common Integer. Minimum common genes before switching to
#'   pairwise study selection. Default 100.
#' @param min_overlap Integer. Minimum genes a study must share with the
#'   selected subset to be included. Default 500.
#' @param verbose Logical. Print progress messages. Default TRUE.
#' @return A list with components:
#'   \item{df_combined}{data.frame with gene columns, y (factor), and study}
#'   \item{X_combined}{numeric matrix of gene expression values}
#'   \item{combined_y}{integer vector of binary outcomes}
#'   \item{combined_study}{character vector of study labels}
#'   \item{study_list}{list of per-study metadata}
#'   \item{common_genes}{character vector of gene names used}
load_ovarian_casestudy <- function(max_genes = 2000, min_common = 100,
                                   min_overlap = 500, verbose = TRUE) {
  if (!requireNamespace("curatedOvarianData", quietly = TRUE))
    stop("curatedOvarianData package is required.")

  if (verbose) cat("Loading curatedOvarianData...\n")
  all_ds <- data(package = "curatedOvarianData")$results[, "Item"]
  all_ds <- gsub("\\s.*", "", all_ds)

  ## Load each dataset and check eligibility
  study_list <- list()
  for (ds_name in all_ds) {
    tryCatch({
      data(list = ds_name, package = "curatedOvarianData", envir = environment())
      eset <- get(ds_name, envir = environment())
      if (!inherits(eset, "ExpressionSet")) next

      pd <- Biobase::pData(eset)
      has_os <- "days_to_death" %in% names(pd) && "vital_status" %in% names(pd)
      if (!has_os) next

      os_days <- as.numeric(pd$days_to_death)
      vital   <- pd$vital_status

      valid <- !is.na(os_days) & !is.na(vital)
      if (sum(valid) < 50) next

      os_binary <- ifelse(os_days >= 1095, 1,
                          ifelse(vital == "deceased" & os_days < 1095, 0, NA))
      valid2 <- valid & !is.na(os_binary)
      if (sum(valid2) < 50) next
      if (length(unique(os_binary[valid2])) < 2) next

      study_list[[ds_name]] <- list(
        eset = eset, valid_idx = which(valid2),
        os_binary = os_binary[valid2], n = sum(valid2)
      )
      if (verbose)
        cat(sprintf("  %s: n=%d (class 0: %d, class 1: %d)\n",
                    ds_name, sum(valid2),
                    sum(os_binary[valid2] == 0), sum(os_binary[valid2] == 1)))
    }, error = function(e) NULL)
  }

  if (verbose) cat(sprintf("\nEligible studies: %d\n", length(study_list)))

  ## Find common genes
  gene_lists <- lapply(study_list, function(x) Biobase::featureNames(x$eset))
  common_genes <- Reduce(intersect, gene_lists)
  if (verbose) cat(sprintf("Common genes across all studies: %d\n", length(common_genes)))

  ## If too few common genes, select studies with maximum overlap
  if (length(common_genes) < min_common) {
    if (verbose) cat("Too few common genes, selecting studies with maximum overlap...\n")
    study_names <- names(study_list)
    best_pair <- NULL; best_common <- 0
    for (i in seq_along(study_names)) {
      for (j in seq_along(study_names)) {
        if (j <= i) next
        cg <- length(intersect(gene_lists[[i]], gene_lists[[j]]))
        if (cg > best_common) { best_common <- cg; best_pair <- c(i, j) }
      }
    }
    selected <- study_names[best_pair]
    common_genes <- intersect(gene_lists[[best_pair[1]]], gene_lists[[best_pair[2]]])

    for (sn in setdiff(study_names, selected)) {
      cg <- intersect(common_genes, gene_lists[[sn]])
      if (length(cg) >= min_overlap) { common_genes <- cg; selected <- c(selected, sn) }
    }
    study_list <- study_list[selected]
    if (verbose) cat(sprintf("Selected %d studies with %d common genes\n",
                             length(study_list), length(common_genes)))
  }

  ## Limit to top most-variable genes
  if (length(common_genes) > max_genes) {
    first_eset <- study_list[[1]]$eset
    gene_var <- apply(Biobase::exprs(first_eset)[common_genes, ], 1, var, na.rm = TRUE)
    common_genes <- names(sort(gene_var, decreasing = TRUE))[seq_len(max_genes)]
  }

  ## Build combined matrix
  combined_X <- list(); combined_y <- c(); combined_study <- c(); sample_ids <- c()
  for (sn in names(study_list)) {
    sl <- study_list[[sn]]
    expr_mat <- t(Biobase::exprs(sl$eset)[common_genes, sl$valid_idx, drop = FALSE])
    combined_X[[sn]] <- expr_mat
    combined_y <- c(combined_y, sl$os_binary)
    combined_study <- c(combined_study, rep(sn, sl$n))
    sample_ids <- c(sample_ids, paste0(sn, "_", seq_len(sl$n)))
  }

  X_combined <- do.call(rbind, combined_X)
  rownames(X_combined) <- sample_ids

  if (verbose) {
    cat(sprintf("\nCombined dataset: N=%d samples, S=%d studies, G=%d genes\n",
                nrow(X_combined), length(unique(combined_study)), ncol(X_combined)))
    cat(sprintf("Class distribution: 0=%d (%.1f%%), 1=%d (%.1f%%)\n",
                sum(combined_y == 0), 100 * mean(combined_y == 0),
                sum(combined_y == 1), 100 * mean(combined_y == 1)))
  }

  df_combined <- data.frame(X_combined, check.names = TRUE)
  df_combined$y <- factor(combined_y)
  df_combined$study <- combined_study

  list(
    df_combined = df_combined,
    X_combined = X_combined,
    combined_y = combined_y,
    combined_study = combined_study,
    study_list = study_list,
    common_genes = common_genes
  )
}
