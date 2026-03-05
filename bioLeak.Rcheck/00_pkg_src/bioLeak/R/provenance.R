# Provenance capture for reproducibility tracking ---------------------

#' Capture session provenance information
#'
#' Internal helper that records R version, loaded packages, platform,
#' timestamp, and optionally git SHA and hardware details.
#'
#' @param capture_git Logical; if TRUE, attempts to capture the current
#'   git short SHA via \code{system2("git", ...)}. Default FALSE.
#' @param capture_hardware Logical; if TRUE, captures basic hardware info
#'   from \code{Sys.info()}. Default TRUE.
#' @return A list with elements \code{r_version}, \code{packages} (data.frame),
#'   \code{platform}, \code{timestamp}, \code{git_sha}, and \code{hardware}.
#' @keywords internal
.bio_capture_provenance <- function(capture_git = TRUE, capture_hardware = TRUE) {
  si <- utils::sessionInfo()

  # Extract package versions from sessionInfo
  base_pkgs <- data.frame(
    package = si$basePkgs,
    version = vapply(si$basePkgs, function(p) {
      tryCatch(as.character(utils::packageVersion(p)), error = function(e) NA_character_)
    }, character(1)),
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  other_pkgs <- if (length(si$otherPkgs)) {
    data.frame(
      package = vapply(si$otherPkgs, function(p) p$Package, character(1)),
      version = vapply(si$otherPkgs, function(p) p$Version, character(1)),
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  } else {
    data.frame(package = character(0), version = character(0),
               stringsAsFactors = FALSE)
  }

  packages <- rbind(base_pkgs, other_pkgs)

  # Git SHA
 git_sha <- NA_character_
  if (isTRUE(capture_git)) {
    git_sha <- tryCatch({
      res <- system2("git", c("rev-parse", "--short", "HEAD"),
                      stdout = TRUE, stderr = TRUE)
      if (length(res) && !inherits(res, "try-error") && nchar(res[1]) > 0) {
        res[1]
      } else {
        NA_character_
      }
    }, error = function(e) NA_character_,
       warning = function(w) NA_character_)
  }

  # Hardware info
  hardware <- NULL
  if (isTRUE(capture_hardware)) {
    sinfo <- Sys.info()
    hardware <- list(
      sysname = unname(sinfo["sysname"]),
      nodename = unname(sinfo["nodename"]),
      machine = unname(sinfo["machine"])
    )
  }

  list(
    r_version = paste0(R.version$major, ".", R.version$minor),
    packages = packages,
    platform = si$platform,
    timestamp = Sys.time(),
    git_sha = git_sha,
    hardware = hardware
  )
}
