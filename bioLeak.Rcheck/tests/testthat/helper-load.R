# During R CMD check, bioLeak is loaded via tests/testthat.R (library(bioLeak)).
# During interactive development, load via pkgload if not already attached.
if (!"package:bioLeak" %in% search()) {
  pkg_root <- normalizePath("../..", winslash = "/", mustWork = FALSE)
  if (file.exists(file.path(pkg_root, "DESCRIPTION"))) {
    pkgload::load_all(pkg_root, export_all = FALSE)
  }
}
