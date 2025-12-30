if (!"package:bioLeak" %in% search()) {
  pkg_root <- normalizePath("../..", winslash = "/", mustWork = FALSE)
  pkgload::load_all(pkg_root, export_all = FALSE)
}
