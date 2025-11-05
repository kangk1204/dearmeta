#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
target_version <- if (length(args) >= 1 && nzchar(args[[1]])) args[[1]] else NA_character_

cran_install <- function(pkg) {
  install.packages(pkg, repos = "https://cloud.r-project.org")
}

ensure_remotes <- function() {
  if (!requireNamespace("remotes", quietly = TRUE)) {
    cran_install("remotes")
  }
}

needs_install <- !requireNamespace("renv", quietly = TRUE)
if (!is.na(target_version) && !needs_install) {
  installed <- as.character(packageVersion("renv"))
  needs_install <- !identical(installed, target_version)
}

if (needs_install) {
  if (!is.na(target_version)) {
    ensure_remotes()
    remotes::install_version(
      "renv",
      version = target_version,
      repos = "https://cloud.r-project.org",
      upgrade = "never"
    )
  } else {
    cran_install("renv")
  }
  message("Installed renv", if (!is.na(target_version)) paste0("@", target_version))
} else {
  message("renv@", packageVersion("renv"), " already present; skipping install")
}
