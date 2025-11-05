#!/usr/bin/env Rscript

# Install R and Bioconductor dependencies required by DearMeta.

cran_packages <- c(
  optparse = "1.7.5",
  jsonlite = "2.0.0",
  "data.table" = "1.17.8",
  ggplot2 = "4.0.0",
  plotly = "4.11.0",
  htmlwidgets = "1.6.4",
  DT = "0.34.0",
  RColorBrewer = "1.1-3",
  VennDiagram = "1.7.3"
)

bioc_packages <- c(
  minfi = "1.56.0",
  sesame = "1.28.0",
  sesameData = "1.27.1",
  limma = "3.66.0",
  sva = "3.58.0",
  DMRcate = "3.6.0",
  IlluminaHumanMethylationEPICanno.ilm10b4.hg19 = "0.6.0",
  IlluminaHumanMethylationEPICv2anno.20a1.hg38 = "1.0.0",
  IlluminaHumanMethylation450kanno.ilmn12.hg19 = "0.6.1",
  IlluminaHumanMethylation450kmanifest = "0.4.0",
  FlowSorted.Blood.EPIC = "2.13.0",
  FlowSorted.Blood.450k = "1.47.0"
)

ensure_cran <- function(packages) {
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes", repos = "https://cloud.r-project.org")
  }
  for (pkg in names(packages)) {
    target_version <- packages[[pkg]]
    installed <- if (requireNamespace(pkg, quietly = TRUE)) as.character(packageVersion(pkg)) else NULL
    if (!identical(installed, target_version)) {
      message(sprintf("Installing %s@%s from CRAN", pkg, target_version))
      remotes::install_version(pkg, version = target_version, repos = "https://cloud.r-project.org", upgrade = "never", quiet = TRUE)
    }
  }
}

ensure_bioc <- function(packages, bioc_version = "3.22") {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }
  for (pkg in names(packages)) {
    target_version <- packages[[pkg]]
    installed <- if (requireNamespace(pkg, quietly = TRUE)) as.character(packageVersion(pkg)) else NULL
    if (!identical(installed, target_version)) {
      message(sprintf("Installing %s@%s from Bioconductor %s", pkg, target_version, bioc_version))
      BiocManager::install(pkg, ask = FALSE, update = FALSE, version = bioc_version)
    }
    installed_now <- if (requireNamespace(pkg, quietly = TRUE)) as.character(packageVersion(pkg)) else NULL
    if (!identical(installed_now, target_version)) {
      stop(sprintf("Pinned version %s@%s could not be installed", pkg, target_version))
    }
  }
}

ensure_cran(cran_packages)
ensure_bioc(bioc_packages)

prefetch_sesame_resources <- function() {
  if (!requireNamespace("sesameData", quietly = TRUE)) {
    message("sesameData package is unavailable; skipping resource prefetch.")
    return(invisible(FALSE))
  }
  required_titles <- c("KYCG.EPIC.Mask.20220123")
  for (title in required_titles) {
    message("Prefetching sesameData resource: ", title)
    tryCatch(
      {
        invisible(sesameData::sesameDataGet(title))
        message("  cached ", title)
      },
      error = function(e) {
        message("  failed to cache ", title, ": ", conditionMessage(e))
      }
    )
  }
  invisible(TRUE)
}

prefetch_sesame_resources()

message("All R dependencies installed.")
