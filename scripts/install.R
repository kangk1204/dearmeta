#!/usr/bin/env Rscript

# Install R and Bioconductor dependencies required by DearMeta.

cran_packages <- c(
  "optparse",
  "jsonlite",
  "data.table",
  "ggplot2",
  "plotly",
  "htmlwidgets",
  "DT",
  "RColorBrewer",
  "VennDiagram"
)

bioc_packages <- c(
  "minfi",
  "sesame",
  "sesameData",
  "limma",
  "sva",
  "DMRcate",
  "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
  "IlluminaHumanMethylationEPICv2anno.20a1.hg38",
  "IlluminaHumanMethylation450kanno.ilmn12.hg19",
  "IlluminaHumanMethylation450kmanifest"
)

ensure_cran <- function(packages) {
  to_install <- setdiff(packages, rownames(installed.packages()))
  if (length(to_install) > 0) {
    install.packages(to_install, repos = "https://cloud.r-project.org")
  }
}

ensure_bioc <- function(packages) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }
  to_install <- setdiff(packages, rownames(installed.packages()))
  if (length(to_install) > 0) {
    BiocManager::install(to_install, ask = FALSE, update = TRUE)
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
