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

message("All R dependencies installed.")
