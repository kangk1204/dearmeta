#!/usr/bin/env Rscript

cran_repo <- "https://cloud.r-project.org"

ensure_cran <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    message("Installing CRAN packages: ", paste(missing, collapse = ", "))
    install.packages(missing, repos = cran_repo)
  }
}

ensure_cran(c("BiocManager", "png"))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  stop("BiocManager could not be loaded even after installation")
}

message("Setting Bioconductor version 3.22")
BiocManager::install(version = "3.22", ask = FALSE, update = FALSE)

preload <- c(
  "XVector",
  "Biostrings",
  "Rhtslib",
  "Rsamtools",
  "KEGGREST",
  "Rhdf5lib",
  "rhdf5",
  "rhdf5filters",
  "SparseArray",
  "DelayedArray",
  "DelayedMatrixStats",
  "HDF5Array",
  "SummarizedExperiment",
  "AnnotationDbi",
  "GenomicAlignments",
  "GEOquery",
  "GenomicFeatures",
  "AnnotationHub",
  "ExperimentHub",
  "BSgenome",
  "ShortRead",
  "cigarillo",
  "XML",
  "restfulr"
)

message("Installing targeted Bioconductor prerequisites")
BiocManager::install(preload, dependencies = TRUE, ask = FALSE, update = FALSE)

message("Bioconductor priming complete.")
