#!/usr/bin/env Rscript

# Install R and Bioconductor dependencies required by DearMeta.

cran_packages <- c(
  "optparse",
  "jsonlite",
  "data.table",
  "remotes",
  "msigdbr",
  "fgsea",
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
  "fgsea",
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

install_dmrff <- function(ref = "perishky/dmrff@8e0469e5238c4c2de84746af143a733600537be4") {
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes", repos = "https://cloud.r-project.org")
  }
  needs_install <- !requireNamespace("dmrff", quietly = TRUE)
  if (!needs_install) {
    sha <- tryCatch(utils::packageDescription("dmrff", fields = c("RemoteSha", "GithubSHA1")), error = function(...) NA_character_)
    sha_str <- sha[!is.na(sha) & nzchar(sha)]
    sha_str <- if (length(sha_str) > 0) sha_str[1] else NA_character_
    needs_install <- is.na(sha_str) || !grepl("8e0469e", sha_str, fixed = TRUE)
  }
  if (needs_install) {
    remotes::install_github(ref, upgrade = "never")
  }
}

install_dmrff()

download_msigdb_fgsea <- function(cache_dir) {
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  cache_file <- file.path(cache_dir, "msigdb_fgsea_hs.rds")
  if (file.exists(cache_file)) {
    return(invisible(cache_file))
  }
  message("Downloading MSigDB gene sets (C2/C5/Reactome) via msigdbr...")
  if (!requireNamespace("msigdbr", quietly = TRUE)) {
    stop("msigdbr is required to download MSigDB gene sets; ensure install succeeded.")
  }
  fetch_set <- function(category, subcategory = NULL) {
    msigdbr::msigdbr(
      species = "Homo sapiens",
      category = category,
      subcategory = subcategory
    )
  }
  c2 <- fetch_set("C2")
  c5 <- fetch_set("C5")
  reactome <- fetch_set("C2", "CP:REACTOME")
  to_list <- function(df) {
    split(df$gene_symbol, paste(df$gs_name, df$gs_cat, df$gs_subcat, sep = "|"))
  }
  payload <- list(
    C2 = to_list(c2),
    C5 = to_list(c5),
    REACTOME = to_list(reactome)
  )
  saveRDS(payload, cache_file)
  message("Saved MSigDB gene sets to ", cache_file)
  invisible(cache_file)
}

script_dir <- tryCatch(normalizePath(dirname(sys.frame(1)$ofile %||% "."), winslash = "/", mustWork = FALSE), error = function(...) normalizePath(getwd(), winslash = "/", mustWork = FALSE))
repo_root <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = FALSE)
msigdb_cache_dir <- file.path(repo_root, ".dearmeta_cache", "msigdb")
tryCatch(
  download_msigdb_fgsea(msigdb_cache_dir),
  error = function(e) {
    warning("Failed to prefetch MSigDB gene sets: ", conditionMessage(e))
  }
)

message("All R dependencies installed.")
