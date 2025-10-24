#!/usr/bin/env Rscript

# DearMeta analysis pipeline

suppressPackageStartupMessages({
  library(optparse)
  library(jsonlite)
  library(data.table)
  library(minfi)
  library(sesame)
  library(sesameData)
  library(limma)
  library(sva)
  library(DMRcate)
  library(ggplot2)
  library(plotly)
  library(htmlwidgets)
  library(DT)
  library(RColorBrewer)
  library(VennDiagram)
  library(grid)
})

set.seed(1234)
options(SESAMEDATA_USE_ALT = TRUE)

# ---- Option parsing ---------------------------------------------------------

option_list <- list(
  make_option("--gse", type = "character", dest = "gse", help = "GEO series accession"),
  make_option("--project-root", type = "character", dest = "project_root", help = "Project directory for the GSE workspace"),
  make_option("--config", type = "character", dest = "config", help = "Path to configure.tsv"),
  make_option("--output-root", type = "character", dest = "output_root", help = "Root directory for outputs"),
  make_option("--min-group-size", type = "integer", default = 3, help = "Minimum samples per group"),
  make_option("--fdr-threshold", type = "double", default = 0.05, help = "Adjusted p-value threshold"),
  make_option("--delta-beta-threshold", type = "double", default = 0.05, help = "Absolute delta-beta threshold"),
  make_option("--top-n-cpgs", type = "integer", default = 10000, help = "Number of CpGs to retain for plots/tables")
)

opt <- parse_args(OptionParser(option_list = option_list))

mandatory <- c("gse", "project_root", "config", "output_root")
missing_opts <- mandatory[!mandatory %in% names(opt) | vapply(opt[mandatory], is.null, logical(1))]
if (length(missing_opts) > 0) {
  stop("Missing required arguments: ", paste(missing_opts, collapse = ", "))
}

if (is.null(opt$top_n_cpgs) || length(opt$top_n_cpgs) == 0) {
  opt$top_n_cpgs <- 10000L
} else {
  opt$top_n_cpgs <- as.integer(opt$top_n_cpgs)
  if (is.na(opt$top_n_cpgs) || opt$top_n_cpgs <= 0) {
    stop("--top-n-cpgs must be a positive integer")
  }
}

project_root <- normalizePath(opt$project_root, mustWork = TRUE)
config_path <- normalizePath(opt$config, mustWork = TRUE)
output_root <- normalizePath(opt$output_root, mustWork = TRUE)

# ---- Paths -----------------------------------------------------------------

paths <- list(
  preprocess = file.path(project_root, "02_preprocess"),
  analysis = file.path(project_root, "03_analysis"),
  figures = file.path(project_root, "04_figures"),
  interactive = file.path(project_root, "05_interactive"),
  runtime = file.path(project_root, "runtime")
)

for (dir_path in paths) {
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
  }
}

log_message <- function(...) {
  msg <- paste(...)
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg))
  flush.console()
  log_file <- file.path(paths$runtime, "pipeline.log")
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg), file = log_file, append = TRUE)
}

log_message("DearMeta analysis launched for", opt$gse)

# ---- Helper functions -------------------------------------------------------

read_configure <- function(path, project_root) {
  lines <- readLines(path)
  lines <- lines[!grepl("^\\s*#", lines)]
  if (length(lines) == 0) {
    stop("configure.tsv contains no data rows.")
  }
  cfg <- fread(text = paste(lines, collapse = "\n"), sep = "\t", na.strings = c("", "NA"))
  if (!"dear_group" %in% names(cfg)) {
    stop("configure.tsv must contain a dear_group column.")
  }
  cfg[, dear_group := trimws(as.character(dear_group))]
  cfg <- cfg[dear_group != ""]
  if (nrow(cfg) == 0) {
    stop("After filtering for dear_group assignments no samples remain.")
  }
  base_cols <- c("gsm_id", "sample_name", "idat_red", "idat_grn", "platform_version", "species")
  missing_base <- base_cols[!base_cols %in% names(cfg)]
  if (length(missing_base) > 0) {
    stop("configure.tsv is missing columns: ", paste(missing_base, collapse = ", "))
  }
  cfg[, idat_red := normalizePath(file.path(project_root, idat_red), mustWork = TRUE)]
  cfg[, idat_grn := normalizePath(file.path(project_root, idat_grn), mustWork = TRUE)]
  cfg
}

to_basename <- function(red_path, grn_path) {
  # Determine IDAT basename (without _Red.idat/_Grn.idat)
  pattern <- "_(Red|Grn)\\.idat$"
  if (grepl(pattern, basename(red_path), ignore.case = TRUE)) {
    sub(pattern, "", red_path, ignore.case = TRUE)
  } else if (grepl(pattern, basename(grn_path), ignore.case = TRUE)) {
    sub(pattern, "", grn_path, ignore.case = TRUE)
  } else {
    stop("Unable to derive IDAT basename from files: ", red_path, " and ", grn_path)
  }
}

detect_candidate_batches <- function(cfg) {
  ignore_cols <- c("dear_group", "gsm_id", "sample_name", "idat_red", "idat_grn", "platform_version", "species")
  candidates <- character(0)
  keywords <- c("slide", "barcode", "sentrix", "array", "plate", "batch", "center", "processing", "chip", "position")
  for (col in setdiff(names(cfg), ignore_cols)) {
    unique_ratio <- length(unique(na.omit(cfg[[col]]))) / nrow(cfg)
    if (any(grepl(paste(keywords, collapse = "|"), col, ignore.case = TRUE)) && unique_ratio < 0.9) {
      candidates <- c(candidates, col)
    } else if (unique_ratio > 0 && unique_ratio <= 0.2) {
      candidates <- c(candidates, col)
    }
  }
  unique(candidates)
}

split_covariates <- function(cfg, batch_cols) {
  ignore_cols <- c("dear_group", "gsm_id", "sample_name", "idat_red", "idat_grn", "platform_version", "species")
  numeric_cols <- character(0)
  factor_cols <- character(0)
  for (col in setdiff(names(cfg), ignore_cols)) {
    if (col %in% batch_cols) {
      next
    }
    if (is.numeric(cfg[[col]])) {
      numeric_cols <- c(numeric_cols, col)
    } else {
      factor_cols <- c(factor_cols, col)
    }
  }
  list(numeric = unique(numeric_cols), factor = unique(factor_cols))
}

safe_factor <- function(x) {
  x <- as.character(x)
  x[x == ""] <- NA_character_
  factor(x)
}

construct_design <- function(groups, numeric_covars, factor_covars, data, surrogate = NULL) {
  df <- data.table(group = groups)
  for (col in numeric_covars) {
    df[[col]] <- as.numeric(data[[col]])
  }
  for (col in factor_covars) {
    df[[col]] <- safe_factor(data[[col]])
  }
  if (!is.null(surrogate)) {
    df <- cbind(df, surrogate)
  }
  model.matrix(~ ., data = df)
}

save_matrix <- function(mat, path, compress = TRUE) {
  dt <- as.data.table(mat, keep.rownames = "probe_id")
  if (compress) {
    fwrite(dt, file = path, sep = "\t", compress = "gzip")
  } else {
    fwrite(dt, file = path, sep = "\t")
  }
}

write_plot <- function(plot_obj, filename, width = 7, height = 5) {
  tryCatch(
    {
      ggsave(filename = paste0(filename, ".png"), plot = plot_obj, width = width, height = height, dpi = 300)
      ggsave(filename = paste0(filename, ".pdf"), plot = plot_obj, width = width, height = height)
    },
    error = function(e) {
      log_message(sprintf("Skipping plot %s due to: %s", filename, conditionMessage(e)))
    }
  )
}

create_interactive <- function(plot_obj, filename) {
  html_path <- paste0(filename, ".html")
  tryCatch(
    {
      saveWidget(plot_obj, file = html_path, selfcontained = TRUE)
      html_path
    },
    error = function(e) {
      log_message(sprintf("Skipping interactive %s due to: %s", filename, conditionMessage(e)))
      NULL
    }
  )
}

# ---- Load configuration -----------------------------------------------------

config <- read_configure(config_path, project_root)
group_counts <- table(config$dear_group)
if (any(group_counts < opt$min_group_size)) {
  stop("Insufficient samples per group: ", paste(names(group_counts[group_counts < opt$min_group_size]), collapse = ", "))
}

log_message("Loaded", nrow(config), "samples for analysis across groups:", paste(names(group_counts), group_counts, sep = "=", collapse = "; "))

# Determine candidate batches and covariates
candidate_batches <- detect_candidate_batches(config)
covars <- split_covariates(config, candidate_batches)

platform_values <- unique(config$platform_version)
if (length(platform_values) != 1) {
  stop("Mixed platform_version values detected; DearMeta expects a single array type per run.")
}
platform_label <- platform_values[1]

array_type_map <- list(
  EPICv1 = list(dmr = "EPIC", sesame = "EPIC"),
  EPICv2 = list(dmr = "EPICv2", sesame = "EPICv2"),
  HM450 = list(dmr = "450K", sesame = "HM450.address")
)

platform_info <- array_type_map[[platform_label]]
if (is.null(platform_info)) {
  stop("Unsupported platform_version detected: ", platform_label)
}
log_message("Detected array platform:", platform_label)

# ---- Load IDATs with minfi --------------------------------------------------

log_message("Reading IDATs with minfi...")
basenames <- vapply(seq_len(nrow(config)), function(i) to_basename(config$idat_red[i], config$idat_grn[i]), character(1))
targets <- data.frame(
  Sample_Name = config$gsm_id,
  Sample_Group = config$dear_group,
  Basename = basenames,
  stringsAsFactors = FALSE
)
RGset <- read.metharray.exp(targets = targets, force = TRUE)
detP <- detectionP(RGset)
failed_probes <- rowMeans(detP > 0.01)
keep_probes <- failed_probes <= 0.1
log_message("minfi: retaining", sum(keep_probes), "probes out of", length(keep_probes))
RGset_filtered <- RGset[keep_probes, ]
MSet <- preprocessNoob(RGset_filtered)
beta_minfi <- getBeta(MSet)
if (is.null(rownames(beta_minfi))) {
  rownames(beta_minfi) <- featureNames(MSet)
}
M_minfi <- getM(MSet)

preprocess_dir <- paths$preprocess
saveRDS(MSet, file = file.path(preprocess_dir, "minfi_MSet.rds"))
save_matrix(detP, file.path(preprocess_dir, "detection_p_minfi.tsv.gz"))

# ---- Sesame pipeline --------------------------------------------------------

log_message("Running sesame pipeline...")
sesame_available <- TRUE
sesame_error <- NULL
sesame_betas <- list()
sesame_detp <- list()
if (!is.null(platform_info$sesame)) {
  tryCatch(
    {
      sesameDataGet(platform_info$sesame, verbose = TRUE)
      for (i in seq_len(nrow(config))) {
        base <- basenames[i]
        sig <- openSesame(base, prep = "Noob", platform = platform_label)
        sesame_betas[[config$gsm_id[i]]] <- getBetas(sig)
        sesame_detp[[config$gsm_id[i]]] <- sig$poobah
      }
    },
    error = function(e) {
      sesame_available <<- FALSE
      sesame_error <<- conditionMessage(e)
    }
  )
} else {
  sesame_available <- FALSE
  sesame_error <- "No sesame dataset configured for this platform."
}

if (sesame_available && length(sesame_betas) > 0) {
  beta_sesame <- do.call(cbind, sesame_betas)
  beta_sesame <- pmax(pmin(beta_sesame, 1 - 1e-6), 1e-6)
  colnames(beta_sesame) <- names(sesame_betas)
  sesame_poobah <- do.call(cbind, sesame_detp)
  colnames(sesame_poobah) <- names(sesame_detp)
  M_sesame <- log2(beta_sesame / (1 - beta_sesame))
  save_matrix(sesame_poobah, file.path(preprocess_dir, "pOOBAH_sesame.tsv.gz"))
} else {
  sesame_available <- FALSE
  beta_sesame <- NULL
  M_sesame <- NULL
  sesame_poobah <- NULL
  log_message("Skipping sesame pipeline:", if (!is.null(sesame_error)) sesame_error else "unknown error")
}

# ---- Batch effect assessment ------------------------------------------------

log_message("Assessing batch effects...")

assess_pca <- function(M_matrix, metadata, variables, prefix) {
  pcs <- prcomp(t(M_matrix), center = TRUE, scale. = TRUE)
  pc_count <- min(10, ncol(pcs$x))
  scores <- as.data.table(pcs$x[, seq_len(pc_count), drop = FALSE])
  scores[, sample := metadata$gsm_id]
  assessments <- list()
  pc_names <- grep("^PC", names(scores), value = TRUE)
  for (var in variables) {
    value <- metadata[[var]]
    if (length(unique(na.omit(value))) < 2) {
      next
    }
    var_type <- if (is.numeric(value)) "numeric" else "factor"
    pvals <- c()
    r2 <- c()
    for (pc in pc_names) {
      df <- data.frame(pc = scores[[pc]], variable = value)
      if (var_type == "numeric") {
        fit <- lm(pc ~ variable, data = df)
      } else {
        fit <- lm(pc ~ factor(variable), data = df)
      }
      an <- anova(fit)
      p <- an$`Pr(>F)`[1]
      pvals <- c(pvals, p)
      r2 <- c(r2, summary(fit)$r.squared)
    }
    assessments[[var]] <- data.table(
      variable = var,
      pc = pc_names[seq_len(length(pvals))],
      p_value = p.adjust(pvals, method = "BH"),
      r_squared = r2,
      type = var_type,
      pipeline = prefix
    )
  }
  if (length(assessments) > 0) {
    assessment_dt <- rbindlist(assessments)
  } else {
    assessment_dt <- data.table(variable = character(), pc = character(), p_value = numeric(), r_squared = numeric(), type = character(), pipeline = character())
  }
  list(scores = scores, assessments = assessment_dt)
}

empty_pca_result <- function() {
  list(
    scores = data.table(sample = character(), PC1 = numeric(), PC2 = numeric()),
    assessments = data.table(
      variable = character(),
      pc = character(),
      p_value = numeric(),
      r_squared = numeric(),
      type = character(),
      pipeline = character()
    )
  )
}

score_batch_effect <- function(pca_result, batch_vars) {
  stats <- pca_result$assessments[variable %in% batch_vars]
  if (nrow(stats) == 0) {
    return(NA_real_)
  }
  median(stats$p_value, na.rm = TRUE)
}

evaluate_correction <- function(method, M_matrix, apply_combat_fn, label) {
  if (method == "combat") {
    corrected <- apply_combat_fn(M_matrix)
  } else {
    corrected <- M_matrix
  }
  pca <- assess_pca(corrected, metadata_dt, variables_to_check, label)
  score <- score_batch_effect(pca, candidate_batches)
  list(name = method, M = corrected, pca = pca, score = score)
}

choose_best_method <- function(methods) {
  scores <- vapply(methods, function(x) x$score, numeric(1))
  scores_clean <- ifelse(is.na(scores), -Inf, scores)
  idx <- which.max(scores_clean)
  list(selected = methods[[idx]], all = methods)
}

metadata_dt <- as.data.table(config)
variables_to_check <- unique(c("dear_group", candidate_batches, covars$numeric, covars$factor))

pca_minfi_pre <- assess_pca(M_minfi, metadata_dt, variables_to_check, "minfi_pre")
pca_minfi <- pca_minfi_pre

if (sesame_available && !is.null(M_sesame)) {
  pca_sesame_pre <- assess_pca(M_sesame, metadata_dt, variables_to_check, "sesame_pre")
  pca_sesame <- pca_sesame_pre
} else {
  pca_sesame_pre <- empty_pca_result()
  pca_sesame <- pca_sesame_pre
}

should_run_combat <- FALSE
batch_to_use <- NULL
if (length(candidate_batches) > 0) {
  batch_stats <- pca_minfi$assessments[variable %in% candidate_batches]
  if (nrow(batch_stats) > 0) {
    significant <- batch_stats[p_value < 0.05 & r_squared >= 0.05]
    if (nrow(significant) > 0) {
      should_run_combat <- TRUE
      batch_to_use <- candidate_batches[1]
      log_message("ComBat will be applied using batch column", batch_to_use)
    }
  }
}

run_sva <- FALSE
if (!should_run_combat) {
  log_message("Evaluating need for SVA surrogate variables.")
  mod <- model.matrix(~ factor(dear_group), data = metadata_dt)
  n_sv <- tryCatch(num.sv(M_minfi, mod, method = "leek"), error = function(e) 0)
  if (!is.na(n_sv) && n_sv > 0) {
    run_sva <- TRUE
    log_message("SVA will add", n_sv, "surrogate variables.")
    sva_obj <- sva(M_minfi, mod, mod0 = model.matrix(~ 1, data = metadata_dt), n.sv = n_sv)
    surrogate_vars <- sva_obj$sv
  } else {
    surrogate_vars <- NULL
  }
} else {
  surrogate_vars <- NULL
}

apply_combat <- function(M_matrix, metadata, batch_col, design) {
  batch <- factor(metadata[[batch_col]])
  combat_res <- ComBat(dat = M_matrix, batch = batch, mod = design, par.prior = TRUE)
  rownames(combat_res) <- rownames(M_matrix)
  colnames(combat_res) <- colnames(M_matrix)
  combat_res
}

design_matrix <- construct_design(config$dear_group, covars$numeric, covars$factor, metadata_dt, surrogate_vars)

minfi_methods <- list(
  evaluate_correction("none", M_minfi, function(x) x, "minfi_none")
)
if (should_run_combat) {
  minfi_methods <- c(minfi_methods, list(evaluate_correction("combat", M_minfi, function(x) apply_combat(x, metadata_dt, batch_to_use, design_matrix), "minfi_combat")))
}
minfi_selection <- choose_best_method(minfi_methods)
selected_minfi_method <- minfi_selection$selected$name
M_minfi_corrected <- minfi_selection$selected$M
pca_minfi_post <- minfi_selection$selected$pca
minfi_score <- minfi_selection$selected$score
log_message("Selected batch correction for minfi: %s (median batch p-value: %s)", selected_minfi_method, ifelse(is.na(minfi_score), "NA", sprintf("%.4f", minfi_score)))

beta_from_M <- function(M_matrix) {
  beta <- 2^M_matrix / (1 + 2^M_matrix)
  beta <- pmax(pmin(beta, 1 - 1e-6), 1e-6)
  beta
}
beta_minfi <- beta_from_M(M_minfi_corrected)
save_matrix(beta_minfi, file.path(preprocess_dir, "beta_minfi.tsv.gz"))
save_matrix(M_minfi_corrected, file.path(preprocess_dir, "mvals_minfi.tsv.gz"))

if (sesame_available && !is.null(M_sesame)) {
  sesame_methods <- list(
    evaluate_correction("none", M_sesame, function(x) x, "sesame_none")
  )
  if (should_run_combat) {
    sesame_methods <- c(sesame_methods, list(evaluate_correction("combat", M_sesame, function(x) apply_combat(x, metadata_dt, batch_to_use, design_matrix), "sesame_combat")))
  }
  sesame_selection <- choose_best_method(sesame_methods)
  selected_sesame_method <- sesame_selection$selected$name
  M_sesame_corrected <- sesame_selection$selected$M
  pca_sesame_post <- sesame_selection$selected$pca
  sesame_score <- sesame_selection$selected$score
  log_message("Selected batch correction for sesame: %s (median batch p-value: %s)", selected_sesame_method, ifelse(is.na(sesame_score), "NA", sprintf("%.4f", sesame_score)))
  beta_sesame <- beta_from_M(M_sesame_corrected)
  save_matrix(beta_sesame, file.path(preprocess_dir, "beta_sesame.tsv.gz"))
  save_matrix(M_sesame_corrected, file.path(preprocess_dir, "mvals_sesame.tsv.gz"))
} else {
  selected_sesame_method <- if (sesame_available) "none" else "not_available"
  M_sesame_corrected <- NULL
  pca_sesame_post <- empty_pca_result()
  beta_sesame <- NULL
}

# ---- Differential methylation ----------------------------------------------

log_message("Running limma differential methylation...")
run_limma <- function(M_matrix, beta_matrix, metadata, design) {
  fit <- lmFit(M_matrix, design)
  fit <- eBayes(fit)
  top <- topTable(fit, coef = "groupdear_group", number = nrow(M_matrix), sort.by = "P")
  top <- data.table(probe_id = rownames(top), top)
  row_idx <- match(top$probe_id, rownames(beta_matrix))
  missing_idx <- is.na(row_idx)
  if (any(missing_idx)) {
    warning("Omitting ", sum(missing_idx), " probes absent from beta matrix during limma aggregation.")
    top <- top[!missing_idx]
    row_idx <- row_idx[!missing_idx]
  }
  beta_df <- beta_matrix[row_idx, , drop = FALSE]
  rownames(beta_df) <- top$probe_id
  group_levels <- levels(factor(metadata$dear_group))
  if (length(group_levels) != 2) {
    stop("Limma model currently supports two-group comparison.")
  }
  beta_group1 <- rowMeans(beta_df[, metadata$dear_group == group_levels[1], drop = FALSE], na.rm = TRUE)
  beta_group2 <- rowMeans(beta_df[, metadata$dear_group == group_levels[2], drop = FALSE], na.rm = TRUE)
  top[, delta_beta := beta_group2 - beta_group1]
  top[, direction := ifelse(delta_beta > 0, "hypermethylated", "hypomethylated")]
  top
}

design_for_limma <- construct_design(config$dear_group, covars$numeric, covars$factor, metadata_dt, surrogate_vars)
colnames(design_for_limma)[2] <- "groupdear_group"  # rename for clarity

limit_top_hits <- function(dt, n = 10000) {
  if (is.null(dt) || length(dt) == 0) {
    return(dt)
  }
  if (!inherits(dt, "data.frame") && !inherits(dt, "data.table")) {
    return(dt)
  }
  dt_names <- names(dt)
  if (is.null(dt_names) || !("P.Value" %in% dt_names)) {
    return(dt)
  }
  rows <- suppressWarnings(as.integer(tryCatch(nrow(dt), error = function(e) NA_integer_)))
  if (is.na(rows) || rows <= 0) {
    return(dt)
  }
  if (rows <= n) {
    return(dt)
  }
  dt <- dt[order(P.Value)]
  dt[seq_len(min(n, nrow(dt)))]
}

results_minfi <- run_limma(M_minfi_corrected, beta_minfi, metadata_dt, design_for_limma)
results_minfi <- limit_top_hits(results_minfi, opt$top_n_cpgs)
if (sesame_available && !is.null(M_sesame_corrected) && !is.null(beta_sesame)) {
  results_sesame <- run_limma(M_sesame_corrected, beta_sesame, metadata_dt, design_for_limma)
  results_sesame <- limit_top_hits(results_sesame, opt$top_n_cpgs)
} else {
  results_sesame <- results_minfi[0]
}

filter_significant <- function(res, fdr, delta_thresh) {
  res[adj.P.Val <= fdr & abs(delta_beta) >= delta_thresh]
}

sig_minfi <- filter_significant(results_minfi, opt$fdr_threshold, opt$delta_beta_threshold)
sig_sesame <- filter_significant(results_sesame, opt$fdr_threshold, opt$delta_beta_threshold)

merge_results <- function(minfi, sesame) {
  minfi[, pipeline := "minfi"]
  sesame[, pipeline := "sesame"]
  union <- rbindlist(list(minfi, sesame), fill = TRUE, use.names = TRUE)
  union[, shared_significant := probe_id %in% intersect(minfi$probe_id, sesame$probe_id)]
  intersection <- minfi[probe_id %in% sesame$probe_id]
  list(union = union, intersection = intersection)
}

integrated <- merge_results(sig_minfi, sig_sesame)
integrated$union <- limit_top_hits(integrated$union, opt$top_n_cpgs)
integrated$intersection <- limit_top_hits(integrated$intersection, opt$top_n_cpgs)

analysis_dir <- paths$analysis
fwrite(results_minfi, file.path(analysis_dir, "dmp_minfi.tsv.gz"), sep = "\t", compress = "gzip")
fwrite(results_sesame, file.path(analysis_dir, "dmp_sesame.tsv.gz"), sep = "\t", compress = "gzip")
fwrite(integrated$intersection, file.path(analysis_dir, "cpg_intersection.tsv.gz"), sep = "\t", compress = "gzip")
fwrite(integrated$union, file.path(analysis_dir, "cpg_union.tsv.gz"), sep = "\t", compress = "gzip")

# ---- DMRcate ----------------------------------------------------------------

log_message("Running DMRcate...")
array_type <- platform_info$dmr
if (nrow(sig_minfi) > 0) {
  annotation <- cpg.annotate(
    object = M_minfi_corrected,
    datatype = "array",
    what = "M",
    arraytype = array_type,
    analysis.type = "differential",
    design = design_for_limma,
    coef = "groupdear_group"
  )
  dmr_minfi <- dmrcate(annotation)
  dmr_ranges <- extractRanges(dmr_minfi, genome = "hg38")
} else {
  log_message("Skipping DMRcate for minfi: no significant CpGs detected at FDR threshold.")
  dmr_ranges <- data.frame(seqnames = character(), start = integer(), end = integer(), width = integer(), n.sites = integer(), meanstat = numeric())
}
fwrite(as.data.frame(dmr_ranges), file.path(analysis_dir, "dmr_minfi.tsv.gz"), sep = "\t", compress = "gzip")

if (sesame_available && !is.null(M_sesame_corrected)) {
  annotation_sesame <- cpg.annotate(
    object = M_sesame_corrected,
    datatype = "array",
    what = "M",
    arraytype = array_type,
    analysis.type = "differential",
    design = design_for_limma,
    coef = "groupdear_group"
  )
  dmr_sesame <- dmrcate(annotation_sesame)
  dmr_ranges_sesame <- extractRanges(dmr_sesame, genome = "hg38")
  fwrite(as.data.frame(dmr_ranges_sesame), file.path(analysis_dir, "dmr_sesame.tsv.gz"), sep = "\t", compress = "gzip")
} else {
  dmr_ranges_sesame <- data.frame(seqnames = character(), start = integer(), end = integer(), width = integer(), n.sites = integer(), meanstat = numeric())
  fwrite(dmr_ranges_sesame, file.path(analysis_dir, "dmr_sesame.tsv.gz"), sep = "\t", compress = "gzip")
}

# ---- Annotation -------------------------------------------------------------

log_message("Annotating CpGs...")
anno_df <- as.data.frame(getAnnotation(MSet))
if (!"probe_id" %in% names(anno_df)) {
  anno_df$probe_id <- rownames(anno_df)
}
annotation_data <- as.data.table(anno_df)
setnames(annotation_data, names(annotation_data), tolower(names(annotation_data)))
setkey(annotation_data, probe_id)

pick_column <- function(dt, options, default = NA_character_) {
  for (opt in options) {
    if (opt %in% names(dt)) {
      return(dt[[opt]])
    }
  }
  rep(default, nrow(dt))
}

annotations_clean <- data.table(
  probe_id = annotation_data$probe_id,
  chr = as.character(pick_column(annotation_data, c("chr", "chromosome", "seqnames"))),
  mapinfo = as.numeric(pick_column(annotation_data, c("mapinfo", "pos", "position"), default = NA_real_)),
  ucsc_refgene_name = pick_column(annotation_data, c("ucsc_refgene_name", "gene", "gene_symbol")),
  ucsc_refgene_group = pick_column(annotation_data, c("ucsc_refgene_group", "gene_group")),
  relation_to_island = pick_column(annotation_data, c("relation_to_island", "relationtoisland"))
)
format_gene_field <- function(x) {
  if (is.null(x) || all(is.na(x))) {
    return(rep(NA_character_, length(x)))
  }
  sapply(x, function(item) {
    if (is.na(item) || item == "") {
      return(NA_character_)
    }
    genes <- unique(trimws(strsplit(item, ";")[[1]]))
    genes <- genes[genes != ""]
    if (length(genes) == 0) NA_character_ else paste(genes, collapse = ";")
  }, USE.NAMES = FALSE)
}
annotations_clean[, gene_symbols := format_gene_field(ucsc_refgene_name)]
annotations_clean[, gene_regions := format_gene_field(ucsc_refgene_group)]

annotate_results <- function(res) {
  merge(res, annotations_clean, by = "probe_id", all.x = TRUE, sort = FALSE)
}
annotated_minfi <- annotate_results(results_minfi)
annotated_sesame <- annotate_results(results_sesame)
annotated_intersection <- annotate_results(integrated$intersection)
fwrite(annotated_minfi, file.path(analysis_dir, "annotated_cpg_minfi.tsv.gz"), sep = "\t", compress = "gzip")
fwrite(annotated_sesame, file.path(analysis_dir, "annotated_cpg_sesame.tsv.gz"), sep = "\t", compress = "gzip")
fwrite(annotated_intersection, file.path(analysis_dir, "annotated_cpg_intersection.tsv.gz"), sep = "\t", compress = "gzip")

# ---- Figures (static) -------------------------------------------------------

log_message("Generating static plots...")
fig_dir <- paths$figures

plot_pca <- function(pca_scores, metadata, color_var, title) {
  if (is.null(pca_scores) || nrow(pca_scores) == 0 || !"sample" %in% names(pca_scores)) {
    return(ggplot() + theme_void() + labs(title = paste(title, "(not available)")))
  }
  df <- merge(pca_scores, metadata, by.x = "sample", by.y = "gsm_id")
  df <- as.data.frame(df)
  ggplot(df, aes(x = PC1, y = PC2, color = .data[[color_var]])) +
    geom_point(size = 3) +
    theme_minimal() +
    labs(title = title, color = color_var)
}

pca_plots <- list(
  plot_pca(pca_minfi$scores, metadata_dt, "dear_group", "PCA minfi (pre)"),
  plot_pca(pca_minfi_post$scores, metadata_dt, "dear_group", "PCA minfi (post)"),
  plot_pca(pca_sesame$scores, metadata_dt, "dear_group", "PCA sesame (pre)"),
  plot_pca(pca_sesame_post$scores, metadata_dt, "dear_group", "PCA sesame (post)")
)
write_plot(pca_plots[[1]], file.path(fig_dir, "pca_minfi_pre"))
write_plot(pca_plots[[2]], file.path(fig_dir, "pca_minfi_post"))
write_plot(pca_plots[[3]], file.path(fig_dir, "pca_sesame_pre"))
write_plot(pca_plots[[4]], file.path(fig_dir, "pca_sesame_post"))

volcano_plot <- function(results, title) {
  df <- as.data.frame(results)
  if (!all(c("delta_beta", "adj.P.Val") %in% names(df))) {
    return(ggplot() + theme_void() + labs(title = paste(title, "(not available)")))
  }
  df$adj.P.Val <- as.numeric(df$adj.P.Val)
  df$delta_beta <- as.numeric(df$delta_beta)
  ggplot(df, aes(x = delta_beta, y = -log10(adj.P.Val), color = adj.P.Val <= opt$fdr_threshold)) +
    geom_point(alpha = 0.6) +
    scale_color_manual(values = c("grey60", "firebrick")) +
    theme_minimal() +
    labs(title = title, x = "Delta beta", y = "-log10(FDR)", color = paste0("FDR<=", opt$fdr_threshold))
}

write_plot(volcano_plot(results_minfi, "Volcano minfi"), file.path(fig_dir, "volcano_minfi"))
write_plot(volcano_plot(results_sesame, "Volcano sesame"), file.path(fig_dir, "volcano_sesame"))
write_plot(volcano_plot(integrated$intersection, "Volcano intersection"), file.path(fig_dir, "volcano_intersection"))

density_plot <- function(beta_matrix, title) {
  df <- melt(as.data.table(beta_matrix, keep.rownames = "probe_id"), id.vars = "probe_id", variable.name = "sample", value.name = "beta")
  df <- as.data.frame(df)
  ggplot(df, aes(beta, color = sample)) +
    geom_density() +
    theme_minimal() +
    labs(title = title, x = "Beta value")
}

write_plot(density_plot(beta_minfi, "Density betas minfi"), file.path(fig_dir, "density_betas_minfi"))
if (sesame_available && !is.null(beta_sesame)) {
  write_plot(density_plot(beta_sesame, "Density betas sesame"), file.path(fig_dir, "density_betas_sesame"))
}

boxplot_plot <- function(beta_matrix, title) {
  df <- melt(as.data.table(beta_matrix, keep.rownames = "probe_id"), id.vars = "probe_id", variable.name = "sample", value.name = "beta")
  df <- as.data.frame(df)
  ggplot(df, aes(x = sample, y = beta)) +
    geom_boxplot(outlier.size = 0.4) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    labs(title = title, x = "Sample", y = "Beta value")
}

write_plot(boxplot_plot(beta_minfi, "Beta boxplot minfi"), file.path(fig_dir, "boxplot_betas_minfi"))
if (sesame_available && !is.null(beta_sesame)) {
  write_plot(boxplot_plot(beta_sesame, "Beta boxplot sesame"), file.path(fig_dir, "boxplot_betas_sesame"))
}

manhattan_plot <- function(res, title) {
  df <- as.data.table(res)
  setnames(df, names(df), tolower(names(df)))
  required <- c("p.value", "mapinfo", "chr")
  if (!all(required %in% names(df))) {
    return(ggplot() + theme_void() + labs(title = paste(title, "(not available)")))
  }
  df <- df[!is.na(p.value) & !is.na(mapinfo) & !is.na(chr)]
  if (nrow(df) == 0) {
    return(ggplot() + theme_void() + labs(title = paste(title, "(insufficient data)")))
  }
  df[, chr_num := as.numeric(gsub("chr", "", chr, ignore.case = TRUE))]
  df[is.na(chr_num) & grepl("x", chr, ignore.case = TRUE), chr_num := 23]
  df[is.na(chr_num) & grepl("y", chr, ignore.case = TRUE), chr_num := 24]
  df[is.na(chr_num), chr_num := 25]
  setorder(df, chr_num, mapinfo)
  df[, index := .I]
  df[, neg_log10_p := -log10(p.value)]
  df[, tooltip := paste0(
    "CpG: ", probe_id,
    ifelse(!is.null(gene_symbols) & !is.na(gene_symbols), paste0("<br>Gene: ", gene_symbols), ""),
    ifelse(!is.na(p.value), sprintf("<br>p-value: %.4g", p.value), ""),
    ifelse("adj.p.val" %in% names(df) && !is.na(adj.p.val), sprintf("<br>adj.p-value: %.4g", adj.p.val), ""),
    ifelse("delta_beta" %in% names(df) && !is.na(delta_beta), sprintf("<br>delta-beta: %.3f", delta_beta), ""),
    ifelse("logfc" %in% names(df) && !is.na(logfc), sprintf("<br>logFC: %.3f", logfc), "")
  )]
  centers <- df[, .(center = mean(index)), by = chr_num]
  ggplot(df, aes(x = index, y = neg_log10_p, color = factor(chr_num %% 2), text = tooltip)) +
    geom_point(alpha = 0.6, size = 0.7) +
    scale_x_continuous(breaks = centers$center, labels = centers$chr_num) +
    scale_color_manual(values = c("#1f77b4", "#ff7f0e")) +
    theme_minimal() +
    theme(legend.position = "none") +
    labs(title = title, x = "Chromosome", y = "-log10(p-value)")
}

write_plot(manhattan_plot(annotated_minfi, "Manhattan minfi"), file.path(fig_dir, "manhattan_minfi"))
write_plot(manhattan_plot(annotated_sesame, "Manhattan sesame"), file.path(fig_dir, "manhattan_sesame"))
write_plot(manhattan_plot(annotated_intersection, "Manhattan intersection"), file.path(fig_dir, "manhattan_intersection"))

qq_plot <- function(res, title) {
  df <- as.data.frame(res)
  colnames(df) <- tolower(colnames(df))
  if (!"p.value" %in% names(df)) {
    return(ggplot() + theme_void() + labs(title = paste(title, "(not available)")))
  }
  observed <- sort(df$p.value)
  expected <- -log10(ppoints(length(observed)))
  df <- data.table(expected = expected, observed = -log10(observed))
  ggplot(df, aes(expected, observed)) +
    geom_line(color = "steelblue") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    theme_minimal() +
    labs(title = title, x = "Expected -log10(p)", y = "Observed -log10(p)")
}

write_plot(qq_plot(results_minfi, "QQ minfi"), file.path(fig_dir, "qq_minfi"))
write_plot(qq_plot(results_sesame, "QQ sesame"), file.path(fig_dir, "qq_sesame"))

save_venn <- function(sets, filename_base) {
  venn <- venn.diagram(
    sets,
    filename = NULL,
    fill = c("#3b8bba", "#e24a33"),
    alpha = 0.5,
    cat.cex = 1.2,
    cex = 1.4,
    lty = "dashed"
  )
  png(paste0(filename_base, ".png"), width = 1600, height = 1600, res = 220)
  grid.newpage()
  grid.draw(venn)
  dev.off()
  pdf(paste0(filename_base, ".pdf"), width = 8, height = 8)
  grid.newpage()
  grid.draw(venn)
  dev.off()
}

if (nrow(sig_minfi) > 0 || nrow(sig_sesame) > 0) {
  save_venn(list(minfi = sig_minfi$probe_id, sesame = sig_sesame$probe_id), file.path(fig_dir, "venn_minfi_sesame"))
}

# ---- Interactive outputs ----------------------------------------------------

log_message("Creating interactive visuals...")
interactive_dir <- paths$interactive

interactive_files <- list()
path <- create_interactive(ggplotly(pca_plots[[1]]), file.path(interactive_dir, "pca_pre"))
if (!is.null(path)) interactive_files$pca_pre <- path
path <- create_interactive(ggplotly(pca_plots[[2]]), file.path(interactive_dir, "pca_post"))
if (!is.null(path)) interactive_files$pca_post <- path
path <- create_interactive(ggplotly(volcano_plot(annotated_minfi, "Volcano minfi")), file.path(interactive_dir, "volcano_minfi"))
if (!is.null(path)) interactive_files$volcano_minfi <- path
path <- create_interactive(ggplotly(volcano_plot(integrated$intersection, "Volcano intersection")), file.path(interactive_dir, "volcano_intersection"))
if (!is.null(path)) interactive_files$volcano_intersection <- path
path <- create_interactive(ggplotly(manhattan_plot(annotated_minfi, "Manhattan minfi")), file.path(interactive_dir, "manhattan_minfi"))
if (!is.null(path)) interactive_files$manhattan_minfi <- path
path <- create_interactive(ggplotly(manhattan_plot(annotated_intersection, "Manhattan intersection")), file.path(interactive_dir, "manhattan_intersection"))
if (!is.null(path)) interactive_files$manhattan_intersection <- path

datatable_widget <- datatable(as.data.frame(annotated_minfi), options = list(pageLength = 25, scrollX = TRUE))
path <- create_interactive(datatable_widget, file.path(interactive_dir, "table_minfi"))
if (!is.null(path)) interactive_files$table_minfi <- path
if (sesame_available && nrow(results_sesame) > 0) {
  path <- create_interactive(ggplotly(volcano_plot(annotated_sesame, "Volcano sesame")), file.path(interactive_dir, "volcano_sesame"))
  if (!is.null(path)) interactive_files$volcano_sesame <- path
  path <- create_interactive(ggplotly(manhattan_plot(annotated_sesame, "Manhattan sesame")), file.path(interactive_dir, "manhattan_sesame"))
  if (!is.null(path)) interactive_files$manhattan_sesame <- path
  datatable_widget2 <- datatable(as.data.frame(annotated_sesame), options = list(pageLength = 25, scrollX = TRUE))
  path <- create_interactive(datatable_widget2, file.path(interactive_dir, "table_sesame"))
  if (!is.null(path)) interactive_files$table_sesame <- path
}
datatable_widget3 <- datatable(as.data.frame(annotated_intersection), options = list(pageLength = 25, scrollX = TRUE))
path <- create_interactive(datatable_widget3, file.path(interactive_dir, "table_intersection"))
if (!is.null(path)) interactive_files$table_intersection <- path

index_path <- file.path(interactive_dir, "index.html")
index_html <- "<html><head><title>DearMeta Interactive</title></head><body><h1>Interactive Outputs</h1><ul>"
for (name in names(interactive_files)) {
  rel <- basename(interactive_files[[name]])
  index_html <- paste0(index_html, sprintf("<li><a href=\"%s\">%s</a></li>", rel, rel))
}
index_html <- paste0(index_html, "</ul></body></html>")
writeLines(index_html, con = index_path)

# ---- Summary ----------------------------------------------------------------

manifest_entries <- function(dir_path) {
  files <- list.files(dir_path, full.names = TRUE)
  lapply(files, function(p) {
    info <- file.info(p)
    list(path = p, size = unname(info$size))
  })
}

top_cpg_table <- function(dt) {
  if (nrow(dt) == 0) {
    return(list())
  }
  dt_sorted <- dt[order(adj.P.Val)]
  if (nrow(dt_sorted) > 10) {
    dt_sorted <- dt_sorted[1:10]
  }
  preferred_cols <- c("probe_id", "logFC", "delta_beta", "adj.P.Val", "P.Value", "chr", "mapinfo", "gene_symbols", "gene_regions", "relation_to_island")
  cols <- intersect(preferred_cols, names(dt_sorted))
  df <- as.data.frame(dt_sorted[, ..cols])
  split(df, seq_len(nrow(df)))
}

top_cpgs <- list(
  minfi = top_cpg_table(annotated_minfi),
  sesame = top_cpg_table(annotated_sesame),
  intersection = top_cpg_table(annotated_intersection)
)

metric_entry <- function(dt) {
  if (is.null(dt) || nrow(dt) == 0) {
    return(list())
  }
  as.list(dt[1])
}

batch_metrics <- list(
  minfi_pre = metric_entry(pca_minfi$assessments[variable %in% candidate_batches, .(median_p = median(p_value, na.rm = TRUE), median_r2 = median(r_squared, na.rm = TRUE))]),
  minfi_post = metric_entry(pca_minfi_post$assessments[variable %in% candidate_batches, .(median_p = median(p_value, na.rm = TRUE), median_r2 = median(r_squared, na.rm = TRUE))]),
  sesame_pre = metric_entry(pca_sesame$assessments[variable %in% candidate_batches, .(median_p = median(p_value, na.rm = TRUE), median_r2 = median(r_squared, na.rm = TRUE))]),
  sesame_post = metric_entry(pca_sesame_post$assessments[variable %in% candidate_batches, .(median_p = median(p_value, na.rm = TRUE), median_r2 = median(r_squared, na.rm = TRUE))])
)

surrogate_count <- if (!is.null(surrogate_vars)) ncol(as.matrix(surrogate_vars)) else 0

summary <- list(
  gse = opt$gse,
  samples = nrow(config),
  platform = platform_label,
  array_type = platform_info$dmr,
  groups = as.list(group_counts),
  candidate_batches = candidate_batches,
  covariates = covars,
  batch_methods = list(minfi = selected_minfi_method, sesame = selected_sesame_method),
  combat_applied = should_run_combat,
  batch_column = batch_to_use,
  sva_surrogates = surrogate_count,
  top_n_cpgs = opt$top_n_cpgs,
  significant_cpgs = list(
    minfi = nrow(sig_minfi),
    sesame = nrow(sig_sesame),
    intersection = nrow(integrated$intersection)
  ),
  top_cpgs = top_cpgs,
  batch_metrics = batch_metrics,
  files = list(
    preprocess = manifest_entries(paths$preprocess),
    analysis = manifest_entries(paths$analysis),
    figures = manifest_entries(paths$figures),
    interactive = manifest_entries(paths$interactive)
  )
)
write_json(summary, file.path(paths$runtime, "analysis_summary.json"), pretty = TRUE, auto_unbox = TRUE)

log_message("DearMeta analysis completed successfully.")
