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
  library(htmltools)
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
  make_option("--min-group-size", type = "integer", dest = "min_group_size", default = 2, help = "Minimum samples per group"),
  make_option("--fdr-threshold", type = "double", dest = "fdr_threshold", default = 0.05, help = "Adjusted p-value threshold"),
  make_option("--delta-beta-threshold", type = "double", dest = "delta_beta_threshold", default = 0.05, help = "Absolute delta-beta threshold"),
  make_option("--top-n-cpgs", type = "integer", dest = "top_n_cpgs", default = 10000, help = "Number of CpGs to retain for plots/tables")
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

options(dearmeta.runtime_dir = paths$runtime)

log_message <- function(...) {
  args <- list(...)
  msg <- NULL
  if (length(args) >= 1 && is.character(args[[1]]) && grepl("%", args[[1]], fixed = FALSE)) {
    msg <- tryCatch(do.call(sprintf, args), error = function(e) NULL)
  }
  if (is.null(msg)) {
    msg <- paste(vapply(args, function(x) {
      if (is.null(x)) {
        "NULL"
      } else if (is.atomic(x) && length(x) == 1) {
        as.character(x)
      } else {
        paste(capture.output(print(x)), collapse = " ")
      }
    }, character(1)), collapse = " ")
  }
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg))
  flush.console()
  log_file <- file.path(paths$runtime, "pipeline.log")
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg), file = log_file, append = TRUE)
}

`%||%` <- function(x, y) {
  if (is.null(x)) {
    return(y)
  }
  if (is.character(x) && length(x) == 1 && identical(trimws(x), "")) {
    return(y)
  }
  if (length(x) == 0) {
    return(y)
  }
  x
}

DEFAULT_PROTECTED_BATCH <- c("sex", "gender", "sex_at_birth", "biological_sex")

split_directive_values <- function(value) {
  parts <- trimws(unlist(strsplit(value, "[,;|]+")))
  parts <- parts[parts != ""]
  unique(parts)
}

parse_config_directives <- function(directive_lines) {
  directives <- list(
    batch = character(),
    protect = character(),
    drop = character(),
    numeric = character(),
    factor = character()
  )
  if (length(directive_lines) == 0) {
    return(directives)
  }
  pattern <- "^#\\s*dear[_-]?([a-zA-Z0-9_]+)\\s*:\\s*(.+)$"
  for (line in directive_lines) {
    matches <- regexec(pattern, line)
    capture <- regmatches(line, matches)[[1]]
    if (length(capture) == 3) {
      key <- tolower(capture[2])
      values <- split_directive_values(capture[3])
      if (length(values) == 0) next
      if (key %in% c("batch", "batches", "batch_columns")) {
        directives$batch <- unique(c(directives$batch, values))
      } else if (key %in% c("protect", "protect_columns", "avoid_batch")) {
        directives$protect <- unique(c(directives$protect, values))
      } else if (key %in% c("drop", "drop_columns", "exclude")) {
        directives$drop <- unique(c(directives$drop, values))
      } else if (key %in% c("numeric_covariates", "covariates_numeric", "numeric")) {
        directives$numeric <- unique(c(directives$numeric, values))
      } else if (key %in% c("factor_covariates", "covariates_factor", "factor", "categorical")) {
        directives$factor <- unique(c(directives$factor, values))
      }
    }
  }
  directives
}

match_columns <- function(cfg_names, requested) {
  if (length(requested) == 0) {
    return(character())
  }
  matched <- cfg_names[tolower(cfg_names) %in% tolower(requested)]
  unique(matched)
}

log_message("DearMeta analysis launched for", opt$gse)

# ---- Helper functions -------------------------------------------------------

read_configure <- function(path, project_root) {
  raw_lines <- readLines(path)
  directive_lines <- raw_lines[grepl("^\\s*#", raw_lines)]
  data_lines <- raw_lines[!grepl("^\\s*#", raw_lines)]
  directives <- parse_config_directives(directive_lines)
  if (length(data_lines) == 0) {
    stop("configure.tsv contains no data rows.")
  }
  cfg <- fread(text = paste(data_lines, collapse = "\n"), sep = "\t", na.strings = c("", "NA"))
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
  drop_cols <- match_columns(names(cfg), directives$drop)
  if (length(drop_cols) > 0) {
    invalid <- intersect(drop_cols, base_cols)
    if (length(invalid) > 0) {
      log_message("Cannot drop required columns: %s", paste(invalid, collapse = ", "))
      drop_cols <- setdiff(drop_cols, invalid)
    }
    if (length(drop_cols) > 0) {
      cfg[, (drop_cols) := NULL]
      log_message("Dropped columns per directives: %s", paste(drop_cols, collapse = ", "))
    }
  }
  cfg[, idat_red := normalizePath(file.path(project_root, idat_red), mustWork = TRUE)]
  cfg[, idat_grn := normalizePath(file.path(project_root, idat_grn), mustWork = TRUE)]
  attr(cfg, "directives") <- directives
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

load_sesame_manifest <- function(candidates, verbose = FALSE) {
  entries <- unique(as.character(unlist(candidates, use.names = FALSE)))
  entries <- entries[nzchar(entries)]
  attempts <- list()
  if (length(entries) == 0) {
    return(list(ordering = NULL, controls = NULL, source = NULL, attempts = attempts))
  }
  for (entry in entries) {
    result <- tryCatch(
      sesameDataGet(entry, verbose = verbose),
      error = function(e) {
        attempts[[entry]] <<- conditionMessage(e)
        NULL
      }
    )
    if (is.list(result) && all(c("ordering", "controls") %in% names(result))) {
      return(list(ordering = result$ordering, controls = result$controls, source = entry, attempts = attempts))
    }
    if (!is.null(result) && is.null(attempts[[entry]])) {
      attempts[[entry]] <- "Dataset lacks ordering/controls manifest."
    }
  }
  list(ordering = NULL, controls = NULL, source = NULL, attempts = attempts)
}

sanitize_covariate <- function(x) {
  if (is.character(x)) {
    x <- trimws(x)
    x[x == ""] <- NA_character_
  }
  x
}

has_useful_variation <- function(x, min_non_missing = 2) {
  x <- sanitize_covariate(x)
  non_missing <- sum(!is.na(x))
  if (non_missing < min_non_missing) {
    return(FALSE)
  }
  length(unique(x[!is.na(x)])) >= 2
}

duplicates_dear_group <- function(group_values, candidate_values) {
  groups <- sanitize_covariate(group_values)
  values <- sanitize_covariate(candidate_values)
  shared_idx <- !is.na(values) & !is.na(groups)
  if (sum(shared_idx) < 2) {
    return(FALSE)
  }
  shared_groups <- groups[shared_idx]
  shared_values <- values[shared_idx]
  mapping_value_to_group <- tapply(shared_groups, shared_values, function(x) length(unique(x)))
  if (!is.null(mapping_value_to_group) &&
    all(mapping_value_to_group == 1) &&
    length(mapping_value_to_group) <= length(unique(na.omit(shared_groups)))) {
    return(TRUE)
  }
  mapping_group_to_value <- tapply(shared_values, shared_groups, function(x) length(unique(x)))
  if (!is.null(mapping_group_to_value) &&
    all(mapping_group_to_value == 1) &&
    length(mapping_group_to_value) <= length(unique(na.omit(shared_values)))) {
    return(TRUE)
  }
  FALSE
}

remove_group_duplicates <- function(columns, cfg, group_col = "dear_group", log_prefix = "") {
  if (length(columns) == 0) {
    return(columns)
  }
  keep <- logical(length(columns))
  for (i in seq_along(columns)) {
    col <- columns[i]
    if (!col %in% names(cfg)) {
      keep[i] <- FALSE
      next
    }
    if (duplicates_dear_group(cfg[[group_col]], cfg[[col]])) {
      if (nzchar(log_prefix)) {
        log_message("%s column %s mirrors dear_group and will be excluded.", log_prefix, col)
      }
      keep[i] <- FALSE
    } else {
      keep[i] <- TRUE
    }
  }
  columns[keep]
}

detect_candidate_batches <- function(cfg, protect_cols = character()) {
  ignore_cols <- c("dear_group", "gsm_id", "sample_name", "idat_red", "idat_grn", "platform_version", "species")
  candidates <- character(0)
  keywords <- c("slide", "barcode", "sentrix", "array", "plate", "batch", "center", "processing", "chip", "position")
  protected_lower <- tolower(protect_cols)
  for (col in setdiff(names(cfg), ignore_cols)) {
    if (tolower(col) %in% protected_lower) {
      next
    }
    values <- sanitize_covariate(cfg[[col]])
    if (!has_useful_variation(values)) {
      next
    }
    unique_ratio <- length(unique(na.omit(values))) / nrow(cfg)
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
  dropped <- list()
  group_values <- sanitize_covariate(cfg$dear_group)
  group_levels <- unique(na.omit(group_values))
  for (col in setdiff(names(cfg), ignore_cols)) {
    if (col %in% batch_cols) {
      next
    }
    values <- sanitize_covariate(cfg[[col]])
    non_missing <- sum(!is.na(values))
    if (non_missing < 2) {
      dropped[[col]] <- "less than two observations"
      next
    }
    if (!has_useful_variation(values)) {
      dropped[[col]] <- "no variation after cleaning"
      next
    }
    if (duplicates_dear_group(group_values, values)) {
      dropped[[col]] <- "duplicates dear_group"
      next
    }
    if (is.numeric(cfg[[col]])) {
      numeric_cols <- c(numeric_cols, col)
    } else if (is.logical(cfg[[col]])) {
      factor_cols <- c(factor_cols, col)
    } else if (is.character(cfg[[col]])) {
      numeric_candidate <- suppressWarnings(as.numeric(values))
      if (sum(!is.na(numeric_candidate)) >= 2 && length(unique(na.omit(numeric_candidate))) >= 2) {
        numeric_cols <- c(numeric_cols, col)
      } else {
        factor_cols <- c(factor_cols, col)
      }
    } else {
      factor_cols <- c(factor_cols, col)
    }
  }
  result <- list(numeric = unique(numeric_cols), factor = unique(factor_cols))
  attr(result, "dropped_covariates") <- dropped
  result
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

beta_from_M <- function(M_matrix) {
  beta <- 2^M_matrix / (1 + 2^M_matrix)
  beta <- pmax(pmin(beta, 1 - 1e-6), 1e-6)
  beta
}

apply_combat <- function(M_matrix, metadata, batch_col, design) {
  if (is.null(batch_col) || !nzchar(batch_col)) {
    stop("apply_combat requires a batch column")
  }
  batch <- safe_factor(metadata[[batch_col]])
  if (length(levels(batch)) < 2) {
    stop(sprintf("batch column %s has fewer than two levels", batch_col))
  }
  design_rank <- qr(design)$rank
  if (design_rank < ncol(design)) {
    stop(sprintf("design matrix is rank deficient (rank %s < %s)", design_rank, ncol(design)))
  }
  combat_res <- ComBat(dat = M_matrix, batch = batch, mod = design, par.prior = TRUE)
  rownames(combat_res) <- rownames(M_matrix)
  colnames(combat_res) <- colnames(M_matrix)
  combat_res
}

covariate_association_stats <- function(metadata, group_col, numeric_covars = character(), factor_covars = character()) {
  stats <- list()
  group <- safe_factor(metadata[[group_col]])
  for (col in numeric_covars) {
    values <- metadata[[col]]
    if (all(is.na(values))) {
      next
    }
    df <- data.frame(value = as.numeric(values), group = group)
    df <- df[complete.cases(df), , drop = FALSE]
    if (nrow(df) < 3 || length(unique(df$value)) < 2) {
      next
    }
    fit <- tryCatch(lm(value ~ group, data = df), error = function(e) NULL)
    if (is.null(fit)) {
      next
    }
    an <- tryCatch(anova(fit), error = function(e) NULL)
    if (is.null(an)) {
      next
    }
    p_val <- suppressWarnings(an$`Pr(>F)`[1])
    if (is.na(p_val)) {
      p_val <- 1
    }
    stats[[length(stats) + 1L]] <- data.table(
      name = col,
      type = "numeric",
      p_value = p_val,
      score = -log10(p_val + 1e-12)
    )
  }
  for (col in factor_covars) {
    values <- safe_factor(metadata[[col]])
    if (length(levels(values)) < 2) {
      next
    }
    tbl <- table(group, values)
    if (!all(dim(tbl) > 0)) {
      next
    }
    test <- tryCatch(suppressWarnings(chisq.test(tbl)), error = function(e) NULL)
    if (is.null(test) || !is.finite(test$p.value)) {
      test <- tryCatch(fisher.test(tbl), error = function(e) NULL)
    }
    p_val <- if (!is.null(test) && is.finite(test$p.value)) test$p.value else 1
    stats[[length(stats) + 1L]] <- data.table(
      name = col,
      type = "factor",
      p_value = p_val,
      score = -log10(p_val + 1e-12)
    )
  }
  if (length(stats) == 0) {
    return(data.table(name = character(), type = character(), p_value = numeric(), score = numeric()))
  }
  rbindlist(stats)
}

build_combos <- function(values, max_size) {
  if (length(values) == 0 || max_size <= 0) {
    return(list(character()))
  }
  combos <- list(character())
  max_k <- min(length(values), max_size)
  for (k in seq_len(max_k)) {
    if (k == 0) {
      next
    }
    cmb <- combn(values, k, simplify = FALSE)
    combos <- c(combos, cmb)
  }
  combos
}

generate_covariate_sets <- function(required_numeric, required_factor, optional_numeric_stats, optional_factor_stats, max_optional = 3, max_combo = 2) {
  optional_numeric <- character()
  optional_factor <- character()
  if (nrow(optional_numeric_stats) > 0) {
    optional_numeric <- optional_numeric_stats[order(-score, name)]$name
    optional_numeric <- head(optional_numeric, max_optional)
  }
  if (nrow(optional_factor_stats) > 0) {
    optional_factor <- optional_factor_stats[order(-score, name)]$name
    optional_factor <- head(optional_factor, max_optional)
  }
  numeric_combos <- build_combos(optional_numeric, max_combo)
  factor_combos <- build_combos(optional_factor, max_combo)
  sets <- list()
  seen <- character()
  for (num in numeric_combos) {
    for (fac in factor_combos) {
      numeric_covars <- unique(c(required_numeric, num))
      factor_covars <- unique(c(required_factor, fac))
      key <- paste(
        paste(sort(numeric_covars), collapse = "|"),
        paste(sort(factor_covars), collapse = "|"),
        sep = "||"
      )
      if (key %in% seen) {
        next
      }
      seen <- c(seen, key)
      sets[[length(sets) + 1L]] <- list(
        numeric = numeric_covars,
        factor = factor_covars,
        optional_numeric = setdiff(numeric_covars, required_numeric),
        optional_factor = setdiff(factor_covars, required_factor)
      )
    }
  }
  if (length(sets) == 0) {
    sets <- list(list(
      numeric = unique(required_numeric),
      factor = unique(required_factor),
      optional_numeric = character(),
      optional_factor = character()
    ))
  }
  for (i in seq_along(sets)) {
    sets[[i]]$id <- sprintf("set_%03d", i)
  }
  sets
}

build_batch_options <- function(metadata, group_col, manual_batches, candidate_batches, max_auto = 3) {
  manual <- unique(manual_batches)
  auto <- setdiff(candidate_batches, manual)
  is_dup <- function(col) {
    if (is.null(col) || is.na(col) || !nzchar(col) || !col %in% names(metadata)) {
      return(TRUE)
    }
    duplicates_dear_group(metadata[[group_col]], metadata[[col]])
  }
  if (length(manual) > 0) {
    manual <- manual[!vapply(manual, is_dup, logical(1))]
  }
  if (length(auto) > 0) {
    auto <- auto[!vapply(auto, is_dup, logical(1))]
  }
  stats <- covariate_association_stats(metadata, group_col, factor_covars = auto)
  auto_ranked <- if (nrow(stats) > 0) stats[order(-score, name)]$name else character()
  auto_limited <- head(auto_ranked, max_auto)
  unique(c(NA_character_, manual, auto_limited))
}

format_model_id <- function(set_id, batch_col, use_combat, use_sva) {
  batch_label <- if (is.null(batch_col) || !nzchar(batch_col)) "none" else batch_col
  sprintf("%s|batch=%s|combat=%s|sva=%s", set_id, batch_label, use_combat, use_sva)
}

execute_model_configuration <- function(params, data, opt, candidate_batches, log_fn = NULL, keep_outputs = FALSE) {
  metadata <- data$metadata
  M_minfi <- data$M_minfi
  M_sesame <- data$M_sesame
  sesame_available <- isTRUE(data$sesame_available)
  stage <- "initialisation"
  result <- tryCatch({
    stage <- "resolve_covariates"
    factor_final <- unique(c(params$factor_covars, if (!params$use_combat && !is.null(params$batch_col)) params$batch_col else NULL))
    numeric_final <- unique(params$numeric_covars)

    old_model_id <- getOption("dearmeta.current_model_id")
    options(dearmeta.current_model_id = params$model_id %||% "<unnamed>")
    on.exit(options(dearmeta.current_model_id = old_model_id), add = TRUE)

    stage <- "design_base"
    design_base <- construct_design(metadata$dear_group, numeric_final, factor_final, metadata, surrogate = NULL)
    if (qr(design_base)$rank < ncol(design_base)) {
      stop("design matrix is rank deficient")
    }

    stage <- "sva_setup"
    surrogate_vars <- NULL
    n_sv <- 0L
    if (isTRUE(params$use_sva)) {
      mod_full <- design_base
      mod0 <- model.matrix(~ 1, data = metadata)
      n_sv <- tryCatch(num.sv(M_minfi, mod_full, method = "leek"), error = function(e) NA_integer_)
      if (is.na(n_sv) || n_sv <= 0) {
        n_sv <- 0L
      } else {
        stage <- "compute_sva"
        sva_obj <- sva(M_minfi, mod_full, mod0 = mod0, n.sv = n_sv)
        surrogate_vars <- sva_obj$sv
        surrogate_vars <- as.matrix(surrogate_vars)
        colnames(surrogate_vars) <- paste0("SV", seq_len(ncol(surrogate_vars)))
      }
    }

    stage <- "apply_sva_minfi"
    M_minfi_sva <- if (!is.null(surrogate_vars) && ncol(surrogate_vars) > 0) {
      removeBatchEffect(M_minfi, covariates = surrogate_vars, design = design_base)
    } else {
      M_minfi
    }

    stage <- "apply_sva_sesame"
    M_sesame_sva <- if (sesame_available && !is.null(M_sesame)) {
      if (!is.null(surrogate_vars) && ncol(surrogate_vars) > 0) {
        removeBatchEffect(M_sesame, covariates = surrogate_vars, design = design_base)
      } else {
        M_sesame
      }
    } else {
      NULL
    }

    stage <- "design_with_sv"
    design_with_sv <- construct_design(metadata$dear_group, numeric_final, factor_final, metadata, surrogate = surrogate_vars)
    if (qr(design_with_sv)$rank < ncol(design_with_sv)) {
      stop("design+SV matrix is rank deficient")
    }

    stage <- "apply_combat"
    M_minfi_corrected <- M_minfi_sva
    M_sesame_corrected <- M_sesame_sva
    if (!is.null(params$batch_col)) {
      levels_count <- length(unique(na.omit(metadata[[params$batch_col]])))
      if (levels_count < 2) {
        stop(sprintf("batch column %s has fewer than two levels", params$batch_col))
      }
      if (isTRUE(params$use_combat)) {
        M_minfi_corrected <- apply_combat(M_minfi_sva, metadata, params$batch_col, design_with_sv)
        if (sesame_available && !is.null(M_sesame_sva)) {
          M_sesame_corrected <- apply_combat(M_sesame_sva, metadata, params$batch_col, design_with_sv)
        }
      }
    } else if (isTRUE(params$use_combat)) {
      stop("ComBat requested without a batch column")
    }

    stage <- "pca_assessment"
    variables_to_check <- unique(c("dear_group", candidate_batches, numeric_final, factor_final, params$batch_col))
    pca_minfi_post <- assess_pca(M_minfi_corrected, metadata, variables_to_check, "minfi_post")
    pca_sesame_post <- if (sesame_available && !is.null(M_sesame_corrected)) {
      assess_pca(M_sesame_corrected, metadata, variables_to_check, "sesame_post")
    } else {
      empty_pca_result()
    }

    stage <- "beta_conversion"
    beta_minfi_corrected <- beta_from_M(M_minfi_corrected)
    beta_sesame_corrected <- if (sesame_available && !is.null(M_sesame_corrected)) beta_from_M(M_sesame_corrected) else NULL

    stage <- "limma_minfi"
    if (!is.null(log_fn)) {
      log_fn(
        "run_limma input (model %s): M_minfi dims=%s×%s, beta_minfi dims=%s×%s, rownames(M)=%s, rownames(beta)=%s",
        params$model_id %||% "<unnamed>",
        nrow(M_minfi_corrected),
        ncol(M_minfi_corrected),
        nrow(beta_minfi_corrected),
        ncol(beta_minfi_corrected),
        if (is.null(rownames(M_minfi_corrected))) "NULL" else "set",
        if (is.null(rownames(beta_minfi_corrected))) "NULL" else "set"
      )
    }
    results_minfi <- run_limma(M_minfi_corrected, beta_minfi_corrected, metadata, design_with_sv)
    results_minfi <- limit_top_hits(results_minfi, opt$top_n_cpgs)

    stage <- "limma_sesame"
    if (sesame_available && !is.null(M_sesame_corrected) && !is.null(beta_sesame_corrected)) {
      if (!is.null(log_fn)) {
        log_fn(
          "run_limma sesame input (model %s): M_sesame dims=%s×%s, beta_sesame dims=%s×%s, rownames(M)=%s, rownames(beta)=%s",
          params$model_id %||% "<unnamed>",
          nrow(M_sesame_corrected),
          ncol(M_sesame_corrected),
          nrow(beta_sesame_corrected),
          ncol(beta_sesame_corrected),
          if (is.null(rownames(M_sesame_corrected))) "NULL" else "set",
          if (is.null(rownames(beta_sesame_corrected))) "NULL" else "set"
        )
      }
      results_sesame <- run_limma(M_sesame_corrected, beta_sesame_corrected, metadata, design_with_sv)
      results_sesame <- limit_top_hits(results_sesame, opt$top_n_cpgs)
    } else {
      results_sesame <- results_minfi[0]
    }

    stage <- "metric_collection"
    sig_minfi <- filter_significant(results_minfi, opt$fdr_threshold, opt$delta_beta_threshold)
    sig_sesame <- filter_significant(results_sesame, opt$fdr_threshold, opt$delta_beta_threshold)
    sig_intersection <- if (nrow(sig_minfi) > 0 && nrow(sig_sesame) > 0) length(intersect(sig_minfi$probe_id, sig_sesame$probe_id)) else 0L
    batch_median_p_minfi <- score_batch_effect(pca_minfi_post, candidate_batches)
    batch_median_p_sesame <- if (sesame_available && !is.null(M_sesame_corrected)) score_batch_effect(pca_sesame_post, candidate_batches) else NA_real_
    group_stats <- pca_minfi_post$assessments[variable == "dear_group"]
    group_median_p <- if (nrow(group_stats) > 0) median(group_stats$p_value, na.rm = TRUE) else NA_real_

    metrics <- list(
      numeric_covars = numeric_final,
      factor_covars = factor_final,
      use_sva = params$use_sva,
      use_combat = params$use_combat,
      batch_col = params$batch_col,
      batch_in_design = !params$use_combat && !is.null(params$batch_col),
      n_surrogates = if (!is.null(surrogate_vars)) ncol(surrogate_vars) else 0L,
      batch_median_p_minfi = batch_median_p_minfi,
      batch_median_p_sesame = batch_median_p_sesame,
      group_median_p = group_median_p,
      n_sig_minfi = nrow(sig_minfi),
      n_sig_sesame = nrow(sig_sesame),
      n_sig_intersection = sig_intersection
    )

    outputs <- NULL
    if (isTRUE(keep_outputs)) {
      outputs <- list(
        design_base = design_base,
        design_with_sv = design_with_sv,
        surrogate_vars = surrogate_vars,
        M_minfi = M_minfi_corrected,
        beta_minfi = beta_minfi_corrected,
        M_sesame = M_sesame_corrected,
        beta_sesame = beta_sesame_corrected,
        pca_minfi_post = pca_minfi_post,
        pca_sesame_post = pca_sesame_post,
        results_minfi = results_minfi,
        results_sesame = results_sesame,
        covariates_final = list(numeric = numeric_final, factor = factor_final),
        n_surrogates = metrics$n_surrogates
      )
    }

    list(status = "ok", metrics = metrics, outputs = outputs, message = NULL)
  }, error = function(e) {
    message <- conditionMessage(e)
    detail <- sprintf("%s (stage: %s)", message, stage)
    if (!is.null(log_fn)) {
      log_fn("Model %s failed: %s", params$model_id %||% "<unnamed>", detail)
      tb <- utils::capture.output(traceback())
      if (length(tb) > 0) {
        log_fn("Traceback for %s: %s", params$model_id %||% "<unnamed>", paste(tb, collapse = " | "))
      }
    }
    list(status = "error", metrics = NULL, outputs = NULL, message = detail)
  })
  result
}

optimize_design_matrix <- function(metadata, covariates, manual_covariates, protected_columns, candidate_batches, manual_batch_columns, M_minfi, M_sesame, sesame_available, opt, log_fn = log_message) {
  required_numeric <- unique(manual_covariates$numeric)
  required_factor <- unique(c(manual_covariates$factor, intersect(protected_columns, covariates$factor)))
  optional_numeric <- setdiff(covariates$numeric, required_numeric)
  optional_factor <- setdiff(covariates$factor, required_factor)
  optional_stats <- covariate_association_stats(metadata, "dear_group", numeric_covars = optional_numeric, factor_covars = optional_factor)
  stats_numeric <- optional_stats[type == "numeric"]
  stats_factor <- optional_stats[type == "factor"]
  covariate_sets <- generate_covariate_sets(required_numeric, required_factor, stats_numeric, stats_factor)
  batch_options <- build_batch_options(metadata, "dear_group", manual_batch_columns, candidate_batches)
  if (!is.null(log_fn)) {
    log_fn("Evaluating %s covariate sets across %s batch options (including none).", length(covariate_sets), length(batch_options))
  }
  data_payload <- list(
    metadata = metadata,
    M_minfi = M_minfi,
    M_sesame = M_sesame,
    sesame_available = sesame_available
  )
  configs <- list()
  for (cov_idx in seq_along(covariate_sets)) {
    cov_set <- covariate_sets[[cov_idx]]
    for (batch_val in batch_options) {
      batch_col <- if (is.na(batch_val)) NULL else batch_val
      combat_opts <- if (is.null(batch_col)) c(FALSE) else c(FALSE, TRUE)
      for (use_combat in combat_opts) {
        for (use_sva in c(FALSE, TRUE)) {
          model_id <- format_model_id(cov_set$id, batch_col, use_combat, use_sva)
          configs[[length(configs) + 1L]] <- list(
            model_id = model_id,
            numeric_covars = cov_set$numeric,
            factor_covars = cov_set$factor,
            use_combat = use_combat,
            use_sva = use_sva,
            batch_col = batch_col
          )
        }
      }
    }
  }
  if (length(configs) == 0) {
    stop("No model configurations generated for optimisation")
  }
  if (!is.null(log_fn)) {
    log_fn("Generated %s model configurations for optimisation.", length(configs))
  }
  evaluation <- vector("list", length(configs))
  for (i in seq_along(configs)) {
    params <- configs[[i]]
    result <- execute_model_configuration(params, data_payload, opt, candidate_batches, log_fn = log_fn, keep_outputs = FALSE)
    metrics <- result$metrics
    evaluation[[i]] <- list(
      model_id = params$model_id,
      status = result$status,
      message = result$message %||% "",
      batch_col = params$batch_col %||% NA_character_,
      use_combat = params$use_combat,
      use_sva = params$use_sva,
      batch_in_design = if (!is.null(metrics)) metrics$batch_in_design else NA,
      n_surrogates = if (!is.null(metrics)) metrics$n_surrogates else NA_integer_,
      n_covariates_numeric = if (!is.null(metrics)) length(metrics$numeric_covars) else NA_integer_,
      n_covariates_factor = if (!is.null(metrics)) length(metrics$factor_covars) else NA_integer_,
      batch_median_p_minfi = if (!is.null(metrics)) metrics$batch_median_p_minfi else NA_real_,
      batch_median_p_sesame = if (!is.null(metrics)) metrics$batch_median_p_sesame else NA_real_,
      group_median_p = if (!is.null(metrics)) metrics$group_median_p else NA_real_,
      n_sig_minfi = if (!is.null(metrics)) metrics$n_sig_minfi else NA_integer_,
      n_sig_sesame = if (!is.null(metrics)) metrics$n_sig_sesame else NA_integer_,
      n_sig_intersection = if (!is.null(metrics)) metrics$n_sig_intersection else NA_integer_
    )
  }
  evaluation_dt <- rbindlist(evaluation, fill = TRUE)
  evaluation_dt[, n_covariates_total := n_covariates_numeric + n_covariates_factor]
  valid_models <- evaluation_dt[status == "ok"]
  if (nrow(valid_models) == 0) {
    if (!is.null(log_fn)) {
      log_fn("All model configurations failed during optimisation; falling back will be required.")
    }
    return(list(
      status = "failed",
      best = NULL,
      evaluated = evaluation_dt,
      covariate_sets = covariate_sets,
      covariate_stats = list(optional_numeric = stats_numeric, optional_factor = stats_factor),
      batch_options = batch_options,
      error = "All model configurations failed"
    ))
  }
  valid_models[, batch_ok := is.na(batch_median_p_minfi) | batch_median_p_minfi >= 0.2]
  valid_models[, batch_score := ifelse(is.na(batch_median_p_minfi), 0, batch_median_p_minfi)]
  setorder(
    valid_models,
    -batch_ok,
    -n_sig_intersection,
    -n_sig_minfi,
    -n_sig_sesame,
    -batch_score,
    n_covariates_total,
    n_surrogates
  )
  best_row <- valid_models[1]
  idx_map <- vapply(configs, function(x) x$model_id, character(1))
  best_params <- configs[[match(best_row$model_id, idx_map)]]
  best_result <- execute_model_configuration(best_params, data_payload, opt, candidate_batches, log_fn = log_fn, keep_outputs = TRUE)
  if (!identical(best_result$status, "ok")) {
    stop("Failed to execute best model configuration: ", best_result$message %||% "unknown error")
  }
  numeric_desc <- if (length(best_result$metrics$numeric_covars) == 0) "none" else paste(best_result$metrics$numeric_covars, collapse = ", ")
  factor_desc <- if (length(best_result$metrics$factor_covars) == 0) "none" else paste(best_result$metrics$factor_covars, collapse = ", ")
  if (!is.null(log_fn)) {
    log_fn(
      "Selected design %s (batch=%s, combat=%s, sva=%s, numeric=%s, factor=%s, sig_minfi=%s, sig_sesame=%s, sig_intersection=%s)",
      best_row$model_id,
      best_params$batch_col %||% "none",
      best_params$use_combat,
      best_params$use_sva,
      numeric_desc,
      factor_desc,
      best_result$metrics$n_sig_minfi,
      best_result$metrics$n_sig_sesame,
      best_result$metrics$n_sig_intersection
    )
  }
  list(
    status = "ok",
    best = list(
      params = best_params,
      metrics = best_result$metrics,
      outputs = best_result$outputs
    ),
    evaluated = evaluation_dt,
    covariate_sets = covariate_sets,
    covariate_stats = list(optional_numeric = stats_numeric, optional_factor = stats_factor),
    batch_options = batch_options
  )
}

fallback_design_selection <- function(metadata, covariates, manual_covariates, protected_columns, candidate_batches, M_minfi, M_sesame, sesame_available, opt, log_fn = log_message) {
  if (!is.null(log_fn)) {
    log_fn("Attempting fallback design selection with simplified configurations.")
  }
  data_payload <- list(
    metadata = metadata,
    M_minfi = M_minfi,
    M_sesame = M_sesame,
    sesame_available = sesame_available
  )
  protected_factors <- intersect(protected_columns, covariates$factor)
  fallback_configs <- list(
    list(
      model_id = "fallback_full",
      numeric_covars = unique(covariates$numeric),
      factor_covars = unique(covariates$factor),
      use_combat = FALSE,
      use_sva = FALSE,
      batch_col = NULL
    ),
    list(
      model_id = "fallback_full_sva",
      numeric_covars = unique(covariates$numeric),
      factor_covars = unique(covariates$factor),
      use_combat = FALSE,
      use_sva = TRUE,
      batch_col = NULL
    ),
    list(
      model_id = "fallback_manual",
      numeric_covars = unique(manual_covariates$numeric),
      factor_covars = unique(c(manual_covariates$factor, protected_factors)),
      use_combat = FALSE,
      use_sva = FALSE,
      batch_col = NULL
    ),
    list(
      model_id = "fallback_minimal",
      numeric_covars = character(),
      factor_covars = unique(protected_factors),
      use_combat = FALSE,
      use_sva = FALSE,
      batch_col = NULL
    )
  )
  evaluations <- list()
  for (cfg in fallback_configs) {
    res <- execute_model_configuration(cfg, data_payload, opt, candidate_batches, log_fn = log_fn, keep_outputs = TRUE)
    evaluations[[length(evaluations) + 1L]] <- list(
      model_id = cfg$model_id,
      status = res$status,
      message = res$message %||% "",
      batch_col = cfg$batch_col %||% NA_character_,
      use_combat = cfg$use_combat,
      use_sva = cfg$use_sva
    )
    if (identical(res$status, "ok")) {
      if (!is.null(log_fn)) {
        log_fn("Fallback design %s succeeded.", cfg$model_id)
      }
      evaluation_dt <- rbindlist(evaluations, fill = TRUE)
      return(list(
        status = "ok",
        best = list(params = cfg, metrics = res$metrics, outputs = res$outputs),
        evaluated = evaluation_dt,
        covariate_sets = list(),
        covariate_stats = list(optional_numeric = data.table(), optional_factor = data.table()),
        batch_options = c(NA_character_),
        note = "fallback"
      ))
    }
  }
  evaluation_dt <- if (length(evaluations) > 0) rbindlist(evaluations, fill = TRUE) else data.table()
  stop("Fallback design selection failed; no viable configuration executed successfully.")
}

run_limma <- function(M_matrix, beta_matrix, metadata, design) {
  model_id <- getOption("dearmeta.current_model_id", "<unknown>")
  runtime_dir <- getOption("dearmeta.runtime_dir")
  tryCatch({
    M_matrix <- as.matrix(M_matrix)
    beta_matrix <- as.matrix(beta_matrix)
    log_message("run_limma[%s]: initial dims M=%s×%s, beta=%s×%s", model_id, nrow(M_matrix), ncol(M_matrix), nrow(beta_matrix), ncol(beta_matrix))
    if (!identical(rownames(M_matrix), rownames(beta_matrix))) {
      shared <- intersect(rownames(M_matrix), rownames(beta_matrix))
      if (length(shared) == 0) {
        stop("No shared probes between M matrix and beta matrix")
      }
      missing_m <- setdiff(rownames(M_matrix), shared)
      missing_b <- setdiff(rownames(beta_matrix), shared)
      if (length(missing_m) > 0) {
        log_message("run_limma: dropping %s probes absent from beta matrix prior to alignment.", length(missing_m))
      }
      if (length(missing_b) > 0) {
        log_message("run_limma: dropping %s probes absent from M matrix prior to alignment.", length(missing_b))
      }
      M_matrix <- M_matrix[shared, , drop = FALSE]
      beta_matrix <- beta_matrix[shared, , drop = FALSE]
    }
    keep_m <- rowSums(is.finite(M_matrix)) == ncol(M_matrix)
    keep_beta <- rowSums(is.finite(beta_matrix)) == ncol(beta_matrix)
    keep <- keep_m & keep_beta
    if (!all(keep)) {
      removed <- sum(!keep)
      log_message("run_limma[%s]: dropping %s probes with missing values prior to limma.", model_id, removed)
      M_matrix <- M_matrix[keep, , drop = FALSE]
      beta_matrix <- beta_matrix[keep, , drop = FALSE]
    }
    if (nrow(M_matrix) == 0) {
      stop("No probes remain after filtering missing values; aborting limma step.")
    }
    var_m <- apply(M_matrix, 1, var, na.rm = TRUE)
    keep_var <- is.finite(var_m) & (var_m > 0)
    if (!all(keep_var)) {
      removed_var <- sum(!keep_var)
      log_message("run_limma[%s]: dropping %s probes with zero variance prior to limma.", model_id, removed_var)
      M_matrix <- M_matrix[keep_var, , drop = FALSE]
      beta_matrix <- beta_matrix[keep_var, , drop = FALSE]
    }
    log_message("run_limma[%s]: proceeding with %s probes", model_id, nrow(M_matrix))
    fit <- lmFit(M_matrix, design)
    fit <- eBayes(fit)
    group_coef <- grep("^group", colnames(design), value = TRUE)
    if (length(group_coef) == 0) {
      stop("Design matrix lacks a dear_group coefficient")
    }
    group_coef <- group_coef[1]
    log_message("run_limma[%s]: using coefficient %s for group comparison", model_id, group_coef)
    top <- topTable(fit, coef = group_coef, number = nrow(M_matrix), sort.by = "P")
    top <- data.table(probe_id = rownames(top), top)
    available <- top$probe_id[top$probe_id %in% rownames(beta_matrix)]
    if (length(available) == 0) {
      stop("No overlapping probes between limma results and beta matrix")
    }
    if (length(available) < nrow(top)) {
      dropped <- setdiff(top$probe_id, available)
      log_message("run_limma[%s]: omitting %s probes absent from beta matrix during aggregation.", model_id, length(dropped))
      top <- top[probe_id %in% available]
    }
    beta_df <- tryCatch(
      beta_matrix[top$probe_id, , drop = FALSE],
      error = function(e) {
        log_message(
          "run_limma[%s]: beta subset failure with %s probes; head ids: %s",
          model_id,
          length(top$probe_id),
          paste(head(top$probe_id, 5), collapse = ",")
        )
        log_message(
          "run_limma[%s]: beta rowname head: %s",
          model_id,
          paste(head(rownames(beta_matrix), 5), collapse = ",")
        )
        stop(e)
      }
    )
    log_message("run_limma[%s]: beta subset succeeded", model_id)
    log_message("run_limma[%s]: beta_df dims=%s×%s", model_id, nrow(beta_df), ncol(beta_df))
    group_levels <- levels(factor(metadata$dear_group))
    if (length(group_levels) != 2) {
      stop("Limma model currently supports two-group comparison.")
    }
    beta_group1 <- tryCatch(rowMeans(beta_df[, metadata$dear_group == group_levels[1], drop = FALSE], na.rm = TRUE), error = function(e) {
      log_message("run_limma[%s]: rowMeans group1 failed", model_id)
      stop(e)
    })
    beta_group2 <- tryCatch(rowMeans(beta_df[, metadata$dear_group == group_levels[2], drop = FALSE], na.rm = TRUE), error = function(e) {
      log_message("run_limma[%s]: rowMeans group2 failed", model_id)
      stop(e)
    })
    top[, delta_beta := beta_group2 - beta_group1]
    top[, direction := ifelse(delta_beta > 0, "hypermethylated", "hypomethylated")]
    top
  }, error = function(e) {
    if (!is.null(runtime_dir) && dir.exists(runtime_dir)) {
      safe_model <- gsub("[^A-Za-z0-9_]+", "_", model_id)
      debug_path <- file.path(runtime_dir, sprintf("run_limma_debug_%s_%s.rds", safe_model, format(Sys.time(), "%Y%m%d%H%M%S")))
      debug_payload <- list(
        error = conditionMessage(e),
        model_id = model_id,
        M_matrix = M_matrix,
        beta_matrix = beta_matrix,
        metadata = metadata,
        design = design,
        available_top = if (exists("top")) top else NULL
      )
      tryCatch(saveRDS(debug_payload, debug_path), error = function(save_err) {
        log_message("run_limma[%s]: failed to write debug snapshot: %s", model_id, conditionMessage(save_err))
      })
      log_message("run_limma[%s]: saved debug snapshot to %s", model_id, debug_path)
    }
    stop(e)
  })
}

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

filter_significant <- function(res, fdr, delta_thresh) {
  res[adj.P.Val <= fdr & abs(delta_beta) >= delta_thresh]
}

merge_results <- function(minfi, sesame) {
  minfi[, pipeline := "minfi"]
  sesame[, pipeline := "sesame"]
  union <- rbindlist(list(minfi, sesame), fill = TRUE, use.names = TRUE)
  union[, shared_significant := probe_id %in% intersect(minfi$probe_id, sesame$probe_id)]
  intersection <- minfi[probe_id %in% sesame$probe_id]
  list(union = union, intersection = intersection)
}

top_cpg_table_name <- function(dt, n = 100) {
  head(dt[order(P.Value)], n = n)$probe_id
}

write_tsv_gz <- function(data, path) {
  dt <- as.data.table(data)
  if (nrow(dt) == 0) {
    con <- gzfile(path, "w")
    on.exit(close(con), add = TRUE)
    header <- paste(names(dt), collapse = "\t")
    writeLines(header, con = con)
  } else {
    fwrite(dt, file = path, sep = "\t", compress = "gzip")
  }
}

save_matrix <- function(mat, path, compress = TRUE) {
  dt <- as.data.table(mat, keep.rownames = "probe_id")
  if (compress) {
    write_tsv_gz(dt, path)
  } else {
    fwrite(dt, file = path, sep = "\t")
  }
}

write_plot <- function(plot_obj, filename, width = 7, height = 5) {
  tryCatch(
    {
      build <- ggplot_build(plot_obj)
      layer_rows <- vapply(build$data, function(layer) nrow(layer), integer(1))
      total_rows <- sum(layer_rows)
      log_message("Rendering plot %s (%s data rows)", filename, total_rows)
      ggsave(filename = paste0(filename, ".png"), plot = plot_obj, width = width, height = height, dpi = 300)
      ggsave(filename = paste0(filename, ".pdf"), plot = plot_obj, width = width, height = height)
      invisible(TRUE)
    },
    error = function(e) {
      trace <- paste(capture.output(traceback()), collapse = " | ")
      log_message(sprintf("Skipping plot %s due to: %s | Trace: %s", filename, conditionMessage(e), trace))
      invisible(FALSE)
    }
  )
}

resolve_widget_title <- function(widget) {
  if (inherits(widget, "plotly")) {
    layout_title <- widget$x$layout$title
    if (is.list(layout_title)) {
      return(layout_title$text %||% layout_title)
    }
    if (is.character(layout_title) && length(layout_title) == 1) {
      return(layout_title)
    }
  }
  if (!is.null(widget$attrs) && is.list(widget$attrs) && !is.null(widget$attrs$title)) {
    title <- widget$attrs$title
    if (is.character(title) && length(title) == 1) {
      return(title)
    }
  }
  NULL
}

format_identifier_title <- function(identifier) {
  if (is.null(identifier) || !nzchar(identifier)) {
    return("DearMeta Interactive")
  }
  cleaned <- gsub("[^[:alnum:]]+", " ", identifier)
  parts <- unlist(strsplit(cleaned, "\\s+"))
  parts <- parts[parts != ""]
  if (length(parts) == 0) {
    return("DearMeta Interactive")
  }
  paste(
    vapply(
      parts,
      function(word) {
        first <- substr(word, 1, 1)
        if (!nzchar(first)) {
          return("")
        }
        rest <- ""
        if (nchar(word) > 1) {
          rest <- substr(word, 2, nchar(word))
        }
        paste0(toupper(first), tolower(rest))
      },
      character(1),
      USE.NAMES = FALSE
    ),
    collapse = " "
  )
}

default_widget_subtitle <- function(widget) {
  if (inherits(widget, "plotly")) {
    return("Drag to zoom, hover to inspect values, and use the toolbar to export high-resolution images.")
  }
  if ("datatables" %in% class(widget)) {
    return("Use the search box, column sorting, and export buttons to explore the data interactively.")
  }
  NULL
}

default_widget_description <- function(widget) {
  if (inherits(widget, "plotly")) {
    return("Tip: double-click anywhere on the chart to reset the view.")
  }
  if ("datatables" %in% class(widget)) {
    return("Tip: combine header filters with the export buttons to take data offline.")
  }
  NULL
}

interactive_output_metadata <- function(key) {
  switch(
    key,
    pca_pre = list(
      title = "PCA (Pre-correction)",
      subtitle = "Principal component analysis before correction. Hover to inspect sample IDs.",
      description = "Use the legend to toggle sample groups and assess clustering before batch correction.",
      category = "Quality Control",
      icon = "QC"
    ),
    pca_post = list(
      title = "PCA (Post-correction)",
      subtitle = "Principal component analysis after correction. Toggle legend entries to focus on a group.",
      description = "Compare post-correction group separation and confirm batch mitigation.",
      category = "Quality Control",
      icon = "QC"
    ),
    volcano_minfi = list(
      title = "Volcano · minfi",
      subtitle = "-log10(p-value) versus delta beta. Firebrick markers pass the FDR threshold.",
      description = "Hover for CpG IDs, delta beta, and p-values derived from the minfi pipeline.",
      category = "Differential Analysis",
      icon = "DA"
    ),
    volcano_sesame = list(
      title = "Volcano · sesame",
      subtitle = "-log10(p-value) versus delta beta for sesame results.",
      description = "Hold shift to zoom along an axis and use the toolbar to export high-resolution images.",
      category = "Differential Analysis",
      icon = "DA"
    ),
    volcano_intersection = list(
      title = "Volcano · intersection",
      subtitle = "Shared signals across pipelines visualised in a volcano plot.",
      description = "Inspect consistent CpGs that meet the filtering thresholds in both pipelines.",
      category = "Differential Analysis",
      icon = "DA"
    ),
    manhattan_minfi = list(
      title = "Manhattan · minfi",
      subtitle = "Genome-wide distribution of -log10(p-value) signals for minfi.",
      description = "Drag horizontally to zoom into genomic regions and inspect CpG annotations.",
      category = "Genomic Context",
      icon = "GX"
    ),
    manhattan_sesame = list(
      title = "Manhattan · sesame",
      subtitle = "Sesame-based Manhattan plot with interactive zoom and hover.",
      description = "Toggle legend entries to highlight one group or significance band at a time.",
      category = "Genomic Context",
      icon = "GX"
    ),
    manhattan_intersection = list(
      title = "Manhattan · intersection",
      subtitle = "Intersection hits highlighted across the genome.",
      description = "Spot reproducible genomic hotspots with consistent differential methylation.",
      category = "Genomic Context",
      icon = "GX"
    ),
    table_minfi = list(
      title = "Annotated CpGs · minfi",
      subtitle = "Search, filter, sort, and export the table using the controls above.",
      description = "Interactive table containing annotated minfi CpGs retained for reporting.",
      category = "Annotated Tables",
      icon = "TB"
    ),
    table_sesame = list(
      title = "Annotated CpGs · sesame",
      subtitle = "Search, filter, sort, and export the table using the controls above.",
      description = "Interactive table containing annotated sesame CpGs retained for reporting.",
      category = "Annotated Tables",
      icon = "TB"
    ),
    table_intersection = list(
      title = "Annotated CpGs · intersection",
      subtitle = "Search, filter, sort, and export the table using the controls above.",
      description = "Annotated CpGs that are shared between pipelines after filtering.",
      category = "Annotated Tables",
      icon = "TB"
    ),
    table_intersection_dual = list(
      title = "Shared CpGs · dual pipeline view",
      subtitle = "Compare minfi and sesame effect sizes, FDR, and direction for overlapping CpGs.",
      description = "Inspect concordance and discrepancies across pipelines; export for downstream validation.",
      category = "Annotated Tables",
      icon = "TB"
    ),
    model_selection_table = list(
      title = "Design Optimiser · Summary",
      subtitle = "Ranked view of candidate batch/covariate models and evaluation metrics.",
      description = "Review how batch choices, covariates, and SVA/ComBat decisions were scored during automatic selection.",
      category = "Batch & Covariates",
      icon = "BC"
    ),
    default = list(
      title = format_identifier_title(key),
      subtitle = "Interactive visual overview.",
      description = "Explore this interactive output.",
      category = "Additional Outputs",
      icon = "**"
    )
  )
}

escape_regex <- function(text) {
  gsub("([][{}()+*^$|\\\\.?-])", "\\\\\\1", text)
}

relative_to_root <- function(target, root) {
  root_abs <- normalizePath(root, winslash = "/", mustWork = FALSE)
  target_abs <- normalizePath(target, winslash = "/", mustWork = FALSE)
  if (!nzchar(root_abs) || !nzchar(target_abs)) {
    return(basename(target))
  }
  pattern <- paste0("^", escape_regex(root_abs), "/?")
  rel <- sub(pattern, "", target_abs)
  rel <- gsub("^/+", "", rel)
  if (!nzchar(rel)) {
    rel <- basename(target_abs)
  }
  rel
}

format_number <- function(value, digits = 0) {
  if (is.null(value) || length(value) == 0 || all(is.na(value))) {
    return("0")
  }
  formatted <- format(round(as.numeric(value), digits), big.mark = ",", trim = TRUE, scientific = FALSE, nsmall = digits)
  if (length(formatted) == 0 || is.na(formatted)) {
    return("0")
  }
  formatted
}

format_percent <- function(value, digits = 1) {
  if (is.null(value) || length(value) == 0 || all(is.na(value))) {
    return("0%")
  }
  paste0(sprintf(paste0("%.", digits, "f"), as.numeric(value) * 100), "%")
}

format_bytes <- function(bytes) {
  if (is.null(bytes) || length(bytes) == 0 || all(is.na(bytes))) {
    return("0 B")
  }
  units <- c("B", "KB", "MB", "GB", "TB")
  idx <- 1
  value <- as.numeric(bytes)
  while (value >= 1024 && idx < length(units)) {
    value <- value / 1024
    idx <- idx + 1
  }
  sprintf("%.1f %s", value, units[idx])
}

style_plotly_widget <- function(widget, export_name) {
  if (!is.null(widget$x$data) && length(widget$x$data) > 0) {
    widget$x$data <- lapply(widget$x$data, function(trace) {
      mode <- trace$mode %||% ""
      type <- trace$type %||% ""
      if (grepl("markers", mode, fixed = TRUE) || identical(type, "scattergl")) {
        marker <- trace$marker %||% list()
        if (!is.null(marker$size)) {
          if (length(marker$size) <= 1) {
            marker$size <- max(as.numeric(marker$size %||% 0) * 2, 20)
          } else {
            marker$size <- pmax(as.numeric(marker$size) * 2, 20)
          }
        } else {
          marker$size <- 20
        }
        marker$opacity <- marker$opacity %||% 0.82
        line <- marker$line %||% list()
        line$width <- line$width %||% 0
        line$color <- line$color %||% "rgba(15,23,42,0.35)"
        marker$line <- line
        trace$marker <- marker
      }
      trace
    })
  }
  existing_layout <- widget$x$layout %||% list()
  existing_title <- existing_layout$title %||% list()
  if (is.list(existing_title)) {
    existing_title$font <- modifyList(
      existing_title$font %||% list(),
      list(family = "Inter, 'Segoe UI', 'Helvetica Neue', sans-serif", color = "#F8FAFC", size = 20)
    )
  } else if (is.character(existing_title) && length(existing_title) == 1) {
    existing_title <- list(
      text = existing_title,
      font = list(family = "Inter, 'Segoe UI', 'Helvetica Neue', sans-serif", color = "#F8FAFC", size = 20)
    )
  }

  existing_xaxis <- existing_layout$xaxis %||% list()
  existing_x_title <- NULL
  if (!is.null(existing_xaxis$title)) {
    if (is.list(existing_xaxis$title)) {
      existing_x_title <- existing_xaxis$title$text %||% existing_xaxis$title
    } else if (is.character(existing_xaxis$title)) {
      existing_x_title <- existing_xaxis$title
    }
  }
  xaxis_defaults <- list(
    zeroline = FALSE,
    showline = TRUE,
    linecolor = "rgba(148, 163, 184, 0.4)",
    tickfont = list(color = "#CBD5E1", size = 13),
    gridcolor = "rgba(148, 163, 184, 0.18)",
    gridwidth = 0.7,
    ticks = "outside",
    tickcolor = "rgba(148, 163, 184, 0.45)",
    title = list(font = list(color = "#F8FAFC", size = 17))
  )
  if (!is.null(existing_x_title)) {
    xaxis_defaults$title$text <- existing_x_title
  }
  new_xaxis <- modifyList(existing_xaxis, xaxis_defaults)

  existing_yaxis <- existing_layout$yaxis %||% list()
  existing_y_title <- NULL
  if (!is.null(existing_yaxis$title)) {
    if (is.list(existing_yaxis$title)) {
      existing_y_title <- existing_yaxis$title$text %||% existing_yaxis$title
    } else if (is.character(existing_yaxis$title)) {
      existing_y_title <- existing_yaxis$title
    }
  }
  yaxis_defaults <- list(
    zeroline = FALSE,
    showline = TRUE,
    linecolor = "rgba(148, 163, 184, 0.4)",
    tickfont = list(color = "#CBD5E1", size = 13),
    gridcolor = "rgba(148, 163, 184, 0.18)",
    gridwidth = 0.7,
    ticks = "outside",
    tickcolor = "rgba(148, 163, 184, 0.45)",
    title = list(font = list(color = "#F8FAFC", size = 17))
  )
  if (!is.null(existing_y_title)) {
    yaxis_defaults$title$text <- existing_y_title
  }
  new_yaxis <- modifyList(existing_yaxis, yaxis_defaults)

  legend_defaults <- list(
    orientation = "h",
    bgcolor = "rgba(15, 23, 42, 0.6)",
    bordercolor = "rgba(148, 163, 184, 0.25)",
    borderwidth = 1,
    xanchor = "center",
    x = 0.5,
    y = -0.18,
    font = list(color = "#E2E8F0", size = 13)
  )
  new_legend <- modifyList(existing_layout$legend %||% list(), legend_defaults)

  widget %>%
    plotly::layout(
      margin = list(l = 70, r = 40, t = 65, b = 70),
      font = list(
        family = "Inter, 'Segoe UI', 'Helvetica Neue', sans-serif",
        color = "#E2E8F0",
        size = 14
      ),
      paper_bgcolor = "rgba(15, 23, 42, 0)",
      plot_bgcolor = "rgba(15, 23, 42, 0.92)",
      hovermode = "closest",
      hoverlabel = list(
        bgcolor = "#111827",
        bordercolor = "#1f2937",
        font = list(color = "#F8FAFC", size = 13)
      ),
      xaxis = new_xaxis,
      yaxis = new_yaxis,
      legend = new_legend,
      title = existing_title
    ) %>%
    plotly::config(
      displaylogo = FALSE,
      responsive = TRUE,
      modeBarButtonsToRemove = c("lasso2d", "select2d", "autoScale2d"),
      toImageButtonOptions = list(
        format = "png",
        filename = export_name,
        width = 1400,
        height = 900,
        scale = 2
      )
    )
}

style_datatable_widget <- function(widget) {
  buttons_enabled <- "Buttons" %in% (widget$x$extensions %||% character())
  default_options <- list(
    pageLength = 25,
    lengthMenu = list(c(10, 25, 50, 100, -1), c("10", "25", "50", "100", "All")),
    scrollX = TRUE,
    autoWidth = TRUE,
    deferRender = TRUE,
    dom = if (buttons_enabled) "<'dm-toolbar'Bfr>t<'dm-footer'ip>" else "<'dm-toolbar'f>t<'dm-footer'ip>",
    language = list(
      search = "Search:",
      lengthMenu = "Show _MENU_ rows per page",
      info = "Showing _START_ to _END_ of _TOTAL_ rows",
      infoEmpty = "No rows to display",
      paginate = list(previous = "Prev", `next` = "Next")
    )
  )
  if (buttons_enabled) {
    default_options$buttons <- list(
      list(extend = "copy", className = "btn-copy", text = "Copy"),
      list(extend = "csv", className = "btn-csv", text = "CSV"),
      list(extend = "excel", className = "btn-excel", text = "Excel")
    )
  } else {
    widget$x$options$buttons <- NULL
  }
  widget$x$options <- modifyList(default_options, widget$x$options %||% list())
  if (buttons_enabled) {
    widget$x$extensions <- unique(c(widget$x$extensions %||% character(), "Buttons"))
  }
  existing_class <- widget$x$class %||% ""
  if (!grepl("stripe", existing_class, fixed = TRUE)) {
    widget$x$class <- trimws(paste(existing_class, "stripe hover row-border compact"))
  }
  widget
}

style_widget <- function(widget, export_name) {
  if (inherits(widget, "plotly")) {
    widget <- style_plotly_widget(widget, export_name)
  }
  if ("datatables" %in% class(widget)) {
    widget <- style_datatable_widget(widget)
  }
  widget
}

build_interactive_page <- function(widget, title, subtitle, description) {
  css <- "
:root { color-scheme: dark; }
* { box-sizing: border-box; }
body {
  margin: 0;
  padding: 0;
  font-family: 'Inter', 'Segoe UI', 'Helvetica Neue', sans-serif;
  background: radial-gradient(120% 160% at 50% 0%, rgba(59,130,246,0.25) 0%, rgba(15,23,42,1) 55%);
  color: #E2E8F0;
}
.dm-shell {
  min-height: 100vh;
  display: flex;
  flex-direction: column;
  align-items: center;
  gap: 40px;
  padding: 64px 28px 72px;
}
.dm-header {
  width: min(1360px, 100%);
  text-align: left;
}
.dm-breadcrumb {
  text-transform: uppercase;
  letter-spacing: 0.18em;
  font-size: 0.75rem;
  color: rgba(148, 163, 184, 0.65);
  margin: 0 0 12px;
  font-weight: 600;
}
.dm-title {
  margin: 0;
  font-size: clamp(2.1rem, 4vw, 2.6rem);
  font-weight: 700;
  color: #F8FAFC;
}
.dm-subtitle {
  margin: 10px 0 0;
  font-size: 1rem;
  color: rgba(226, 232, 240, 0.82);
  max-width: 880px;
  line-height: 1.55;
}
.dm-main {
  width: min(1360px, 100%);
}
.dm-card {
  background: rgba(15, 23, 42, 0.82);
  border-radius: 28px;
  padding: clamp(24px, 4vw, 42px);
  box-shadow: 0 38px 78px rgba(15, 23, 42, 0.55);
  border: 1px solid rgba(148, 163, 184, 0.2);
  backdrop-filter: blur(18px);
}
.dm-widget {
  width: 100%;
  position: relative;
}
.dm-widget .plotly.html-widget {
  width: 100% !important;
  min-height: clamp(420px, 60vh, 640px);
}
.dm-widget .datatables.html-widget {
  width: 100% !important;
}
.dm-caption {
  margin-top: 18px;
  font-size: 0.95rem;
  color: rgba(226, 232, 240, 0.75);
}
.dm-footer {
  width: min(1360px, 100%);
  text-align: center;
  font-size: 0.85rem;
  color: rgba(148, 163, 184, 0.7);
  padding-bottom: 16px;
}
.modebar {
  background: rgba(15, 23, 42, 0.86) !important;
  border-radius: 16px !important;
  padding: 6px 8px !important;
  box-shadow: 0 18px 40px rgba(15, 23, 42, 0.55);
}
.modebar-group { border-right: none !important; }
.modebar-btn {
  color: #E2E8F0 !important;
  border-radius: 8px !important;
}
.modebar-btn:hover {
  color: #38BDF8 !important;
  background-color: rgba(56, 189, 248, 0.18) !important;
}
.dataTables_wrapper .dt-buttons .btn {
  background: linear-gradient(135deg, #38BDF8, #6366F1);
  border: none;
  border-radius: 12px;
  padding: 8px 18px;
  color: #0F172A;
  font-weight: 600;
  box-shadow: 0 12px 25px rgba(56, 189, 248, 0.35);
  margin-right: 10px;
}
.dataTables_wrapper .dt-buttons .btn:hover {
  filter: brightness(1.05);
  box-shadow: 0 16px 32px rgba(99, 102, 241, 0.35);
}
.dataTables_wrapper .dataTables_filter input {
  background: rgba(15, 23, 42, 0.6);
  border: 1px solid rgba(148, 163, 184, 0.35);
  border-radius: 12px;
  padding: 8px 12px;
  color: #F8FAFC;
}
.dataTables_wrapper .dataTables_length select {
  background: rgba(15, 23, 42, 0.6);
  border: 1px solid rgba(148, 163, 184, 0.35);
  border-radius: 12px;
  padding: 6px 8px;
  color: #F8FAFC;
}
table.dataTable tbody tr:hover {
  background-color: rgba(59, 130, 246, 0.14) !important;
}
table.dataTable thead th {
  border-bottom: 1px solid rgba(148, 163, 184, 0.4);
}
@media (max-width: 720px) {
  .dm-shell {
    padding: 44px 16px 54px;
  }
  .dm-card {
    padding: 22px 18px 32px;
  }
  .modebar {
    right: 16px !important;
    top: 16px !important;
  }
  .dataTables_wrapper .dt-buttons {
    display: flex;
    flex-wrap: wrap;
    gap: 8px;
  }
}
"
  page_title <- title %||% "DearMeta Interactive"
  subtitle_tag <- if (!is.null(subtitle)) htmltools::tags$p(class = "dm-subtitle", subtitle) else NULL
  caption_tag <- if (!is.null(description)) htmltools::tags$p(class = "dm-caption", description) else NULL
  htmltools::tags$html(
    lang = "en",
    htmltools::tags$head(
      htmltools::tags$meta(charset = "utf-8"),
      htmltools::tags$meta(name = "viewport", content = "width=device-width, initial-scale=1"),
      htmltools::tags$title(page_title),
      htmltools::tags$style(htmltools::HTML(css))
    ),
    htmltools::tags$body(
      htmltools::tags$div(
        class = "dm-shell",
        htmltools::tags$header(
          class = "dm-header",
          htmltools::tags$div(class = "dm-breadcrumb", "DearMeta · Interactive report"),
          htmltools::tags$h1(class = "dm-title", page_title),
          subtitle_tag
        ),
        htmltools::tags$main(
          class = "dm-main",
          htmltools::tags$div(
            class = "dm-card",
            htmltools::tags$div(class = "dm-widget", widget),
            caption_tag
          )
        ),
        htmltools::tags$footer(
          class = "dm-footer",
          "Generated with DearMeta"
        )
      )
    )
  )
}

resolve_widget_metadata <- function(widget, html_path, title = NULL, subtitle = NULL, description = NULL) {
  file_stub <- tools::file_path_sans_ext(basename(html_path))
  list(
    title = title %||% resolve_widget_title(widget) %||% format_identifier_title(file_stub),
    subtitle = subtitle %||% default_widget_subtitle(widget),
    description = description %||% default_widget_description(widget),
    export_name = gsub("[^A-Za-z0-9_-]+", "_", file_stub)
  )
}

save_styled_widget <- function(widget, html_path, title = NULL, subtitle = NULL, description = NULL) {
  metadata <- resolve_widget_metadata(widget, html_path, title, subtitle, description)
  styled_widget <- style_widget(widget, metadata$export_name)
  page <- build_interactive_page(styled_widget, metadata$title, metadata$subtitle, metadata$description)
  htmlwidgets::saveWidget(
    htmltools::browsable(page),
    file = html_path,
    selfcontained = FALSE,
    libdir = paste0(basename(html_path), "_files"),
    title = metadata$title
  )
  html_path
}

create_interactive <- function(plot_obj, filename, title = NULL, subtitle = NULL, description = NULL) {
  html_path <- paste0(filename, ".html")
  raw_widget <- plot_obj
  simple_save <- function(widget_to_save, custom_title = NULL) {
    htmlwidgets::saveWidget(
      widget_to_save,
      file = html_path,
      selfcontained = FALSE,
      libdir = paste0(basename(html_path), "_files"),
      title = custom_title %||% title %||% resolve_widget_title(widget_to_save) %||% format_identifier_title(basename(filename))
    )
    html_path
  }
  styled_requested <- identical(tolower(Sys.getenv("DEARMETA_STYLED_INTERACTIVE", "")), "true")
  if (!styled_requested) {
    return(tryCatch(
      {
        log_message("Rendering interactive %s using simplified widget output.", filename)
        simple_save(raw_widget)
      },
      error = function(e) {
        log_message("Simple interactive save failed for %s: %s", filename, conditionMessage(e))
        NULL
      }
    ))
  }
  tryCatch(
    {
      plot_class <- paste(class(plot_obj), collapse = "/")
      trace_count <- tryCatch(length(plot_obj$x$data), error = function(e) NA_integer_)
      log_message("Rendering interactive %s (%s) (%s traces)", filename, plot_class, trace_count)
      save_styled_widget(plot_obj, html_path, title = title, subtitle = subtitle, description = description)
    },
    error = function(e) {
      trace <- tryCatch(paste(utils::capture.output(traceback()), collapse = " | "), error = function(...) "No traceback available")
      log_message("Skipping interactive %s due to: %s | Trace: %s", filename, conditionMessage(e), trace)
      log_message("Attempting unstyled interactive fallback for %s", filename)
      fallback <- tryCatch(
        simple_save(raw_widget, custom_title = title),
        error = function(e2) {
          log_message("Fallback interactive save failed for %s: %s", filename, conditionMessage(e2))
          NULL
        }
      )
      fallback
    }
  )
}

write_dashboard_index <- function(project_root, interactive_dir, interactive_files, summary) {
  root_index_path <- file.path(project_root, "index.html")
  interactive_index_path <- file.path(interactive_dir, "index.html")

  gse_label <- summary$gse %||% format_identifier_title(basename(project_root))
  sample_count <- summary$samples %||% 0
  group_counts <- unlist(summary$groups %||% list(), use.names = TRUE)
  group_counts <- group_counts[!is.na(group_counts)]
  group_summary <- if (length(group_counts) > 0) paste(format_number(length(group_counts)), "groups") else "Single group"

  batch_label <- if (isTRUE(summary$combat_applied)) {
    label <- "ComBat"
    if (!is.null(summary$batch_column) && nzchar(summary$batch_column)) {
      label <- paste0(label, " · ", summary$batch_column)
    }
    label
  } else if ((summary$sva_surrogates %||% 0) > 0) {
    paste0("SVA · ", format_number(summary$sva_surrogates), " SVs")
  } else {
    "Baseline"
  }
  batch_methods <- summary$batch_methods %||% list()
  batch_method_parts <- character()
  if (!is.null(batch_methods$minfi) && nzchar(batch_methods$minfi)) {
    batch_method_parts <- c(batch_method_parts, paste0("minfi ", batch_methods$minfi))
  }
  if (!is.null(batch_methods$sesame) && nzchar(batch_methods$sesame)) {
    batch_method_parts <- c(batch_method_parts, paste0("sesame ", batch_methods$sesame))
  }
  batch_detail <- if (length(batch_method_parts) > 0) paste(batch_method_parts, collapse = " · ") else "Batch selection: default"

  stats_cards_data <- list(
    list(label = "Samples analysed", value = format_number(sample_count), detail = group_summary),
    list(label = "Platform", value = summary$platform %||% "Unknown platform", detail = summary$array_type %||% ""),
    list(
      label = "FDR threshold",
      value = paste0("\u2264 ", signif(summary$fdr_threshold %||% 0.05, 3)),
      detail = paste("Top", format_number(summary$top_n_cpgs %||% 0), "CpGs retained")
    ),
    list(label = "Batch strategy", value = batch_label, detail = batch_detail)
  )
  stats_cards <- lapply(stats_cards_data, function(stat) {
    htmltools::tags$div(
      class = "dm-stat-card",
      htmltools::tags$span(class = "dm-stat-card__label", stat$label),
      htmltools::tags$span(class = "dm-stat-card__value", stat$value),
      if (!is.null(stat$detail) && nzchar(stat$detail)) htmltools::tags$span(class = "dm-stat-card__detail", stat$detail)
    )
  })

  group_badges <- NULL
  if (length(group_counts) > 0) {
    group_badges <- htmltools::tags$div(
      class = "dm-chip-row",
      lapply(seq_along(group_counts), function(i) {
        name <- names(group_counts)[i]
        label <- if (is.null(name) || !nzchar(name)) sprintf("Group %s", format_number(i)) else format_identifier_title(name)
        htmltools::tags$span(
          class = "dm-chip",
          htmltools::tags$span(class = "dm-chip__label", label),
          htmltools::tags$span(class = "dm-chip__value", format_number(group_counts[[i]]))
        )
      })
    )
  }

  sig_counts <- unlist(summary$significant_cpgs %||% list(), use.names = TRUE)
  sig_counts <- sig_counts[!is.na(sig_counts)]
  sig_section <- NULL
  if (length(sig_counts) > 0) {
    sig_section <- htmltools::tags$div(
      class = "dm-chip-row dm-chip-row--highlight",
      lapply(seq_along(sig_counts), function(i) {
        name <- names(sig_counts)[i]
        label <- format_identifier_title(name)
        htmltools::tags$span(
          class = "dm-chip dm-chip--accent",
          htmltools::tags$span(class = "dm-chip__label", label),
          htmltools::tags$span(class = "dm-chip__value", format_number(sig_counts[[i]]))
        )
      })
    )
  }

  overlap_sentence <- NULL
  overlap_stats <- summary$overlap_top_n %||% list()
  if (!is.null(overlap_stats$overlap)) {
    overlap_sentence <- sprintf(
      "Top %s CpGs overlap: %s (minfi %s, sesame %s).",
      format_number(max(overlap_stats$top_n_minfi %||% 0, overlap_stats$top_n_sesame %||% 0)),
      format_number(overlap_stats$overlap %||% 0),
      format_percent(overlap_stats$overlap_pct_minfi %||% 0),
      format_percent(overlap_stats$overlap_pct_sesame %||% 0)
    )
  }

  render_chip_list <- function(items, empty_label = "None detected") {
    items <- unique(items)
    items <- items[nzchar(items)]
    if (length(items) == 0) {
      return(htmltools::tags$span(class = "dm-empty-note", empty_label))
    }
    htmltools::tags$div(
      class = "dm-chip-stack",
      lapply(items, function(item) {
        htmltools::tags$span(class = "dm-chip dm-chip--subtle", format_identifier_title(item))
      })
    )
  }

  covariate_section <- htmltools::tags$section(
    class = "dm-section",
    htmltools::tags$h2(class = "dm-section__title", "Design covariates"),
    htmltools::tags$div(
      class = "dm-covariate-grid",
      htmltools::tags$div(
        class = "dm-covariate-card",
        htmltools::tags$h3(class = "dm-covariate-card__title", "Candidate batches"),
        render_chip_list(summary$candidate_batches$final %||% character())
      ),
      htmltools::tags$div(
        class = "dm-covariate-card",
        htmltools::tags$h3(class = "dm-covariate-card__title", "Numeric covariates"),
        render_chip_list(summary$covariates$numeric %||% character())
      ),
      htmltools::tags$div(
        class = "dm-covariate-card",
        htmltools::tags$h3(class = "dm-covariate-card__title", "Factor covariates"),
        render_chip_list(summary$covariates$factor %||% character())
      )
    )
  )

  interactive_sections <- list()
  if (length(interactive_files) > 0) {
    entries <- lapply(names(interactive_files), function(key) {
      meta <- interactive_output_metadata(key)
      href <- relative_to_root(interactive_files[[key]], project_root)
      href <- gsub("\\\\", "/", href)
      icon_node <- htmltools::tags$span(class = "dm-link-card__icon", meta$icon)
      if (identical(key, "table_intersection_dual")) {
        shared_stats <- summary$shared_cpgs %||% list()
        shared_count <- shared_stats$count %||% 0
        total_cap <- summary$top_n_cpgs %||% 0
        if (is.null(total_cap) || !is.finite(total_cap) || total_cap <= 0) {
          overlap_stats <- summary$overlap_top_n %||% list()
          total_cap <- max(
            shared_count,
            overlap_stats$top_n_minfi %||% 0,
            overlap_stats$top_n_sesame %||% 0
          )
        }
        percent_val <- if (total_cap > 0) shared_count / total_cap else NA_real_
        count_label <- format_number(shared_count)
        total_label <- if (total_cap > 0) {
          sprintf("of %s total", format_number(total_cap))
        } else {
          "total n/a"
        }
        percent_label <- if (is.na(percent_val)) {
          "shared n/a"
        } else {
          sprintf("%s shared", format_percent(percent_val))
        }
        icon_node <- htmltools::tags$span(
          class = "dm-link-card__icon dm-link-card__icon--metrics",
          htmltools::tags$span(class = "dm-icon-metric__count", count_label),
          htmltools::tags$span(class = "dm-icon-metric__total", total_label),
          htmltools::tags$span(class = "dm-icon-metric__percent", percent_label)
        )
      }
      list(
        category = meta$category,
        card = htmltools::tags$a(
          class = "dm-link-card",
          href = href,
          target = "_blank",
          icon_node,
          htmltools::tags$div(
            class = "dm-link-card__body",
            htmltools::tags$h3(class = "dm-link-card__title", meta$title),
            if (!is.null(meta$description) && nzchar(meta$description)) htmltools::tags$p(class = "dm-link-card__desc", meta$description),
            htmltools::tags$span(class = "dm-link-card__cta", "Open")
          )
        )
      )
    })
    entries <- Filter(function(item) !is.null(item$card), entries)
    if (length(entries) > 0) {
      categories <- split(entries, vapply(entries, function(item) item$category %||% "Additional Outputs", character(1)))
      category_order <- c("Quality Control", "Differential Analysis", "Genomic Context", "Annotated Tables", "Additional Outputs")
      ordered_names <- names(categories)[order(match(names(categories), category_order, nomatch = length(category_order) + 1))]
      interactive_sections <- lapply(ordered_names, function(name) {
        cards <- lapply(categories[[name]], function(item) item$card)
        if (length(cards) == 0) {
          return(NULL)
        }
        htmltools::tags$div(
          class = "dm-interactive-group",
          htmltools::tags$h3(class = "dm-interactive-group__title", name),
          htmltools::tags$div(class = "dm-link-grid", cards)
        )
      })
      interactive_sections <- Filter(Negate(is.null), interactive_sections)
    }
  }

  interactive_block <- NULL
  if (length(interactive_sections) > 0) {
    interactive_block <- htmltools::tags$section(
      class = "dm-section",
      htmltools::tags$h2(class = "dm-section__title", "Interactive reports"),
      htmltools::tags$div(
        class = "dm-section-stack",
        do.call(htmltools::tagList, interactive_sections)
      )
    )
  }

  render_top_cpg_section <- function(label, rows) {
    if (is.null(rows) || length(rows) == 0) {
      return(NULL)
    }
    dt <- data.table::rbindlist(lapply(rows, data.table::as.data.table), fill = TRUE)
    if (nrow(dt) == 0) {
      return(NULL)
    }
    dt <- as.data.frame(dt, stringsAsFactors = FALSE)
    if (nrow(dt) > 10) {
      dt <- dt[seq_len(10), , drop = FALSE]
    }
    column_map <- c(
      probe_id = "CpG ID",
      delta_beta = "\u0394\u03b2",
      adj.P.Val = "FDR",
      P.Value = "p-value",
      gene_symbols = "Gene",
      gene_regions = "Region",
      logFC = "logFC"
    )
    cols <- intersect(names(column_map), colnames(dt))
    if (length(cols) == 0) {
      return(NULL)
    }
    for (col in cols) {
      if (col %in% c("delta_beta", "logFC")) {
        dt[[col]] <- sprintf("%.3f", as.numeric(dt[[col]]))
      } else if (col %in% c("adj.P.Val", "P.Value")) {
        dt[[col]] <- sprintf("%.3g", as.numeric(dt[[col]]))
      } else {
        dt[[col]] <- as.character(dt[[col]])
      }
    }
    header <- htmltools::tags$tr(lapply(cols, function(col) htmltools::tags$th(column_map[[col]])))
    body <- lapply(seq_len(nrow(dt)), function(i) {
      htmltools::tags$tr(lapply(cols, function(col) htmltools::tags$td(dt[[col]][i])))
    })
    htmltools::tags$div(
      class = "dm-table-card",
      htmltools::tags$h3(class = "dm-table-card__title", label),
      htmltools::tags$div(
        class = "dm-table-card__scroll",
        htmltools::tags$table(
          class = "dm-table-card__table",
          htmltools::tags$thead(header),
          htmltools::tags$tbody(body)
        )
      )
    )
  }

  top_cpg_sections <- Filter(
    Negate(is.null),
    list(
      render_top_cpg_section("Top CpGs · minfi", summary$top_cpgs$minfi),
      render_top_cpg_section("Top CpGs · sesame", summary$top_cpgs$sesame),
      render_top_cpg_section("Top CpGs · intersection", summary$top_cpgs$intersection)
    )
  )
  top_cpg_section <- NULL
  if (length(top_cpg_sections) > 0) {
    top_cpg_section <- htmltools::tags$section(
      class = "dm-section",
      htmltools::tags$h2(class = "dm-section__title", "Highlighted CpGs"),
      htmltools::tags$div(class = "dm-table-grid", top_cpg_sections)
    )
  }

  directories <- summary$directories %||% list()
  directory_cards <- list()
  if (length(directories) > 0) {
    directory_cards <- lapply(names(directories), function(name) {
      if (identical(name, "interactive")) {
        return(NULL)
      }
      dir_path <- directories[[name]]
      if (is.null(dir_path) || !nzchar(dir_path)) {
        return(NULL)
      }
      href <- relative_to_root(dir_path, project_root)
      href <- gsub("\\\\", "/", href)
      entries <- summary$files[[name]] %||% list()
      count <- length(entries)
      total_size <- sum(vapply(entries, function(item) item$size %||% 0, numeric(1)), na.rm = TRUE)
      htmltools::tags$a(
        class = "dm-download-card",
        href = href,
        htmltools::tags$h3(class = "dm-download-card__title", format_identifier_title(name)),
        htmltools::tags$p(class = "dm-download-card__meta", sprintf("%s items · %s", format_number(count), format_bytes(total_size))),
        htmltools::tags$span(class = "dm-download-card__cta", "Open folder")
      )
    })
  }

  runtime_files <- summary$runtime_files %||% list()
  runtime_cards <- list()
  if (length(runtime_files) > 0) {
    runtime_cards <- lapply(names(runtime_files), function(name) {
      file_path <- runtime_files[[name]]
      if (is.null(file_path) || !nzchar(file_path)) {
        return(NULL)
      }
      href <- relative_to_root(file_path, project_root)
      href <- gsub("\\\\", "/", href)
      htmltools::tags$a(
        class = "dm-download-card dm-download-card--file",
        href = href,
        target = "_blank",
        htmltools::tags$h3(class = "dm-download-card__title", format_identifier_title(name)),
        htmltools::tags$p(class = "dm-download-card__meta", "Runtime artifact"),
        htmltools::tags$span(class = "dm-download-card__cta", "View file")
      )
    })
  }

  resource_cards <- Filter(Negate(is.null), c(directory_cards, runtime_cards))
  resource_section <- NULL
  if (length(resource_cards) > 0) {
    resource_section <- htmltools::tags$section(
      class = "dm-section",
      htmltools::tags$h2(class = "dm-section__title", "Artifacts"),
      htmltools::tags$div(class = "dm-download-grid", resource_cards)
    )
  }

  hero_description <- sprintf(
    "%s samples processed on %s with FDR threshold %.3f (top %s CpGs retained).",
    format_number(sample_count),
    summary$platform %||% "unknown platform",
    as.numeric(summary$fdr_threshold %||% 0.05),
    format_number(summary$top_n_cpgs %||% 0)
  )

  dashboard_body <- htmltools::tagList(
    htmltools::tags$header(
      class = "dm-hero",
      htmltools::tags$p(class = "dm-hero__eyebrow", "DearMeta · Interactive workspace"),
      htmltools::tags$h1(class = "dm-hero__title", gse_label),
      htmltools::tags$p(class = "dm-hero__subtitle", hero_description)
    ),
    htmltools::tags$div(class = "dm-stats-grid", stats_cards),
    group_badges,
    sig_section,
    if (!is.null(overlap_sentence)) htmltools::tags$p(class = "dm-text-note", overlap_sentence),
    covariate_section,
    interactive_block,
    top_cpg_section,
    resource_section,
    htmltools::tags$footer(
      class = "dm-footer",
      "Generated with DearMeta"
    )
  )

  css <- "
:root { color-scheme: dark; }
* { box-sizing: border-box; }
body {
  margin: 0;
  font-family: 'Inter', 'Segoe UI', 'Helvetica Neue', sans-serif;
  background: radial-gradient(120% 160% at 50% 0%, rgba(59,130,246,0.22) 0%, rgba(15,23,42,1) 60%);
  color: #E2E8F0;
}
a { color: inherit; text-decoration: none; }
.dm-dashboard {
  max-width: 1360px;
  margin: 0 auto;
  padding: 68px 32px 86px;
  display: flex;
  flex-direction: column;
  gap: 48px;
}
.dm-hero__eyebrow {
  text-transform: uppercase;
  letter-spacing: 0.22em;
  color: rgba(148, 163, 184, 0.72);
  font-size: 0.78rem;
  margin: 0;
}
.dm-hero__title {
  margin: 0;
  font-size: clamp(2.4rem, 5vw, 3.3rem);
  font-weight: 700;
  color: #F8FAFC;
}
.dm-hero__subtitle {
  margin: 12px 0 0;
  font-size: 1.08rem;
  color: rgba(226, 232, 240, 0.82);
  max-width: 880px;
  line-height: 1.65;
}
.dm-stats-grid {
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(260px, 1fr));
  gap: 22px;
}
.dm-stat-card {
  background: rgba(15, 23, 42, 0.82);
  border: 1px solid rgba(148, 163, 184, 0.18);
  border-radius: 22px;
  padding: 22px 26px;
  display: flex;
  flex-direction: column;
  gap: 6px;
  box-shadow: 0 32px 68px rgba(15, 23, 42, 0.45);
}
.dm-stat-card__label {
  color: rgba(148, 163, 184, 0.82);
  font-size: 0.8rem;
  letter-spacing: 0.08em;
  text-transform: uppercase;
}
.dm-stat-card__value {
  font-size: 1.9rem;
  font-weight: 700;
  color: #F8FAFC;
}
.dm-stat-card__detail {
  color: rgba(203, 213, 225, 0.82);
  font-size: 0.95rem;
}
.dm-chip-row {
  display: flex;
  flex-wrap: wrap;
  gap: 12px;
}
.dm-chip {
  display: inline-flex;
  align-items: center;
  gap: 6px;
  padding: 8px 14px;
  border-radius: 999px;
  background: rgba(30, 41, 59, 0.78);
  border: 1px solid rgba(148, 163, 184, 0.22);
  font-size: 0.85rem;
}
.dm-chip--accent {
  background: linear-gradient(135deg, rgba(59, 130, 246, 0.6), rgba(139, 92, 246, 0.6));
  border: none;
}
.dm-chip--subtle {
  background: rgba(30, 41, 59, 0.62);
}
.dm-chip__label {
  text-transform: uppercase;
  letter-spacing: 0.08em;
  font-size: 0.72rem;
  color: rgba(226, 232, 240, 0.9);
}
.dm-chip__value {
  font-weight: 600;
  font-size: 0.95rem;
  color: #F8FAFC;
}
.dm-chip-row--highlight {
  margin-top: -6px;
}
.dm-section {
  display: flex;
  flex-direction: column;
  gap: 22px;
}
.dm-section-stack {
  display: flex;
  flex-direction: column;
  gap: 24px;
}
.dm-interactive-group {
  display: flex;
  flex-direction: column;
  gap: 16px;
}
.dm-interactive-group__title {
  margin: 0;
  font-size: 1.2rem;
  font-weight: 600;
  color: rgba(226, 232, 240, 0.92);
  letter-spacing: 0.04em;
  text-transform: uppercase;
}
.dm-section__title {
  margin: 0;
  font-size: 1.55rem;
  font-weight: 600;
  color: #F8FAFC;
}
.dm-text-note {
  margin: 0;
  font-size: 0.95rem;
  color: rgba(226, 232, 240, 0.78);
}
.dm-covariate-grid {
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(260px, 1fr));
  gap: 22px;
}
.dm-covariate-card {
  background: rgba(15, 23, 42, 0.78);
  border: 1px solid rgba(148, 163, 184, 0.18);
  border-radius: 20px;
  padding: 22px 24px;
  box-shadow: 0 18px 40px rgba(15, 23, 42, 0.4);
}
.dm-covariate-card__title {
  margin: 0 0 12px;
  font-size: 1.05rem;
  color: #E2E8F0;
}
.dm-chip-stack {
  display: flex;
  flex-wrap: wrap;
  gap: 8px;
}
.dm-empty-note {
  color: rgba(148, 163, 184, 0.72);
  font-size: 0.9rem;
}
.dm-link-grid {
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(280px, 1fr));
  gap: 22px;
}
.dm-link-card {
  display: flex;
  gap: 16px;
  background: rgba(15, 23, 42, 0.78);
  border: 1px solid rgba(148, 163, 184, 0.22);
  border-radius: 22px;
  padding: 22px 26px;
  box-shadow: 0 28px 54px rgba(15, 23, 42, 0.4);
  transition: transform 0.18s ease, box-shadow 0.18s ease, border-color 0.18s ease;
}
.dm-link-card:hover {
  transform: translateY(-4px);
  border-color: rgba(96, 165, 250, 0.6);
  box-shadow: 0 32px 60px rgba(59, 130, 246, 0.35);
}
.dm-link-card__icon {
  font-size: 0.76rem;
  font-weight: 700;
  color: rgba(148, 163, 184, 0.85);
  text-transform: uppercase;
  letter-spacing: 0.12em;
  padding-top: 4px;
}
.dm-link-card__icon--metrics {
  display: flex;
  flex-direction: column;
  gap: 2px;
  align-items: flex-start;
  padding-top: 0;
  text-transform: none;
  letter-spacing: 0.02em;
  color: rgba(203, 213, 225, 0.78);
}
.dm-link-card__icon--metrics .dm-icon-metric__count {
  font-size: 1.05rem;
  font-weight: 700;
  color: #F8FAFC;
}
.dm-link-card__icon--metrics .dm-icon-metric__total {
  font-size: 0.78rem;
  color: rgba(203, 213, 225, 0.78);
}
.dm-link-card__icon--metrics .dm-icon-metric__percent {
  font-size: 0.78rem;
  font-weight: 600;
  color: #38BDF8;
}
.dm-link-card__body {
  display: flex;
  flex-direction: column;
  gap: 8px;
}
.dm-link-card__title {
  margin: 0;
  font-size: 1.05rem;
  color: #F8FAFC;
}
.dm-link-card__desc {
  margin: 0;
  font-size: 0.92rem;
  color: rgba(203, 213, 225, 0.78);
  line-height: 1.5;
}
.dm-link-card__cta {
  font-size: 0.85rem;
  font-weight: 600;
  color: #38BDF8;
}
.dm-table-grid {
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(320px, 1fr));
  gap: 22px;
}
.dm-table-card {
  background: rgba(15, 23, 42, 0.78);
  border: 1px solid rgba(148, 163, 184, 0.2);
  border-radius: 20px;
  padding: 22px 24px;
}
.dm-table-card__scroll {
  width: 100%;
  overflow-x: auto;
}
.dm-table-card__title {
  margin: 0 0 12px;
  font-size: 1.05rem;
  color: #F8FAFC;
}
.dm-table-card__table {
  width: 100%;
  border-collapse: collapse;
  font-size: 0.85rem;
  table-layout: fixed;
}
.dm-table-card__table thead th {
  text-align: left;
  padding: 6px 8px;
  color: rgba(148, 163, 184, 0.9);
  border-bottom: 1px solid rgba(148, 163, 184, 0.35);
}
.dm-table-card__table tbody td {
  padding: 6px 8px;
  border-bottom: 1px solid rgba(148, 163, 184, 0.1);
  color: rgba(226, 232, 240, 0.9);
}
.dm-download-grid {
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(280px, 1fr));
  gap: 22px;
}
.dm-download-card {
  background: rgba(15, 23, 42, 0.78);
  border: 1px solid rgba(148, 163, 184, 0.22);
  border-radius: 20px;
  padding: 22px 24px;
  display: flex;
  flex-direction: column;
  gap: 10px;
  transition: transform 0.18s ease, border-color 0.18s ease;
}
.dm-download-card:hover {
  transform: translateY(-3px);
  border-color: rgba(96, 165, 250, 0.65);
}
.dm-download-card__title {
  margin: 0;
  font-size: 1.05rem;
  color: #F8FAFC;
}
.dm-download-card__meta {
  margin: 0;
  font-size: 0.85rem;
  color: rgba(148, 163, 184, 0.9);
}
.dm-download-card__cta {
  font-size: 0.85rem;
  font-weight: 600;
  color: #38BDF8;
}
.dm-download-card--file {
  border-style: dashed;
}
.dm-footer {
  text-align: center;
  font-size: 0.88rem;
  color: rgba(148, 163, 184, 0.72);
  padding-top: 16px;
}
@media (max-width: 720px) {
  .dm-dashboard {
    padding: 44px 18px 64px;
    gap: 34px;
  }
  .dm-link-card {
    flex-direction: column;
    align-items: flex-start;
  }
  .dm-link-card__icon {
    padding-top: 0;
  }
}
"

  html_doc <- htmltools::tags$html(
    lang = "en",
    htmltools::tags$head(
      htmltools::tags$meta(charset = "utf-8"),
      htmltools::tags$meta(name = "viewport", content = "width=device-width, initial-scale=1"),
      htmltools::tags$title(sprintf("DearMeta Dashboard · %s", gse_label)),
      htmltools::tags$style(htmltools::HTML(css))
    ),
    htmltools::tags$body(
      htmltools::tags$div(class = "dm-dashboard", dashboard_body)
    )
  )

  htmltools::save_html(htmltools::browsable(html_doc), file = root_index_path)

  redirect_target <- file.path("..", basename(root_index_path))
  redirect_href <- htmltools::htmlEscape(gsub("\\\\", "/", redirect_target))
  redirect_html <- sprintf(
    '<html><head><meta charset="utf-8"><meta http-equiv="refresh" content="0; url=%1$s"></head><body style="background:#0f172a;color:#e2e8f0;font-family:Arial, sans-serif;display:flex;align-items:center;justify-content:center;height:100vh;"><p>Open the <a style="color:#38BDF8;" href="%1$s">DearMeta dashboard</a>.</p></body></html>',
    redirect_href
  )
  writeLines(redirect_html, con = interactive_index_path)
  log_message("Dashboard index created at %s (interactive redirect at %s)", root_index_path, interactive_index_path)
  list(root = root_index_path, interactive = interactive_index_path)
}

# ---- Load configuration -----------------------------------------------------

config <- read_configure(config_path, project_root)
group_counts <- table(config$dear_group)
if (any(group_counts < opt$min_group_size)) {
  stop("Insufficient samples per group: ", paste(names(group_counts[group_counts < opt$min_group_size]), collapse = ", "))
}

log_message("Loaded", nrow(config), "samples for analysis across groups:", paste(names(group_counts), group_counts, sep = "=", collapse = "; "))

# Determine candidate batches and covariates
directives <- attr(config, "directives")
if (is.null(directives)) {
  directives <- list(batch = character(), protect = character(), drop = character(), numeric = character(), factor = character())
}

manual_batch_columns <- match_columns(names(config), directives$batch)
if (length(directives$batch) > 0 && length(manual_batch_columns) < length(directives$batch)) {
  missing <- setdiff(directives$batch, manual_batch_columns)
  log_message("Requested batch columns not found in configure.tsv: %s", paste(missing, collapse = ", "))
}
if (length(manual_batch_columns) > 0) {
  filtered_manual <- remove_group_duplicates(manual_batch_columns, config, log_prefix = "Batch candidate")
  dropped_manual <- setdiff(manual_batch_columns, filtered_manual)
  if (length(dropped_manual) > 0) {
    log_message("Excluded manual batch columns duplicating dear_group: %s", paste(dropped_manual, collapse = ", "))
  }
  manual_batch_columns <- filtered_manual
}

protected_columns <- unique(c(match_columns(names(config), DEFAULT_PROTECTED_BATCH), match_columns(names(config), directives$protect)))
if (length(directives$protect) > 0 && length(match_columns(names(config), directives$protect)) < length(directives$protect)) {
  missing_protect <- setdiff(directives$protect, match_columns(names(config), directives$protect))
  if (length(missing_protect) > 0) {
    log_message("Requested protected columns not found in configure.tsv: %s", paste(missing_protect, collapse = ", "))
  }
}
protected_for_detection <- setdiff(protected_columns, manual_batch_columns)

candidate_batches <- detect_candidate_batches(config, protected_for_detection)
auto_batches <- candidate_batches
if (length(manual_batch_columns) > 0) {
  candidate_batches <- unique(c(manual_batch_columns, candidate_batches))
  log_message("Including manual batch columns: %s", paste(manual_batch_columns, collapse = ", "))
}
if (length(auto_batches) > 0) {
  auto_filtered <- remove_group_duplicates(auto_batches, config, log_prefix = "Batch candidate")
  dropped_auto <- setdiff(auto_batches, auto_filtered)
  if (length(dropped_auto) > 0) {
    log_message("Removed auto-detected batch columns duplicating dear_group: %s", paste(dropped_auto, collapse = ", "))
  }
  auto_batches <- auto_filtered
}
if (length(candidate_batches) > 0) {
  candidate_filtered <- remove_group_duplicates(candidate_batches, config, log_prefix = "Batch candidate")
  candidate_dropped <- setdiff(candidate_batches, candidate_filtered)
  if (length(candidate_dropped) > 0) {
    log_message("Removed batch columns duplicating dear_group from final candidate list: %s", paste(candidate_dropped, collapse = ", "))
  }
  candidate_batches <- candidate_filtered
}

covars <- split_covariates(config, candidate_batches)
dropped_covariates <- attr(covars, "dropped_covariates")
if (length(dropped_covariates) > 0) {
  log_message(
    "Excluded covariates lacking usable variation: %s",
    paste(sprintf("%s (%s)", names(dropped_covariates), unlist(dropped_covariates)), collapse = ", ")
  )
} else {
  dropped_covariates <- list()
}

manual_numeric <- match_columns(names(config), directives$numeric)
if (length(directives$numeric) > 0 && length(manual_numeric) < length(directives$numeric)) {
  missing_numeric <- setdiff(directives$numeric, manual_numeric)
  if (length(missing_numeric) > 0) {
    log_message("Requested numeric covariates not found in configure.tsv: %s", paste(missing_numeric, collapse = ", "))
  }
}
if (length(manual_numeric) > 0) {
  covars$numeric <- unique(c(covars$numeric, manual_numeric))
  log_message("Force-included numeric covariates: %s", paste(manual_numeric, collapse = ", "))
}
manual_factor <- match_columns(names(config), directives$factor)
if (length(directives$factor) > 0 && length(manual_factor) < length(directives$factor)) {
  missing_factor <- setdiff(directives$factor, manual_factor)
  if (length(missing_factor) > 0) {
    log_message("Requested factor covariates not found in configure.tsv: %s", paste(missing_factor, collapse = ", "))
  }
}
if (length(manual_factor) > 0) {
  covars$factor <- unique(c(covars$factor, manual_factor))
  log_message("Force-included factor covariates: %s", paste(manual_factor, collapse = ", "))
}
attr(covars, "dropped_covariates") <- NULL
batch_tracking <- list(auto = auto_batches, manual = manual_batch_columns, protected = protected_columns)
manual_covariates <- list(numeric = manual_numeric, factor = manual_factor)

platform_values <- unique(config$platform_version)
if (length(platform_values) != 1) {
  stop("Mixed platform_version values detected; DearMeta expects a single array type per run.")
}
platform_label <- platform_values[1]

array_type_map <- list(
  EPICv1 = list(
    dmr = "EPICv1",
    sesame_platform = "EPIC",
    sesame = c("EPIC.address", "EPIC.hg38.manifest", "EPIC.hg19.manifest")
  ),
  EPICv2 = list(
    dmr = "EPICv2",
    sesame_platform = "EPICv2",
    sesame = c("EPICv2.address", "EPICv2.hg38.manifest", "EPICv2.hg19.manifest")
  ),
  HM450 = list(
    dmr = "450K",
    sesame_platform = "HM450",
    sesame = c("HM450.address", "HM450.hg38.manifest", "HM450.hg19.manifest")
  )
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
sesame_platform <- platform_info$sesame_platform %||% platform_info$dmr %||% platform_label
sesame_candidates <- unique(c(
  as.character(platform_info$sesame %||% character()),
  paste0(sesame_platform, ".address"),
  paste0(platform_label, ".address"),
  paste0(sesame_platform, ".hg38.manifest"),
  paste0(sesame_platform, ".hg19.manifest")
))
manifest_info <- load_sesame_manifest(sesame_candidates, verbose = FALSE)
sesame_manifest_source <- NULL
if (!is.null(manifest_info$ordering)) {
  log_message("Loaded sesame manifest from %s", manifest_info$source)
  sesame_manifest <- manifest_info$ordering
  sesame_controls <- manifest_info$controls
  sesame_manifest_source <- manifest_info$source
  tryCatch(
    {
      for (i in seq_len(nrow(config))) {
        base <- basenames[i]
        sample_id <- config$gsm_id[i]
        # Apply a sesame recipe that mirrors the original intent (Noob + pOOBAH) while
        # using the official one-letter codes expected by sesame.
        # Force SerialParam to avoid pthread errors on systems that cannot spawn additional workers.
        original_alt_repo <- getOption("SESAMEDATA_USE_ALT")
        attempt <- 1
        sig <- NULL
        suppressWarnings({
          while (attempt <= 2) {
            options(SESAMEDATA_USE_ALT = if (attempt == 1) original_alt_repo else FALSE)
            sig <- tryCatch(
              openSesame(
                base,
                prep = "BP",
                platform = sesame_platform,
                manifest = sesame_manifest,
                controls = sesame_controls,
                BPPARAM = BiocParallel::SerialParam(),
                func = NULL
              ),
              error = function(e) e
            )
            if (!inherits(sig, "error")) {
              break
            }
            err_msg <- conditionMessage(sig)
            if (attempt == 1 && isTRUE(original_alt_repo) &&
              grepl("(cannot be retrieved|찾을 수 없습니다)", err_msg, ignore.case = TRUE)) {
              log_message(
                "Retrying sesame for sample %s without alternate data repo due to missing resource: %s",
                sample_id,
                err_msg
              )
              attempt <- attempt + 1
              next
            }
            stop(
              sprintf("Sesame processing failed for sample %s: %s", sample_id, err_msg),
              call. = FALSE
            )
          }
        })
        options(SESAMEDATA_USE_ALT = original_alt_repo)
        if (inherits(sig, "error")) {
          stop(
            sprintf("Sesame processing failed for sample %s: %s", sample_id, conditionMessage(sig)),
            call. = FALSE
          )
        }
        sesame_betas[[sample_id]] <- getBetas(sig)
        sesame_detp[[sample_id]] <- pOOBAH(sig, return.pval = TRUE)
      }
    },
    error = function(e) {
      sesame_available <<- FALSE
      sesame_error <<- conditionMessage(e)
    }
  )
} else {
  sesame_available <- FALSE
  attempted <- names(manifest_info$attempts)
  if (length(attempted) > 0) {
    log_message(
      "Failed to load sesame manifest for %s; attempted: %s",
      platform_label,
      paste(unique(attempted), collapse = ", ")
    )
    detail <- paste(
      sprintf("%s (%s)", names(manifest_info$attempts), unlist(manifest_info$attempts, use.names = FALSE)),
      collapse = "; "
    )
    if (nzchar(detail)) {
      log_message("Sesame manifest lookup details: %s", detail)
    }
  } else if (length(sesame_candidates) > 0) {
    log_message(
      "Failed to load sesame manifest for %s; attempted fallback list: %s",
      platform_label,
      paste(sesame_candidates, collapse = ", ")
    )
  } else {
    log_message("No sesame manifest candidates defined for platform %s.", platform_label)
  }
  sesame_error <- sprintf("Manifest unavailable for sesame platform %s.", sesame_platform)
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
  if (!is.matrix(M_matrix)) {
    M_matrix <- as.matrix(M_matrix)
  }
  finite_mask <- apply(M_matrix, 1, function(x) all(is.finite(x)))
  M_matrix <- M_matrix[finite_mask, , drop = FALSE]
  if (nrow(M_matrix) == 0 || ncol(M_matrix) < 2) {
    return(empty_pca_result())
  }
  var_mask <- apply(M_matrix, 1, function(x) {
    vals <- x[is.finite(x)]
    if (length(vals) <= 1) {
      return(FALSE)
    }
    sd(vals) > 0
  })
  M_matrix <- M_matrix[var_mask, , drop = FALSE]
  if (nrow(M_matrix) == 0 || ncol(M_matrix) < 2) {
    return(empty_pca_result())
  }
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

metadata_dt <- as.data.table(config)
variables_to_check <- unique(c("dear_group", candidate_batches, covars$numeric, covars$factor))

pca_minfi_pre <- assess_pca(M_minfi, metadata_dt, variables_to_check, "minfi_pre")

if (sesame_available && !is.null(M_sesame)) {
  pca_sesame_pre <- assess_pca(M_sesame, metadata_dt, variables_to_check, "sesame_pre")
} else {
  pca_sesame_pre <- empty_pca_result()
}

design_optimisation <- tryCatch(
  optimize_design_matrix(
    metadata = metadata_dt,
    covariates = covars,
    manual_covariates = manual_covariates,
    protected_columns = protected_columns,
    candidate_batches = candidate_batches,
    manual_batch_columns = manual_batch_columns,
    M_minfi = M_minfi,
    M_sesame = M_sesame,
    sesame_available = sesame_available,
    opt = opt
  ),
  error = function(e) {
    log_message("Design optimisation errored: %s", conditionMessage(e))
    NULL
  }
)

used_fallback_design <- FALSE
if (is.null(design_optimisation) || is.null(design_optimisation$best)) {
  design_optimisation <- tryCatch(
    fallback_design_selection(
      metadata = metadata_dt,
      covariates = covars,
      manual_covariates = manual_covariates,
      protected_columns = protected_columns,
      candidate_batches = candidate_batches,
      M_minfi = M_minfi,
      M_sesame = M_sesame,
      sesame_available = sesame_available,
      opt = opt
    ),
    error = function(e) {
      stop("Design selection failed: ", conditionMessage(e))
    }
  )
  used_fallback_design <- TRUE
}

best_model <- design_optimisation$best
best_params <- best_model$params
best_metrics <- best_model$metrics
best_outputs <- best_model$outputs
model_evaluations <- design_optimisation$evaluated

selected_covariates <- best_outputs$covariates_final
covars$numeric <- selected_covariates$numeric
covars$factor <- selected_covariates$factor

design_without_sv <- best_outputs$design_base
surrogate_vars <- best_outputs$surrogate_vars
design_for_limma <- best_outputs$design_with_sv
if (ncol(design_for_limma) >= 2) {
  colnames(design_for_limma)[2] <- "groupdear_group"
}

M_minfi_corrected <- best_outputs$M_minfi
beta_minfi <- best_outputs$beta_minfi
save_matrix(beta_minfi, file.path(preprocess_dir, "beta_minfi.tsv.gz"))
save_matrix(M_minfi_corrected, file.path(preprocess_dir, "mvals_minfi.tsv.gz"))

if (sesame_available && !is.null(best_outputs$M_sesame) && !is.null(best_outputs$beta_sesame)) {
  M_sesame_corrected <- best_outputs$M_sesame
  beta_sesame <- best_outputs$beta_sesame
  save_matrix(beta_sesame, file.path(preprocess_dir, "beta_sesame.tsv.gz"))
  save_matrix(M_sesame_corrected, file.path(preprocess_dir, "mvals_sesame.tsv.gz"))
  pca_sesame_post <- best_outputs$pca_sesame_post
} else {
  M_sesame_corrected <- NULL
  beta_sesame <- NULL
  pca_sesame_post <- empty_pca_result()
}

pca_minfi_post <- best_outputs$pca_minfi_post

results_minfi <- best_outputs$results_minfi
results_sesame <- best_outputs$results_sesame

batch_to_use <- best_params$batch_col
surrogate_count <- best_outputs$n_surrogates

selected_minfi_method <-
  if (isTRUE(best_params$use_combat) && isTRUE(best_params$use_sva)) {
    "combat+sva"
  } else if (isTRUE(best_params$use_combat)) {
    "combat"
  } else if (isTRUE(best_params$use_sva)) {
    "sva"
  } else {
    "none"
  }

selected_sesame_method <- if (sesame_available && !is.null(M_sesame_corrected)) selected_minfi_method else "not_available"

log_message(
  "Automated design selected (model %s): batch=%s, combat=%s, sva=%s, numeric covariates=%s, factor covariates=%s",
  best_params$model_id,
  batch_to_use %||% "none",
  best_params$use_combat,
  best_params$use_sva,
  if (length(covars$numeric) == 0) "none" else paste(covars$numeric, collapse = ", "),
  if (length(covars$factor) == 0) "none" else paste(covars$factor, collapse = ", ")
)

log_message("Running limma differential methylation with selected design...")

design_selection_summary <- list(
  strategy = if (used_fallback_design || (!is.null(design_optimisation$note) && identical(design_optimisation$note, "fallback"))) "fallback" else "optimisation",
  status = design_optimisation$status %||% "ok",
  model_id = best_params$model_id,
  metrics = best_metrics,
  evaluated = model_evaluations,
  covariate_sets = design_optimisation$covariate_sets,
  covariate_stats = design_optimisation$covariate_stats,
  batch_options = design_optimisation$batch_options
)

model_evaluations_ranked <- NULL
if (!is.null(model_evaluations) && nrow(model_evaluations) > 0) {
  eval_dt <- as.data.table(model_evaluations)
  eval_dt[, n_covariates_total := n_covariates_numeric + n_covariates_factor]
  eval_dt_ok <- eval_dt[status == "ok"]
  if (nrow(eval_dt_ok) > 0) {
    eval_dt_ok[, batch_ok := is.na(batch_median_p_minfi) | batch_median_p_minfi >= 0.2]
    eval_dt_ok[, batch_score := ifelse(is.na(batch_median_p_minfi), 0, batch_median_p_minfi)]
    setorder(
      eval_dt_ok,
      -batch_ok,
      -n_sig_intersection,
      -n_sig_minfi,
      -n_sig_sesame,
      -batch_score,
      n_covariates_total,
      n_surrogates,
      model_id
    )
    eval_dt_ok[, rank := seq_len(.N)]
    eval_dt <- eval_dt[eval_dt_ok[, .(model_id, rank)], on = "model_id"]
  } else {
    eval_dt[, rank := NA_integer_]
  }
  model_evaluations_ranked <- eval_dt
  design_selection_summary$ranking <- model_evaluations_ranked
}

sig_minfi <- filter_significant(results_minfi, opt$fdr_threshold, opt$delta_beta_threshold)
sig_sesame <- filter_significant(results_sesame, opt$fdr_threshold, opt$delta_beta_threshold)

integrated <- merge_results(sig_minfi, sig_sesame)
integrated$union <- limit_top_hits(integrated$union, opt$top_n_cpgs)
integrated$intersection <- limit_top_hits(integrated$intersection, opt$top_n_cpgs)

analysis_dir <- paths$analysis
model_evaluation_table_path <- NULL
model_evaluation_table_export <- NULL
if (!is.null(model_evaluations_ranked) && nrow(model_evaluations_ranked) > 0) {
  model_evaluation_table_path <- file.path(analysis_dir, "model_evaluation_summary.tsv.gz")
  write_tsv_gz(model_evaluations_ranked, model_evaluation_table_path)
  log_message("Model evaluation summary saved to %s", model_evaluation_table_path)
  model_evaluation_table_export <- copy(model_evaluations_ranked)
}
write_tsv_gz(results_minfi, file.path(analysis_dir, "dmp_minfi.tsv.gz"))
write_tsv_gz(results_sesame, file.path(analysis_dir, "dmp_sesame.tsv.gz"))
write_tsv_gz(integrated$intersection, file.path(analysis_dir, "cpg_intersection.tsv.gz"))
write_tsv_gz(integrated$union, file.path(analysis_dir, "cpg_union.tsv.gz"))

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
write_tsv_gz(as.data.frame(dmr_ranges), file.path(analysis_dir, "dmr_minfi.tsv.gz"))

if (sesame_available && !is.null(M_sesame_corrected) && nrow(sig_sesame) > 0) {
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
  write_tsv_gz(as.data.frame(dmr_ranges_sesame), file.path(analysis_dir, "dmr_sesame.tsv.gz"))
} else {
  if (sesame_available && !is.null(M_sesame_corrected)) {
    log_message("Skipping DMRcate for sesame: no significant CpGs detected at FDR threshold.")
  }
  dmr_ranges_sesame <- data.frame(seqnames = character(), start = integer(), end = integer(), width = integer(), n.sites = integer(), meanstat = numeric())
  write_tsv_gz(dmr_ranges_sesame, file.path(analysis_dir, "dmr_sesame.tsv.gz"))
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
write_tsv_gz(annotated_minfi, file.path(analysis_dir, "annotated_cpg_minfi.tsv.gz"))
write_tsv_gz(annotated_sesame, file.path(analysis_dir, "annotated_cpg_sesame.tsv.gz"))
write_tsv_gz(annotated_intersection, file.path(analysis_dir, "annotated_cpg_intersection.tsv.gz"))

shared_probe_ids <- intersect(annotated_minfi$probe_id, annotated_sesame$probe_id)
annotated_shared <- if (length(shared_probe_ids) > 0) {
  minfi_shared <- annotated_minfi[
    probe_id %in% shared_probe_ids,
    .(
      probe_id,
      chr,
      mapinfo,
      gene_symbols,
      gene_regions,
      relation_to_island,
      logFC_minfi = logFC,
      delta_beta_minfi = delta_beta,
      adj.P.Val_minfi = adj.P.Val,
      P.Value_minfi = P.Value,
      direction_minfi = direction
    )
  ]
  sesame_shared <- annotated_sesame[
    probe_id %in% shared_probe_ids,
    .(
      probe_id,
      logFC_sesame = logFC,
      delta_beta_sesame = delta_beta,
      adj.P.Val_sesame = adj.P.Val,
      P.Value_sesame = P.Value,
      direction_sesame = direction
    )
  ]
  shared <- merge(minfi_shared, sesame_shared, by = "probe_id", all = FALSE, sort = FALSE)
  if (nrow(shared) > 0) {
    shared[, adj.P.Val_min := pmin(adj.P.Val_minfi, adj.P.Val_sesame)]
    shared[, adj.P.Val_max := pmax(adj.P.Val_minfi, adj.P.Val_sesame)]
    shared[, P.Value_min := pmin(P.Value_minfi, P.Value_sesame)]
    shared[, logFC_mean := (logFC_minfi + logFC_sesame) / 2]
    shared[, delta_beta_diff := delta_beta_minfi - delta_beta_sesame]
    shared[, delta_beta_mean := (delta_beta_minfi + delta_beta_sesame) / 2]
    shared[, concordant_direction := direction_minfi == direction_sesame]
    setorder(shared, adj.P.Val_min, adj.P.Val_max, P.Value_min)
    metric_cols <- c("logFC_mean", "adj.P.Val_min", "adj.P.Val_max", "P.Value_min", "delta_beta_diff", "delta_beta_mean", "concordant_direction")
    cols_no_metrics <- setdiff(names(shared), metric_cols)
    rel_idx <- match("relation_to_island", cols_no_metrics)
    if (!is.na(rel_idx)) {
      metric_present <- metric_cols[metric_cols %in% names(shared)]
      post_rel <- if (rel_idx < length(cols_no_metrics)) cols_no_metrics[(rel_idx + 1L):length(cols_no_metrics)] else character()
      new_order <- c(
        cols_no_metrics[seq_len(rel_idx)],
        metric_present,
        post_rel
      )
      setcolorder(shared, new_order)
    }
  }
  shared
} else {
  data.table()
}
write_tsv_gz(annotated_shared, file.path(analysis_dir, "annotated_cpg_intersection_dual.tsv.gz"))

# ---- Figures (static) -------------------------------------------------------

log_message("Generating static plots...")
fig_dir <- paths$figures

if (!is.null(model_evaluations_ranked)) {
  eval_ok_plot <- model_evaluations_ranked[status == "ok" & !is.na(rank)]
  if (nrow(eval_ok_plot) > 0) {
    eval_ok_plot <- copy(eval_ok_plot)
    setorder(eval_ok_plot, rank, model_id)
    eval_ok_plot[, model_label := sprintf("#%s %s", rank, model_id)]
    metrics_long <- melt(
      eval_ok_plot,
      id.vars = c("model_id", "rank", "model_label", "use_sva", "use_combat", "batch_col"),
      measure.vars = c("n_sig_intersection", "n_sig_minfi", "n_sig_sesame"),
      variable.name = "metric",
      value.name = "value"
    )
    metrics_long[, metric := factor(
      metric,
      levels = c("n_sig_intersection", "n_sig_minfi", "n_sig_sesame"),
      labels = c("Shared CpGs", "minfi CpGs", "sesame CpGs")
    )]
    metrics_long[, model_label := factor(model_label, levels = unique(model_label))]
    model_metric_plot <- ggplot(metrics_long, aes(x = model_label, y = value, fill = metric)) +
      geom_col(position = position_dodge(width = 0.7), width = 0.65) +
      theme_minimal() +
      scale_fill_brewer(palette = "Set2", name = NULL) +
      labs(
        title = "Design Model Evaluation",
        subtitle = "Shared and pipeline-specific CpGs retained across candidate configurations",
        x = "Model (ranked)",
        y = "CpG count"
      ) +
      theme(
        axis.text.x = element_text(angle = 18, hjust = 1, vjust = 1),
        panel.grid.minor = element_blank()
      )
    write_plot(model_metric_plot, file.path(fig_dir, "model_selection_metrics"), width = 9, height = 6)
  }
}

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
  plot_pca(pca_minfi_pre$scores, metadata_dt, "dear_group", "PCA minfi (pre)"),
  plot_pca(pca_minfi_post$scores, metadata_dt, "dear_group", "PCA minfi (post)"),
  plot_pca(pca_sesame_pre$scores, metadata_dt, "dear_group", "PCA sesame (pre)"),
  plot_pca(pca_sesame_post$scores, metadata_dt, "dear_group", "PCA sesame (post)")
)
write_plot(pca_plots[[1]], file.path(fig_dir, "pca_minfi_pre"))
write_plot(pca_plots[[2]], file.path(fig_dir, "pca_minfi_post"))
write_plot(pca_plots[[3]], file.path(fig_dir, "pca_sesame_pre"))
write_plot(pca_plots[[4]], file.path(fig_dir, "pca_sesame_post"))

prepare_volcano_data <- function(results, fdr_threshold = opt$fdr_threshold, sample_size = 5000) {
  if (is.null(results)) {
    return(data.frame())
  }
  to_numeric <- function(x) {
    if (is.null(x)) return(rep(NA_real_, length(results[[1]] %||% numeric())))
    as.numeric(x)
  }
  delta_beta <- to_numeric(results[["delta_beta"]])
  adj_p <- to_numeric(results[["adj.P.Val"]])
  p_val <- results[["P.Value"]]
  if (is.null(p_val)) {
    p_val <- adj_p
  }
  p_val <- to_numeric(p_val)
  df <- data.frame(
    delta_beta = delta_beta,
    P.Value = p_val,
    adj.P.Val = adj_p,
    stringsAsFactors = FALSE
  )
  if ("probe_id" %in% names(results)) {
    df$probe_id <- as.character(results[["probe_id"]])
  }
  valid_mask <- is.finite(df$delta_beta) & is.finite(df$P.Value) & is.finite(df$adj.P.Val)
  df <- df[valid_mask, , drop = FALSE]
  if (nrow(df) == 0) {
    return(df)
  }
  df$P.Value <- pmax(df$P.Value, .Machine$double.xmin)
  df$adj.P.Val <- pmax(df$adj.P.Val, .Machine$double.xmin)
  if (nrow(df) > sample_size) {
    set.seed(1234)
    df <- df[sample.int(nrow(df), sample_size), , drop = FALSE]
  }
  df$neg_log_p <- -log10(df$P.Value)
  sig_flag <- if (!is.null(fdr_threshold) && length(fdr_threshold) == 1 && is.finite(fdr_threshold)) {
    df$adj.P.Val <= fdr_threshold
  } else {
    rep(FALSE, nrow(df))
  }
  df$significant <- factor(ifelse(sig_flag, "sig", "not_sig"), levels = c("not_sig", "sig"))
  df$hover_id <- if (!is.null(df$probe_id)) df$probe_id else sprintf("Row %s", seq_len(nrow(df)))
  df
}

volcano_plot <- function(results, title, fdr_threshold = opt$fdr_threshold) {
  log_message(
    "volcano_plot input for %s: rows=%s, class=%s",
    title,
    tryCatch(nrow(results), error = function(e) NA_integer_),
    paste(class(results), collapse = "/")
  )
  df <- prepare_volcano_data(results, fdr_threshold)
  log_message("volcano_plot df for %s: rows=%s, class=%s", title, nrow(df), paste(class(df), collapse = "/"))
  if (nrow(df) == 0) {
    return(ggplot() + theme_void() + labs(title = paste(title, "(insufficient data)")))
  }
  log_message("volcano_plot %s: downsampled to %s rows for plotting", title, nrow(df))
  ggplot(df, aes(x = delta_beta, y = neg_log_p, color = significant)) +
    geom_point(alpha = 0.6) +
    scale_color_manual(
      values = c(not_sig = "grey60", sig = "firebrick"),
      limits = c("not_sig", "sig"),
      breaks = c("not_sig", "sig"),
      drop = FALSE,
      na.translate = FALSE
    ) +
    theme_minimal() +
    labs(title = title, x = "Delta beta", y = "-log10(p-value)", color = paste0("FDR<=", signif(fdr_threshold)))
}

write_volcano_base <- function(results, title, filename, width = 7, height = 5) {
  df <- prepare_volcano_data(results, opt$fdr_threshold)
  if (nrow(df) == 0) {
    log_message("Fallback volcano plot for %s skipped: no finite data", title)
    return()
  }
  colors <- ifelse(df$significant == "sig", "firebrick", "grey60")
  png_filename <- paste0(filename, ".png")
  pdf_filename <- paste0(filename, ".pdf")
  png(png_filename, width = width, height = height, units = "in", res = 300)
  plot(
    df$delta_beta,
    df$neg_log_p,
    col = colors,
    pch = 16,
    cex = 0.6,
    main = title,
    xlab = "Delta beta",
    ylab = "-log10(p-value)"
  )
  legend("topright", legend = c(paste0("FDR <=", opt$fdr_threshold), paste0("FDR >", opt$fdr_threshold)), col = c("firebrick", "grey60"), pch = 16, bty = "n")
  dev.off()
  pdf(pdf_filename, width = width, height = height)
  plot(
    df$delta_beta,
    -log10(df$P.Value),
    col = colors,
    pch = 16,
    cex = 0.6,
    main = title,
    xlab = "Delta beta",
    ylab = "-log10(p-value)"
  )
  legend("topright", legend = c(paste0("FDR <=", opt$fdr_threshold), paste0("FDR >", opt$fdr_threshold)), col = c("firebrick", "grey60"), pch = 16, bty = "n")
  dev.off()
  log_message("Fallback volcano plot saved at %s (rows=%s)", filename, nrow(df))
}

build_volcano_plotly <- function(results, title, fdr_threshold = opt$fdr_threshold) {
  df <- prepare_volcano_data(results, fdr_threshold)
  if (nrow(df) == 0) {
    return(NULL)
  }
  colors <- c(not_sig = "grey60", sig = "firebrick")
  tooltip <- paste0(
    "Δβ: ", sprintf("%.3f", df$delta_beta),
    "<br>-log10(p-value): ", sprintf("%.2f", df$neg_log_p),
    "<br>FDR: ", sprintf("%.3g", df$adj.P.Val),
    if (!is.null(df$probe_id)) paste0("<br>CpG: ", df$probe_id) else ""
  )
  plotly::plot_ly(
    df,
    x = ~delta_beta,
    y = ~neg_log_p,
    color = ~significant,
    colors = colors,
    type = "scattergl",
    mode = "markers",
    marker = list(size = 6, opacity = 0.6),
    text = tooltip,
    hoverinfo = "text"
  ) %>%
    plotly::layout(
      title = title,
      xaxis = list(title = "Delta beta"),
      yaxis = list(title = "-log10(p-value)")
    )
}

write_volcano_plotly <- function(results, title, filename) {
  widget <- build_volcano_plotly(results, title, opt$fdr_threshold)
  if (is.null(widget)) {
    log_message("Fallback interactive volcano for %s skipped: no finite data", title)
    return(NULL)
  }
  html_path <- paste0(filename, ".html")
  tryCatch(
    {
      res <- save_styled_widget(
        widget,
        html_path,
        title = title,
        subtitle = "-log10(p-value) versus delta beta. Firebrick markers pass the FDR threshold.",
        description = "Fallback interactive view generated from summarised data."
      )
      log_message("Fallback interactive volcano saved at %s", res)
      res
    },
    error = function(e) {
      log_message("Fallback interactive volcano for %s failed: %s", title, conditionMessage(e))
      NULL
    }
  )
}

write_table_widget <- function(data, filename, title, subtitle = NULL, description = NULL) {
  html_path <- paste0(filename, ".html")
  if (nrow(data) == 0) {
    message_html <- sprintf(
      "<html><head><title>%s</title></head><body><h2>%s</h2><p>No records available.</p></body></html>",
      title,
      title
    )
    writeLines(message_html, con = html_path)
    log_message("Created empty table placeholder at %s", html_path)
    return(html_path)
  }
  write_fallback <- function(reason = NULL, persist_widget = FALSE) {
    if (!is.null(reason)) {
      log_message("Interactive table fallback for %s: %s", title, reason)
    }
    if (persist_widget && !is.null(widget)) {
      tryCatch(
        {
          saveRDS(widget, file = paste0(filename, ".rds"))
          log_message("Interactive table render failed; saved widget object to %s.rds for inspection", filename)
        },
        error = function(e) {
          log_message("Unable to persist failed widget for %s: %s", title, conditionMessage(e))
        }
      )
    }
    message_html <- sprintf(
      "<html><head><title>%s</title></head><body><h2>%s</h2><p>Interactive table could not be rendered. Check pipeline logs for details.</p></body></html>",
      title,
      title
    )
    writeLines(message_html, con = html_path)
    log_message("Fallback empty table written to %s", html_path)
    if (!is.null(ncol(data)) && ncol(data) > 0) {
      diagnostic_path <- paste0(filename, "_sample.tsv")
      tryCatch(
        {
          utils::write.table(
            head(data, 100),
            file = diagnostic_path,
            sep = "\t",
            row.names = FALSE,
            quote = FALSE
          )
          log_message("Saved preview of failed table to %s", diagnostic_path)
        },
        error = function(e) {
          log_message("Unable to write table preview for %s: %s", title, conditionMessage(e))
        }
      )
    }
    return(html_path)
  }
  if (nrow(data) > 2000) {
    log_message("Downsampling interactive table %s from %s to 2000 rows", title, nrow(data))
    data <- head(data, 2000)
  }
  use_buttons <- nrow(data) <= 1000
  widget_error <- NULL
  dt_args <- list(
    data = as.data.frame(data),
    rownames = FALSE,
    filter = "top",
    options = list(pageLength = 25, scrollX = TRUE)
  )
  if (use_buttons) {
    dt_args$extensions <- "Buttons"
  }
  widget <- tryCatch(
    do.call(datatable, dt_args),
    error = function(e) {
      widget_error <<- conditionMessage(e)
      NULL
    }
  )
  if (is.null(widget)) {
    reason <- widget_error %||% "datatable() construction failed with an unknown error."
    return(write_fallback(sprintf("Failed to build datatable widget: %s", reason), persist_widget = FALSE))
  }
  subtitle <- subtitle %||% "Search, filter, sort, and export the table using the controls above."
  description <- description %||% sprintf("Displaying %s rows.", format(nrow(data), big.mark = ","))
  res <- create_interactive(widget, filename, title = title, subtitle = subtitle, description = description)
  if (is.null(res)) {
    return(write_fallback("create_interactive returned NULL", persist_widget = TRUE))
  }
  res
}

volcano_minfi_classes <- paste(sprintf("%s(%s)", names(results_minfi), sapply(results_minfi, function(col) paste(class(col), collapse = "/"))), collapse = ", ")
log_message("volcano_minfi columns: %s", volcano_minfi_classes)
log_message("volcano_minfi rows: %s (finite adj.P.Val: %s)", nrow(results_minfi), sum(is.finite(results_minfi$adj.P.Val)))
if (!write_plot(volcano_plot(results_minfi, "Volcano minfi"), file.path(fig_dir, "volcano_minfi"))) {
  write_volcano_base(results_minfi, "Volcano minfi", file.path(fig_dir, "volcano_minfi"))
}
volcano_sesame_classes <- paste(sprintf("%s(%s)", names(results_sesame), sapply(results_sesame, function(col) paste(class(col), collapse = "/"))), collapse = ", ")
log_message("volcano_sesame columns: %s", volcano_sesame_classes)
log_message("volcano_sesame rows: %s (finite adj.P.Val: %s)", nrow(results_sesame), sum(is.finite(results_sesame$adj.P.Val)))
if (!write_plot(volcano_plot(results_sesame, "Volcano sesame"), file.path(fig_dir, "volcano_sesame"))) {
  write_volcano_base(results_sesame, "Volcano sesame", file.path(fig_dir, "volcano_sesame"))
}
write_plot(volcano_plot(integrated$intersection, "Volcano intersection"), file.path(fig_dir, "volcano_intersection"))

density_plot <- function(beta_matrix, title) {
  df <- melt(as.data.table(beta_matrix, keep.rownames = "probe_id"), id.vars = "probe_id", variable.name = "sample", value.name = "beta")
  df <- as.data.frame(df)
  df <- df[is.finite(df$beta), , drop = FALSE]
  if (nrow(df) == 0) {
    return(ggplot() + theme_void() + labs(title = paste(title, "(insufficient data)")))
  }
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
  df <- df[is.finite(df$beta), , drop = FALSE]
  if (nrow(df) == 0) {
    return(ggplot() + theme_void() + labs(title = paste(title, "(insufficient data)")))
  }
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
  has_gene <- "gene_symbols" %in% names(df)
  gene_vals <- if (has_gene) df$gene_symbols else rep(NA_character_, nrow(df))
  gene_str <- if (has_gene) ifelse(!is.na(gene_vals) & gene_vals != "", paste0("<br>Gene: ", gene_vals), "") else ""
  has_adj <- "adj.p.val" %in% names(df)
  adj_vals <- if (has_adj) df$adj.p.val else rep(NA_real_, nrow(df))
  adj_str <- if (has_adj) ifelse(is.na(adj_vals), "", sprintf("<br>adj.p-value: %.4g", adj_vals)) else ""
  has_delta <- "delta_beta" %in% names(df)
  delta_vals <- if (has_delta) df$delta_beta else rep(NA_real_, nrow(df))
  delta_str <- if (has_delta) ifelse(is.na(delta_vals), "", sprintf("<br>delta-beta: %.3f", delta_vals)) else ""
  has_logfc <- "logfc" %in% names(df)
  logfc_vals <- if (has_logfc) df$logfc else rep(NA_real_, nrow(df))
  logfc_str <- if (has_logfc) ifelse(is.na(logfc_vals), "", sprintf("<br>logFC: %.3f", logfc_vals)) else ""
  pval_str <- ifelse(is.na(df$p.value), "", sprintf("<br>p-value: %.4g", df$p.value))
  df[, tooltip := paste0(
    "CpG: ", probe_id,
    gene_str,
    pval_str,
    adj_str,
    delta_str,
    logfc_str
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

if (!is.null(model_evaluation_table_export) && nrow(model_evaluation_table_export) > 0) {
  eval_display <- copy(model_evaluation_table_export)
  setorder(eval_display, rank, model_id)
  display_cols <- c(
    "rank", "model_id", "status", "batch_col", "batch_in_design",
    "use_combat", "use_sva", "n_covariates_numeric", "n_covariates_factor",
    "n_surrogates", "n_sig_intersection", "n_sig_minfi", "n_sig_sesame",
    "batch_median_p_minfi", "batch_median_p_sesame", "group_median_p", "message"
  )
  existing_cols <- intersect(display_cols, names(eval_display))
  eval_display <- eval_display[, ..existing_cols]
  setnames(eval_display, existing_cols, c(
    "Rank", "Model", "Status", "Batch Column", "Batch in Design",
    "ComBat", "SVA", "Numeric Covars", "Factor Covars",
    "Surrogates", "Shared CpGs", "minfi CpGs", "sesame CpGs",
    "Batch p (minfi)", "Batch p (sesame)", "Group p (minfi)", "Notes"
  )[seq_along(existing_cols)])
  eval_display_df <- as.data.frame(eval_display)
  meta <- interactive_output_metadata("model_selection_table")
  widget <- datatable(
    eval_display_df,
    options = list(pageLength = 10, order = list(list(0, "asc"))),
    rownames = FALSE,
    filter = "top",
    caption = htmltools::tags$caption(
      style = "caption-side: top; text-align: left;",
      htmltools::strong("How the optimiser ranked each model")
    )
  )
  if ("Rank" %in% names(eval_display_df) && any(!is.na(eval_display_df$Rank))) {
    best_rank <- min(eval_display_df$Rank, na.rm = TRUE)
    if (is.finite(best_rank)) {
      widget <- formatStyle(widget, "Rank", target = "row", backgroundColor = DT::styleEqual(best_rank, "rgba(14,165,233,0.18)"))
    }
  }
  path <- create_interactive(
    widget,
    file.path(interactive_dir, "model_selection_table"),
    title = meta$title,
    subtitle = meta$subtitle,
    description = meta$description
  )
  if (!is.null(path)) interactive_files$model_selection_table <- path
}
meta <- interactive_output_metadata("pca_pre")
path <- create_interactive(
  ggplotly(pca_plots[[1]]),
  file.path(interactive_dir, "pca_pre"),
  title = meta$title,
  subtitle = meta$subtitle,
  description = meta$description
)
if (!is.null(path)) interactive_files$pca_pre <- path
meta <- interactive_output_metadata("pca_post")
path <- create_interactive(
  ggplotly(pca_plots[[2]]),
  file.path(interactive_dir, "pca_post"),
  title = meta$title,
  subtitle = meta$subtitle,
  description = meta$description
)
if (!is.null(path)) interactive_files$pca_post <- path
meta <- interactive_output_metadata("volcano_minfi")
volcano_minfi_widget <- build_volcano_plotly(annotated_minfi, "Volcano minfi")
if (!is.null(volcano_minfi_widget)) {
  path <- create_interactive(
    volcano_minfi_widget,
    file.path(interactive_dir, "volcano_minfi"),
    title = meta$title,
    subtitle = meta$subtitle,
    description = meta$description
  )
} else {
  path <- NULL
}
if (!is.null(path)) {
  interactive_files$volcano_minfi <- path
} else {
  fallback <- write_volcano_plotly(annotated_minfi, "Volcano minfi", file.path(interactive_dir, "volcano_minfi"))
  if (!is.null(fallback)) interactive_files$volcano_minfi <- fallback
}
meta <- interactive_output_metadata("volcano_intersection")
volcano_intersection_widget <- build_volcano_plotly(integrated$intersection, "Volcano intersection")
if (!is.null(volcano_intersection_widget)) {
  path <- create_interactive(
    volcano_intersection_widget,
    file.path(interactive_dir, "volcano_intersection"),
    title = meta$title,
    subtitle = meta$subtitle,
    description = meta$description
  )
  if (!is.null(path)) interactive_files$volcano_intersection <- path
}
meta <- interactive_output_metadata("manhattan_minfi")
path <- create_interactive(
  ggplotly(manhattan_plot(annotated_minfi, "Manhattan minfi")),
  file.path(interactive_dir, "manhattan_minfi"),
  title = meta$title,
  subtitle = meta$subtitle,
  description = meta$description
)
if (!is.null(path)) interactive_files$manhattan_minfi <- path
meta <- interactive_output_metadata("manhattan_intersection")
path <- create_interactive(
  ggplotly(manhattan_plot(annotated_intersection, "Manhattan intersection")),
  file.path(interactive_dir, "manhattan_intersection"),
  title = meta$title,
  subtitle = meta$subtitle,
  description = meta$description
)
if (!is.null(path)) interactive_files$manhattan_intersection <- path

meta <- interactive_output_metadata("table_minfi")
path <- write_table_widget(
  annotated_minfi,
  file.path(interactive_dir, "table_minfi"),
  title = meta$title,
  subtitle = meta$subtitle,
  description = meta$description
)
if (!is.null(path)) interactive_files$table_minfi <- path
if (sesame_available && nrow(results_sesame) > 0) {
  meta <- interactive_output_metadata("volcano_sesame")
  volcano_sesame_widget <- build_volcano_plotly(annotated_sesame, "Volcano sesame")
  if (!is.null(volcano_sesame_widget)) {
    path <- create_interactive(
      volcano_sesame_widget,
      file.path(interactive_dir, "volcano_sesame"),
      title = meta$title,
      subtitle = meta$subtitle,
      description = meta$description
    )
  } else {
    path <- NULL
  }
  if (!is.null(path)) {
    interactive_files$volcano_sesame <- path
  } else {
    fallback <- write_volcano_plotly(annotated_sesame, "Volcano sesame", file.path(interactive_dir, "volcano_sesame"))
    if (!is.null(fallback)) interactive_files$volcano_sesame <- fallback
  }
  meta <- interactive_output_metadata("manhattan_sesame")
  path <- create_interactive(
    ggplotly(manhattan_plot(annotated_sesame, "Manhattan sesame")),
    file.path(interactive_dir, "manhattan_sesame"),
    title = meta$title,
    subtitle = meta$subtitle,
    description = meta$description
  )
  if (!is.null(path)) interactive_files$manhattan_sesame <- path
  meta <- interactive_output_metadata("table_sesame")
  path <- write_table_widget(
    annotated_sesame,
    file.path(interactive_dir, "table_sesame"),
    title = meta$title,
    subtitle = meta$subtitle,
    description = meta$description
  )
  if (!is.null(path)) interactive_files$table_sesame <- path
}
meta <- interactive_output_metadata("table_intersection")
path <- write_table_widget(
  annotated_intersection,
  file.path(interactive_dir, "table_intersection"),
  title = meta$title,
  subtitle = meta$subtitle,
  description = meta$description
)
if (!is.null(path)) interactive_files$table_intersection <- path

meta <- interactive_output_metadata("table_intersection_dual")
path <- write_table_widget(
  annotated_shared,
  file.path(interactive_dir, "table_intersection_dual"),
  title = meta$title,
  subtitle = meta$subtitle,
  description = meta$description
)
if (!is.null(path)) interactive_files$table_intersection_dual <- path

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

detailed_overlap_stats <- list()
if (nrow(results_minfi) > 0 && nrow(results_sesame) > 0) {
  overlap_topn <- function(minfi_dt, sesame_dt, n = 1000) {
    n1 <- min(n, nrow(minfi_dt))
    n2 <- min(n, nrow(sesame_dt))
    top_minfi <- head(minfi_dt[order(P.Value)], n1)
    top_sesame <- head(sesame_dt[order(P.Value)], n2)
    if ("probe_id" %in% names(top_minfi) && "probe_id" %in% names(top_sesame)) {
      overlap <- intersect(top_minfi$probe_id, top_sesame$probe_id)
      count_minfi <- length(top_minfi$probe_id)
      count_sesame <- length(top_sesame$probe_id)
      list(
        top_n_minfi = count_minfi,
        top_n_sesame = count_sesame,
        overlap = length(overlap),
        overlap_pct_minfi = if (count_minfi > 0) length(overlap) / count_minfi else 0,
        overlap_pct_sesame = if (count_sesame > 0) length(overlap) / count_sesame else 0
      )
    } else {
      list(
        top_n_minfi = nrow(top_minfi),
        top_n_sesame = nrow(top_sesame),
        overlap = 0,
        overlap_pct_minfi = 0,
        overlap_pct_sesame = 0
      )
    }
  }
  detailed_overlap_stats <- overlap_topn(results_minfi, results_sesame, n = min(1000, opt$top_n_cpgs))
} else {
  detailed_overlap_stats <- list(
    top_n_minfi = nrow(results_minfi),
    top_n_sesame = nrow(results_sesame),
    overlap = 0,
    overlap_pct_minfi = 0,
    overlap_pct_sesame = 0
  )
}

dropped_covariates_summary <- if (length(dropped_covariates) > 0) {
  lapply(names(dropped_covariates), function(name) list(column = name, reason = dropped_covariates[[name]]))
} else {
  list()
}

metric_entry <- function(dt) {
  if (is.null(dt) || nrow(dt) == 0) {
    return(list())
  }
  as.list(dt[1])
}

batch_metrics <- list(
  minfi_pre = metric_entry(pca_minfi_pre$assessments[variable %in% candidate_batches, .(median_p = median(p_value, na.rm = TRUE), median_r2 = median(r_squared, na.rm = TRUE))]),
  minfi_post = metric_entry(pca_minfi_post$assessments[variable %in% candidate_batches, .(median_p = median(p_value, na.rm = TRUE), median_r2 = median(r_squared, na.rm = TRUE))]),
  sesame_pre = metric_entry(pca_sesame_pre$assessments[variable %in% candidate_batches, .(median_p = median(p_value, na.rm = TRUE), median_r2 = median(r_squared, na.rm = TRUE))]),
  sesame_post = metric_entry(pca_sesame_post$assessments[variable %in% candidate_batches, .(median_p = median(p_value, na.rm = TRUE), median_r2 = median(r_squared, na.rm = TRUE))])
)

surrogate_count <- if (!is.null(surrogate_vars)) ncol(as.matrix(surrogate_vars)) else 0

model_evaluations_dt <- if (is.null(model_evaluations)) data.table() else model_evaluations
model_evaluations_export <- if (nrow(model_evaluations_dt) > 0) {
  lapply(seq_len(nrow(model_evaluations_dt)), function(i) as.list(model_evaluations_dt[i]))
} else {
  list()
}
ranking_export <- if (!is.null(model_evaluations_ranked) && nrow(model_evaluations_ranked) > 0) {
  lapply(seq_len(nrow(model_evaluations_ranked)), function(i) as.list(model_evaluations_ranked[i]))
} else {
  list()
}

covariate_sets_export <- lapply(design_selection_summary$covariate_sets, function(set) {
  list(
    id = set$id,
    numeric = set$numeric,
    factor = set$factor,
    optional_numeric = set$optional_numeric,
    optional_factor = set$optional_factor
  )
})

covariate_stats_export <- list(
  optional_numeric = {
    dt <- design_selection_summary$covariate_stats$optional_numeric
    if (!is.null(dt) && nrow(dt) > 0) {
      lapply(seq_len(nrow(dt)), function(i) as.list(dt[i]))
    } else {
      list()
    }
  },
  optional_factor = {
    dt <- design_selection_summary$covariate_stats$optional_factor
    if (!is.null(dt) && nrow(dt) > 0) {
      lapply(seq_len(nrow(dt)), function(i) as.list(dt[i]))
    } else {
      list()
    }
  }
)

design_selection_export <- list(
  strategy = design_selection_summary$strategy,
  status = design_selection_summary$status,
  model_id = design_selection_summary$model_id,
  metrics = design_selection_summary$metrics,
  evaluated = model_evaluations_export,
  ranking = ranking_export,
  covariate_sets = covariate_sets_export,
  covariate_stats = covariate_stats_export,
  batch_options = design_selection_summary$batch_options
)

summary <- list(
  gse = opt$gse,
  samples = nrow(config),
  platform = platform_label,
  array_type = platform_info$dmr,
  sesame_manifest = sesame_manifest_source,
  fdr_threshold = opt$fdr_threshold,
  delta_beta_threshold = opt$delta_beta_threshold,
  groups = as.list(group_counts),
  candidate_batches = list(
    auto = batch_tracking$auto,
    manual = batch_tracking$manual,
    final = candidate_batches,
    selected = batch_to_use
  ),
  protected_batch_columns = batch_tracking$protected,
  covariates = covars,
  manual_covariates = manual_covariates,
  dropped_covariates = dropped_covariates_summary,
  batch_methods = list(minfi = selected_minfi_method, sesame = selected_sesame_method),
  combat_applied = isTRUE(best_params$use_combat),
  batch_column = batch_to_use,
  sva_surrogates = surrogate_count,
  design_selection = design_selection_export,
  top_n_cpgs = opt$top_n_cpgs,
  significant_cpgs = list(
    minfi = nrow(sig_minfi),
    sesame = nrow(sig_sesame),
    intersection = nrow(integrated$intersection)
  ),
  overlap_top_n = detailed_overlap_stats,
  top_cpgs = top_cpgs,
  batch_metrics = batch_metrics,
  shared_cpgs = list(
    count = nrow(annotated_shared),
    concordant = if ("concordant_direction" %in% names(annotated_shared)) sum(annotated_shared$concordant_direction, na.rm = TRUE) else 0
  ),
  directories = list(
    preprocess = paths$preprocess,
    analysis = paths$analysis,
    figures = paths$figures,
    interactive = paths$interactive,
    runtime = paths$runtime
  ),
  runtime_files = list(
    analysis_summary = file.path(paths$runtime, "analysis_summary.json"),
    pipeline_log = file.path(paths$runtime, "pipeline.log"),
    run_config = file.path(paths$runtime, "run_config.json")
  ),
  files = list(
    preprocess = manifest_entries(paths$preprocess),
    analysis = manifest_entries(paths$analysis),
    figures = manifest_entries(paths$figures),
    interactive = manifest_entries(paths$interactive),
    runtime = manifest_entries(paths$runtime)
  )
)
dashboard_paths <- write_dashboard_index(project_root, interactive_dir, interactive_files, summary)
summary$dashboard <- list(
  root = dashboard_paths$root,
  interactive = dashboard_paths$interactive
)
write_json(summary, file.path(paths$runtime, "analysis_summary.json"), pretty = TRUE, auto_unbox = TRUE)

log_message("DearMeta analysis completed successfully.")
