# DEARMETA

DEARMETA is a command-line toolkit that downloads, preprocesses, and analyses Illumina EPIC GEO methylation datasets. A single `dearmeta` command orchestrates metadata collection, IDAT retrieval, processing through an R pipeline, and report generation.

## Why DEARMETA?
- Automatically fetch GEO series metadata and platform annotations.
- Download paired IDAT files with caching to avoid duplicate transfers.
- Create a consistent workspace layout for every dataset.
- Run an opinionated R analysis pipeline (`src/dearmeta/data/analysis.R`) to produce QC plots, differential methylation tables, and HTML reports.
- Preserve run configurations (`configure.tsv` in each workspace root) and logs for reproducibility.

## Prerequisites

| Requirement | Details |
|-------------|---------|
| Python      | 3.10 or newer |
| R           | 4.5 or newer (renv.lock pins Bioconductor 3.22 which requires R ≥ 4.5) |
| Git         | For cloning and keeping the project up to date |
| Internet    | Required for GEO API calls and supplementary file downloads |

**Recommended system packages (Linux/WSL):**
```bash
sudo apt-get update
sudo apt-get install -y build-essential libcurl4-openssl-dev libxml2-dev libssl-dev \\
  libpng-dev libjpeg-dev libtiff-dev libhdf5-dev zlib1g-dev libbz2-dev liblzma-dev
```
These headers/libs are required when compiling CRAN and Bioconductor packages such as `xml2`, `png`, `rhdf5`, and `Rsamtools`.

**Apple Silicon (macOS M1–M4) tips:**
- Install [Homebrew](https://brew.sh/) if it is not already available, then run:
  ```bash
  brew install python@3.11 libxml2 curl openssl
  ```
  If you prefer the Homebrew Python, prepend `/opt/homebrew/opt/python@3.11/bin` to your `PATH` and create a virtual environment (`python3 -m venv .dearmeta && source .dearmeta/bin/activate`) before running any DearMeta scripts.
- Install R 4.5+ using the official Apple Silicon installer from CRAN (the `.pkg` under “macOS arm64”). If you are temporarily pinned to R 4.4/Bioconductor 3.20, the bootstrap script auto-detects that release and installs the matching package matrix, but 4.5+ remains the recommended baseline.
- If R packages complain about headers, export `PKG_CONFIG_PATH="/opt/homebrew/opt/libxml2/lib/pkgconfig:$PKG_CONFIG_PATH"` before re-running `Rscript scripts/install.R` or `bash scripts/bootstrap.sh`.

## Installation (Step by Step)

1. **Clone the repository**
   ```bash
   git clone https://github.com/kangk1204/dearmeta.git
   cd dearmeta
   ```

2. **Create and activate a runtime environment (pick one)**

   **Conda (recommended, keeps Python + R together):**
   ```bash
   conda create -n dearmeta python=3.11 r-base=4.5.1 -c conda-forge
   conda activate dearmeta
   conda install -y -c conda-forge zlib libpng xz
   python -m pip install --upgrade pip
   pip install -r requirements.lock
   ```

   **System Python virtualenv:**
   ```bash
   python -m venv .dearmeta
   source .dearmeta/bin/activate
   python -m pip install --upgrade pip
   pip install -r requirements.lock
   ```
   > When using a standalone virtualenv you must also install R ≥4.5 yourself (e.g., via the CRAN Ubuntu repo or the macOS installer) and ensure `Rscript` points at that build.

3. *(Optional but recommended)* **Prime Bioconductor 3.22 inside the active environment**
   ```bash
   Rscript scripts/prime-bioconductor.R
   ```
   The script installs `BiocManager` + supporting CRAN helpers, pins Bioconductor 3.22, pulls in the "core" group (S4Vectors, IRanges, DelayedArray, SummarizedExperiment, etc.), and then caches the genomics stack (`XVector`, `Biostrings`, `Rsamtools`, `GEOquery`, …) that DearMeta needs before running the main bootstrap. If you prefer to run the commands manually, inspect the script contents for the exact `BiocManager::install()` calls.

4. **Run the bootstrap script (installs pinned Python + R deps)**
   ```bash
   bash scripts/bootstrap.sh
   ```
   This wraps the three manual steps (pip install, `renv::restore`, `scripts/install.R`). You only need to re-run it when dependencies change. The script reads `renv.lock` and refuses to continue if the `Rscript` on your PATH is older than the pinned R version; set `RSCRIPT=/path/to/Rscript-4.5` when you manage multiple builds.

## Quick Start

1. Pick a GEO series ID (e.g. `GSE123456`).
2. Run the download workflow:
   ```bash
   dearmeta download --gse GSE123456
   ```
3. Open `GSE123456/configure.tsv` and fill the `dear_group` column with your biological groups.
   > `configure.tsv` now lives directly under the workspace root (e.g. `GSE123456/configure.tsv`). If you are reusing a workspace created by an older DearMeta version, you might still have `runtime/configure.tsv`; the CLI will warn and fall back to that file, but new runs always write next to the workspace root for easier editing.
4. Launch the R analysis:
   ```bash
   dearmeta analysis --gse GSE123456
   ```
   > Need a specific baseline? Append `--group-ref normoweight` (alias `--group_ref`) so limma treats that group as the reference when building contrasts. Otherwise the first `dear_group` in your configure file becomes the default reference.
5. Reproduce this exact software stack:
   ```bash
   python -m pip install -r requirements.lock
   Rscript -e 'renv::restore(prompt = FALSE)'
   ```
   The `requirements.lock` and `renv.lock` files pin every Python and R dependency (including Bioconductor releases) so collaborators can recreate the DearMeta environment bit-for-bit.

### QC controls
- `--poobah-threshold 0.05` adjusts the sesame pOOBAH failure cutoff (between 0 and 1).
- `--drop-sesame-failed` removes samples whose sesame failure rate exceeds the threshold before continuing.
- Batch/sample processing currently runs serially for stability; there is no parallel worker flag to configure.
- See `docs/statistical_assumptions.md` for citations supporting every QC/batch threshold shipped with DearMeta.

### Cell composition (blood datasets)
- DearMeta can now estimate leukocyte fractions via the Houseman/IDOL method (`minfi::estimateCellCounts2`).
- Enabled by default (`--cell-comp-reference auto`): it auto-detects blood-like metadata and appends `cell_*` columns to `configure.tsv`, which are then available as numeric covariates.
- Force or disable behaviour with `--cell-comp-reference blood` or `--cell-comp-reference none` if auto-detection isn’t appropriate.

DearMeta will create a workspace per series:

```
GSE123456/
├── configure.tsv       # group assignments + covariates
├── 01_download/        # raw IDAT files, series matrix, metadata
├── 02_preprocess/      # preprocessed objects
├── 03_analysis/        # analysis outputs
├── 04_figures/         # static plots
├── 05_interactive/     # HTML/interactive reports
└── runtime/            # logs, run_config.json, pipeline artifacts
```

### Choosing a reference group (`--group-ref`)
- The `dear_group` column in `configure.tsv` tells DearMeta how samples cluster for differential testing; limma builds contrasts by comparing every group to a *reference* (baseline) group.
- If you do **not** pass `--group-ref`, the first non-empty `dear_group` value encountered in `configure.tsv` becomes the baseline. This might be arbitrary (e.g., whichever row you listed first).
- Pin the baseline explicitly with `dearmeta analysis --gse GSE123456 --group-ref control`. The value must match one of the strings already present in `dear_group`, and the CLI will warn if it does not.
- All volcano plots, top-table exports, and HTML dashboards will then label comparisons as `target_vs_control`, making it clear which direction the log-fold changes represent.

## Tips & Troubleshooting
- **Reusing downloads:** Copy `.dearmeta_cache/` from an existing machine to avoid re-downloading long-lived GEO artifacts.
- **Verbose logging:** Append `--verbose` to `dearmeta` commands to print debug logs and capture them in `runtime/pipeline.log`.
- **R package errors:** Ensure the system libraries listed above are installed, then rerun `Rscript scripts/install.R`.
- **R version mismatch:** `scripts/bootstrap.sh` requires R ≥ the version recorded in `renv.lock` (currently 4.5.1). Install a recent R build or point the script at the right binary via `RSCRIPT=/opt/R/4.5/bin/Rscript bash scripts/bootstrap.sh`.
- **Bioconductor bootstrap:** If `install.R` complains about missing packages such as `XVector` or `SparseArray`, rerun `Rscript scripts/prime-bioconductor.R` to preinstall the base Bioconductor stack before invoking the bootstrap script again.
- **libpng / zlib link errors:** When CRAN's `png` package fails with `cannot find -lz`, install the matching Conda libraries inside the active env (`conda install -c conda-forge zlib libpng xz`) so `R CMD INSTALL` can link against them.
- **libxml2 / XML package errors:** `scripts/bootstrap.sh` automatically points `XML_CONFIG` at `/usr/bin/xml2-config` and exports the required include/library paths when it exists. If you are on a non-Debian system, set `XML_CONFIG` to your system `xml2-config` path before running the script so CRAN's `XML`/`xml2` packages can link successfully.
- **Bioconductor version drift:** When Bioconductor only ships a newer release of a pinned package (e.g., `sesameData` or `FlowSorted.*`), the bootstrap script logs a warning and keeps the newer build because it is ≥ the lockfile version.

## Development Notes
- Core Python sources live in `src/dearmeta/`.
- R analysis scripts are stored under `scripts/`.
- Run `ruff` and `pytest` (installed via `pip install -e .[dev]`) before opening a pull request.

## License

Distributed under the Apache License 2.0. See `LICENSE` for full text.
