# DEARMETA

DEARMETA is a command-line toolkit that downloads, preprocesses, and analyses Illumina EPIC GEO methylation datasets. A single `dearmeta` command orchestrates metadata collection, IDAT retrieval, processing through an R pipeline, and report generation.

## Why DEARMETA?
- Automatically fetch GEO series metadata and platform annotations.
- Download paired IDAT files with caching to avoid duplicate transfers.
- Create a consistent workspace layout for every dataset.
- Run an opinionated R analysis pipeline (`scripts/analysis.R`) to produce QC plots, differential methylation tables, and HTML reports.
- Preserve run configurations (`configure.tsv` in each workspace root) and logs for reproducibility.

## Prerequisites

| Requirement | Details |
|-------------|---------|
| Python      | 3.10 or newer |
| R           | 4.0 or newer (needed for the bundled analysis script) |
| Git         | For cloning and keeping the project up to date |
| Internet    | Required for GEO API calls and supplementary file downloads |

**Recommended system packages (Linux):**
```bash
sudo apt-get update
sudo apt-get install -y build-essential libcurl4-openssl-dev libxml2-dev libssl-dev
```
These are commonly needed when compiling R or Python dependencies.

**Apple Silicon (macOS M1–M4) tips:**
- Install [Homebrew](https://brew.sh/) if it is not already available, then run:
  ```bash
  brew install python@3.11 libxml2 curl openssl
  ```
  If you prefer the Homebrew Python, prepend `/opt/homebrew/opt/python@3.11/bin` to your `PATH`.
- Create and activate a dedicated virtual environment *before* running any DearMeta scripts so `python3` resolves to that interpreter:
  ```bash
  python3 -m venv .dearmeta
  source .dearmeta/bin/activate
  ```
  The bootstrap script installs DearMeta in editable mode, so the `dearmeta` command becomes available inside this environment automatically.
- Install R 4.4.x (arm64) from CRAN. That release pairs with Bioconductor 3.20, which the installer script now auto-detects and uses to pull the right package versions.
- If R packages complain about headers, export `PKG_CONFIG_PATH="/opt/homebrew/opt/libxml2/lib/pkgconfig:$PKG_CONFIG_PATH"` before re-running `bash scripts/bootstrap.sh` (which calls `Rscript scripts/install.R` under the hood).

## Installation (Step by Step)

1. **Clone the repository**
   ```bash
   git clone https://github.com/kangk1204/dearmeta.git
   cd dearmeta
   ```

2. **Run the bootstrap script (installs Python + R deps)**
   ```bash
   bash scripts/bootstrap.sh
   ```
   This wraps the three manual steps (pip install, `renv::restore`, `scripts/install.R`). You only need to re-run it when dependencies change.

3. *(Optional)* **Create/activate your own virtual environment**
   ```bash
   python -m venv .venv
   source .venv/bin/activate
   ```

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

## Development Notes
- Core Python sources live in `src/dearmeta/`.
- R analysis scripts are stored under `scripts/`.
- Run `ruff` and `pytest` (installed via `pip install -e .[dev]`) before opening a pull request.

## License

Distributed under the Apache License 2.0. See `LICENSE` for full text.
