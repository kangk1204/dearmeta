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
| R           | 4.5 or newer (matches Bioconductor 3.22 in `renv.lock`) |
| Git         | For cloning/updates |
| Internet    | Needed for GEO downloads |

---

## Ubuntu / Windows (WSL) Installation

### 1. System packages
```bash
sudo apt-get update
sudo apt-get install -y build-essential libcurl4-openssl-dev libxml2-dev libssl-dev \
  libpng-dev libjpeg-dev libtiff-dev libhdf5-dev zlib1g-dev libbz2-dev liblzma-dev
```
These headers are required for CRAN/Bioconductor packages (`xml2`, `png`, `rhdf5`, `Rsamtools`, …). Run the same commands inside Ubuntu WSL if you are on Windows.

### 2. Clone the repository
```bash
git clone https://github.com/kangk1204/dearmeta.git
cd dearmeta
```

### 3. Create a Python environment (pick one)
- **Conda**
  ```bash
  conda create -n dearmeta python=3.11 r-base=4.5.1 -c conda-forge
  conda activate dearmeta
  conda install -y -c conda-forge zlib libpng xz
  python -m pip install --upgrade pip
  pip install -r requirements.lock
  ```
- **Virtualenv**
  ```bash
  python -m venv .dearmeta
  source .dearmeta/bin/activate
  python -m pip install --upgrade pip
  pip install -r requirements.lock
  ```
  > With a plain virtualenv you must have R ≥4.5 installed separately (CRAN Ubuntu repo on Linux, CRAN installer on macOS) and ensure `Rscript` points to that build.

### 4. Prime Bioconductor (optional but speeds up installs)
```bash
Rscript scripts/prime-bioconductor.R
```
This script installs `BiocManager`, pins Bioconductor 3.22, and preloads the heavy genomics stack so Bootstrap finishes faster. Inspect the script if you want to reproduce the commands manually.

### 5. Run the bootstrap script
```bash
bash scripts/bootstrap.sh
```
It performs `pip install -r requirements.lock`, `renv::restore()`, and `Rscript scripts/install.R` in one go. The script verifies that your `Rscript` is new enough; override it via `RSCRIPT=/path/to/Rscript` if you maintain multiple R builds.

### 6. Verify
```bash
dearmeta --help
Rscript -e "sessionInfo()"
```
You should see the `renv` library path in the R output.

---

## macOS (Apple Silicon M1–M4) Installation

### 1. Homebrew prerequisites
```bash
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
brew install python@3.11 libxml2 curl openssl
```
Add Homebrew Python to PATH if desired:
```bash
echo 'export PATH="/opt/homebrew/opt/python@3.11/bin:$PATH"' >> ~/.zprofile
source ~/.zprofile
```

### 2. Install R 4.5+
Download the Apple Silicon `.pkg` from CRAN and install it. R 4.5 is recommended; the bootstrap script can fall back to 4.4/Bioc 3.20 only if unavoidable.

### 3. Create a virtual environment
```bash
python3 -m venv .dearmeta
source .dearmeta/bin/activate
python -m pip install --upgrade pip
pip install -r requirements.lock
```

### 4. Export pkg-config hints
```bash
export PKG_CONFIG_PATH="/opt/homebrew/opt/libxml2/lib/pkgconfig:$PKG_CONFIG_PATH"
```
Add this line to `~/.zprofile` if you rebuild frequently.

### 5. Clone + bootstrap
```bash
git clone https://github.com/kangk1204/dearmeta.git
cd dearmeta
bash scripts/bootstrap.sh
```

### 6. Verify
Same commands as Ubuntu:
```bash
dearmeta --help
Rscript -e "sessionInfo()"
```
Again, confirm the `renv` library path in the R output.

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
