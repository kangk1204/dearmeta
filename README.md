# DearMeta

DearMeta is a command-line toolkit that downloads, preprocesses, and analyses Illumina EPIC GEO methylation datasets. A single `dearmeta` command orchestrates metadata collection, IDAT retrieval, processing through an R pipeline, and report generation.

## Why DearMeta?
- Automatically fetch GEO series metadata and platform annotations.
- Download paired IDAT files with caching to avoid duplicate transfers.
- Create a consistent workspace layout for every dataset.
- Run an opinionated R analysis pipeline (`scripts/analysis.R`) to produce QC plots, differential methylation tables, and HTML reports.
- Preserve run configurations (`runtime/configure.tsv`) and logs for reproducibility.

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
  python3 -m venv .venv
  ```
- Install R 4.0+ using the official Apple Silicon installer from CRAN (the `.pkg` under “macOS arm64”).
- If R packages complain about headers, export `PKG_CONFIG_PATH="/opt/homebrew/opt/libxml2/lib/pkgconfig:$PKG_CONFIG_PATH"` before re-running `Rscript scripts/install.R`.

## Installation (Step by Step)

1. **Clone the repository**
   ```bash
   git clone https://github.com/kangk1204/dearmeta.git
   cd dearmeta
   ```

2. **Create and activate a virtual environment**
   ```bash
   python -m venv .venv
   # Linux / macOS
   source .venv/bin/activate
   # Windows PowerShell
   # .\.venv\Scripts\Activate.ps1
   ```

3. **Install Python dependencies**
   ```bash
   pip install -U pip
   pip install -e .
   # or include developer tools:
   # pip install -e .[dev]
   ```

4. **Install R packages used by the analysis pipeline**
   ```bash
   Rscript scripts/install.R
   ```
   This script installs the required CRAN and Bioconductor packages. Run it once per machine (rerun whenever R is reinstalled).

## Quick Start

1. Pick a GEO series ID (e.g. `GSE123456`).
2. Run the download workflow:
   ```bash
   dearmeta download --gse GSE123456
   ```
3. Open `GSE123456/runtime/configure.tsv` and fill the `dear_group` column with your biological groups.
4. Launch the R analysis:
   ```bash
   dearmeta analysis --gse GSE123456
   ```

DearMeta will create a workspace per series:

```
GSE123456/
├── 01_download/        # raw IDAT files, series matrix, metadata
├── 02_preprocess/      # preprocessed objects
├── 03_analysis/        # analysis outputs
├── 04_figures/         # static plots
├── 05_interactive/     # HTML/interactive reports
└── runtime/            # configure.tsv, logs, run_config.json
```

Large GEO artifacts (`GSE*/`, `.dearmeta*/`, HTML reports, etc.) are excluded by `.gitignore`. You can safely delete them and regenerate with the CLI whenever needed.

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
