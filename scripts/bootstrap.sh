#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(cd -- "${SCRIPT_DIR}/.." && pwd)
PYTHON_BIN=${PYTHON:-python3}
RSCRIPT_BIN=${RSCRIPT:-Rscript}

echo "[Bootstrap] Installing Python dependencies via ${PYTHON_BIN}" >&2
${PYTHON_BIN} -m pip install --upgrade pip
${PYTHON_BIN} -m pip install -r "${REPO_ROOT}/requirements.lock"
echo "[Bootstrap] Installing dearmeta package (editable) so the CLI is on PATH" >&2
${PYTHON_BIN} -m pip install -e "${REPO_ROOT}"

echo "[Bootstrap] Ensuring the renv package is available in R" >&2
${RSCRIPT_BIN} --vanilla -e "if (!requireNamespace('renv', quietly = TRUE)) install.packages('renv', repos = 'https://cloud.r-project.org')"
echo "[Bootstrap] Restoring renv lockfile via ${RSCRIPT_BIN}" >&2
${RSCRIPT_BIN} --vanilla -e "renv::restore(prompt = FALSE)"
echo "[Bootstrap] Installing required R/Bioconductor packages" >&2
${RSCRIPT_BIN} "${REPO_ROOT}/scripts/install.R"

echo "[Bootstrap] Done." >&2
