#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(cd -- "${SCRIPT_DIR}/.." && pwd)
PYTHON_BIN=${PYTHON:-python3}
RSCRIPT_BIN=${RSCRIPT:-Rscript}

if [[ -z "${XML_CONFIG:-}" ]]; then
  SYS_XML_CONFIG=/usr/bin/xml2-config
  if [[ -x "${SYS_XML_CONFIG}" ]]; then
    ACTIVE_XML_CONFIG=$(command -v xml2-config 2>/dev/null || true)
    if [[ -z "${ACTIVE_XML_CONFIG}" || "${ACTIVE_XML_CONFIG}" != "${SYS_XML_CONFIG}" ]]; then
      echo "[Bootstrap] Using system xml2-config at ${SYS_XML_CONFIG}" >&2
    fi
    export XML_CONFIG="${SYS_XML_CONFIG}"
    if [[ -d /usr/lib/x86_64-linux-gnu/pkgconfig ]]; then
      export PKG_CONFIG_PATH="/usr/lib/x86_64-linux-gnu/pkgconfig:${PKG_CONFIG_PATH:-}"
    fi
    if [[ -d /usr/lib/x86_64-linux-gnu ]]; then
      export LD_LIBRARY_PATH="/usr/lib/x86_64-linux-gnu:${LD_LIBRARY_PATH:-}"
      export LIBRARY_PATH="/usr/lib/x86_64-linux-gnu:${LIBRARY_PATH:-}"
      export LIBXML_LIBS="-L/usr/lib/x86_64-linux-gnu -Wl,-rpath,/usr/lib/x86_64-linux-gnu -lxml2 -lz -llzma"
    fi
    if [[ -d /usr/include/libxml2 ]]; then
      export LIBXML_INCDIR="/usr/include/libxml2"
    fi
  fi
fi

LOCK_META=$(${PYTHON_BIN} - <<'PY' "${REPO_ROOT}/renv.lock"
import json
import sys
from pathlib import Path

lock = {}
try:
    with Path(sys.argv[1]).open() as handle:
        lock = json.load(handle)
except Exception:
    pass

print(lock.get("R", {}).get("Version", ""))
print(lock.get("Packages", {}).get("renv", {}).get("Version", ""))
PY
)
REQUIRED_R_VERSION=$(printf '%s\n' "${LOCK_META}" | sed -n '1p')
RENV_VERSION=$(printf '%s\n' "${LOCK_META}" | sed -n '2p')

if [[ -n "${REQUIRED_R_VERSION}" ]]; then
  echo "[Bootstrap] Checking local R version (need >= ${REQUIRED_R_VERSION})" >&2
  CURRENT_R_VERSION=$(${RSCRIPT_BIN} --vanilla -e "cat(as.character(getRversion()))" 2>/dev/null || true)
  if [[ -z "${CURRENT_R_VERSION}" ]]; then
    echo "[Bootstrap] Unable to determine R version via ${RSCRIPT_BIN}. Install R ${REQUIRED_R_VERSION}+ and/or set RSCRIPT=/path/to/Rscript." >&2
    exit 1
  fi
  ${PYTHON_BIN} - <<'PY' "${CURRENT_R_VERSION}" "${REQUIRED_R_VERSION}"
import sys

def normalize(ver: str):
    parts = [int(p) for p in ver.split('.') if p]
    while len(parts) < 3:
        parts.append(0)
    return parts[:3]

current, required = sys.argv[1:3]
if normalize(current) < normalize(required):
    sys.stderr.write(
        "[Bootstrap] Rscript reports version {} but renv.lock requires >= {}.\n"
        "Install a newer R (>= {}) or rerun with RSCRIPT pointing to that version.\n".format(
            current, required, required
        )
    )
    sys.exit(1)
PY
fi

echo "[Bootstrap] Installing Python dependencies via ${PYTHON_BIN}" >&2
${PYTHON_BIN} -m pip install --upgrade pip
${PYTHON_BIN} -m pip install -r "${REPO_ROOT}/requirements.lock"
echo "[Bootstrap] Installing dearmeta package (editable) so the CLI is on PATH" >&2
${PYTHON_BIN} -m pip install -e "${REPO_ROOT}"

echo "[Bootstrap] Ensuring renv is installed" >&2
if [[ -n "${RENV_VERSION}" ]]; then
  ${RSCRIPT_BIN} --vanilla "${REPO_ROOT}/scripts/install-renv.R" "${RENV_VERSION}"
else
  ${RSCRIPT_BIN} --vanilla "${REPO_ROOT}/scripts/install-renv.R"
fi

echo "[Bootstrap] Restoring renv lockfile via ${RSCRIPT_BIN}" >&2
${RSCRIPT_BIN} --vanilla -e "renv::restore(prompt = FALSE)"
echo "[Bootstrap] Installing required R/Bioconductor packages" >&2
${RSCRIPT_BIN} "${REPO_ROOT}/scripts/install.R"

echo "[Bootstrap] Done." >&2
