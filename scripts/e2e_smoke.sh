#!/usr/bin/env bash
# Lightweight, opt-in smoke test harness for DearMeta.
# Usage:
#   GSE_ID=GSE242427 ./scripts/e2e_smoke.sh
# This will run `dearmeta download` and `dearmeta analysis` in the current repo.
# No defaults are run to avoid accidental large downloads.

set -euo pipefail

if [[ -z "${GSE_ID:-}" ]]; then
  echo "Set GSE_ID to run the smoke test. Example:"
  echo "  GSE_ID=GSE242427 ./scripts/e2e_smoke.sh"
  exit 1
fi

echo "[info] Running smoke test for ${GSE_ID}"

if ! command -v dearmeta >/dev/null 2>&1; then
  echo "[error] dearmeta CLI not found on PATH" >&2
  exit 2
fi

dearmeta --help >/dev/null

dearmeta download --gse "${GSE_ID}"

if [[ ! -f "${GSE_ID}/configure.tsv" ]]; then
  echo "[error] Missing ${GSE_ID}/configure.tsv after download." >&2
  exit 3
fi

if ! grep -q "dear_group" "${GSE_ID}/configure.tsv"; then
  echo "[error] configure.tsv missing dear_group column; fill it before analysis." >&2
  exit 4
fi

echo "[info] Reminder: fill ${GSE_ID}/configure.tsv (dear_group) before analysis."
echo "[info] Proceeding with analysis assuming configure.tsv is already populated."

dearmeta analysis --gse "${GSE_ID}"

echo "[info] Smoke test completed for ${GSE_ID}"
