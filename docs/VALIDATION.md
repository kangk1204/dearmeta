# DearMeta validation checklist

This document outlines a lightweight but repeatable validation plan to keep DearMeta publication‑ready across diverse datasets.

## 1) Functional smoke tests (no data needed)
- **CLI option wiring**: `dearmeta --help` and `dearmeta analysis --help` to confirm options show `--intersection-choice`, `--dmr-intersection-top`, QC flags, etc.
- **Unit tests**: `pytest -q` (covers CLI forwarding for `--intersection-choice`, renv detection, GEO parsing).
- **Package presence**: `Rscript -e "sessionInfo()"` after `renv::restore()`; ensure `minfi`, `sesame`, `sesameData`, `sva`, `dmrff`, `VennDiagram` are loaded.

## 2) Controlled dataset panel
Run a panel of representative GEO workspaces (blood/non‑blood, small/medium cohorts, 450K/EPIC):
- Small blood (batch‑sparse): `GSE206418` (EPICv1, n≈20).
- Medium placenta (imbalanced batches): `GSE71678` variants (n=30, 50).
- Optional 450K check: any small 450K series with ≥2 groups.

Commands (per workspace):
```bash
dearmeta analysis --gse GSE71678 --group-ref control --intersection-choice worst
dearmeta analysis --gse GSE206418 --group-ref control --intersection-choice best
```
Environment knobs for speed/repro:
- `DEARMETA_MAX_MODELS=80` to cap model search if needed.
- `DEARMETA_PRESCREEN_PROBES=30000` to enable probe subsampling when matrices are huge.
- `DEARMETA_AUTO_INSTALL_PACKAGES=1` to auto-install FlowSorted.* when blood references are needed.

## 3) What to review after each run
- **Design selection**: `runtime/analysis_summary.json`
  - `design_selection.strategy` (`optimisation` vs `fallback`)
  - `batch_methods.minfi/sesame`, `combat_applied`, `sva_surrogates`
  - `integration.intersection_choice`
  - `batch_filtering.small_cohort_excluded` for excluded batch columns
- **Batch QC**: `05_interactive/batch_qc.html` and `analysis_summary.json -> batch_metrics` (median_p, median_r2 pre/post).
- **PCA plots**: `04_figures/pca_minfi_pre/post`, `pca_sesame_pre/post` for visual drift.
- **Signal**: `significant_cpgs` and `dmr_counts` in the summary; volcano/manhattan PDFs for minfi, sesame, intersection.
- **Logs**: `runtime/pipeline.log` for ComBat/SVA warnings, auto-install messages, and blacklisted batches.

## 4) Acceptance heuristics
- Batch metric median p ≥0.8 preferred; if unmet, ensure the selected model still maximizes batch_score and that ComBat is appropriately disabled/blacklisted.
- Small cohorts (n≤30): verify batch columns with sparse group×batch cells are excluded; SVA capped at ≤5 surrogates.
- Intersection choice: for conservative reporting use `worst`; for sensitivity use `best`; set explicitly in commands and confirm in `analysis_summary.json`.
- NA handling: covariates with <2 non‑missing or <2 unique values are dropped; batch candidates without cross‑group replication are excluded.

## 5) Performance watchpoints
- Slowest stages: sesame preprocessing (pOOBAH, manifest), dmrff DMR calling, large interactive tables/heatmaps.
- Mitigations: lower `--top-n-cpgs`, `--dmr-intersection-top`, or enable prescreen (`DEARMETA_PRESCREEN_PROBES`).

## 6) Reporting
For each dataset, archive:
- `analysis_summary.json` + `analysis_summary.md`
- Key figures (PCA pre/post, volcano/minfi/sesame/intersection, manhattan/intersection)
- `pipeline.log`

This checklist can be automated later (e.g., GitHub Actions with cached demo data), but it is immediately runnable on local workspaces.***
