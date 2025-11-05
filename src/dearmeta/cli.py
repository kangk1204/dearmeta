"""Command line interface for the DearMeta application."""

from __future__ import annotations

import json
import logging
import shutil
import sys
from contextlib import ExitStack
from pathlib import Path
from typing import Optional

import pandas as pd
import typer
from rich.console import Console

from .configure import (
    assemble_metadata_frame,
    export_metadata_tsv,
    generate_configure_tsv,
    harmonise_columns,
    infer_column_types,
    summarize_metadata,
)
from .data import analysis_script
from .files import ensure_dir, write_json, write_text
from .geo import GeoClient
from .logging_utils import ensure_log_file_handler, get_logger, silence_external_loggers
from .run_r import run_r_analysis

app = typer.Typer(add_completion=False, help="Download and analyse Illumina EPIC GEO datasets.")

console = Console()
logger = get_logger()
silence_external_loggers()


def validate_gse(gse: str) -> str:
    gse = gse.upper()
    if not gse.startswith("GSE") or not gse[3:].isdigit():
        raise typer.BadParameter("GSE accession must look like GSE123456.")
    return gse


@app.callback()
def main_callback(verbose: bool = typer.Option(False, "--verbose", "-v", help="Increase log verbosity.")) -> None:
    """Global CLI options."""
    if verbose:
        logger.setLevel(logging.DEBUG)


def build_structure(root: Path) -> dict:
    """Create the required directory structure and return important paths."""
    paths = {
        "download": ensure_dir(root / "01_download"),
        "idat": ensure_dir(root / "01_download" / "idat"),
        "metadata": ensure_dir(root / "01_download" / "metadata"),
        "logs": ensure_dir(root / "01_download" / "logs"),
        "preprocess": ensure_dir(root / "02_preprocess"),
        "analysis": ensure_dir(root / "03_analysis"),
        "figures": ensure_dir(root / "04_figures"),
        "interactive": ensure_dir(root / "05_interactive"),
        "runtime": ensure_dir(root / "runtime"),
    }
    return paths


def _format_mapping(mapping: dict[str, object] | list | tuple | None) -> str:
    if not mapping:
        return "none"
    if isinstance(mapping, dict):
        return "; ".join(f"{key}={value}" for key, value in mapping.items()) or "none"
    if isinstance(mapping, (list, tuple, set)):
        items = [str(item) for item in mapping if str(item)]
        return ", ".join(items) if items else "none"
    return str(mapping)


def generate_markdown_summary(summary_path: Path) -> Path:
    """Render a concise Markdown summary next to analysis_summary.json."""
    if not summary_path.exists():
        raise FileNotFoundError(f"Missing analysis summary at {summary_path}")

    summary = json.loads(summary_path.read_text())
    lines: list[str] = []

    gse = summary.get("gse", "unknown")
    lines.append(f"# DearMeta Summary · {gse}")
    lines.append("")

    samples = summary.get("samples")
    groups = summary.get("groups", {})
    lines.append("## Cohort")
    if samples is not None:
        lines.append(f"- Samples analysed: {samples}")
    if groups:
        lines.append(f"- Group counts: {_format_mapping(groups)}")

    sample_qc = summary.get("sample_qc", {})
    det = sample_qc.get("detection", {}) if isinstance(sample_qc, dict) else {}
    sesame_qc = sample_qc.get("sesame", {}) if isinstance(sample_qc, dict) else {}
    lines.append("- Detection QC threshold: {}".format(det.get("threshold", "n/a")))
    dropped_det = det.get("dropped") or []
    lines.append(f"- Samples dropped (detection): {_format_mapping(dropped_det)}")
    flagged_sesame = sesame_qc.get("flagged") or []
    dropped_sesame = sesame_qc.get("dropped") or []
    lines.append(f"- Sesame flagged (> {sesame_qc.get('threshold', 'n/a')}): {_format_mapping(flagged_sesame)}")
    if dropped_sesame:
        lines.append(f"- Sesame QC dropped: {_format_mapping(dropped_sesame)}")

    probe_filters = summary.get("probe_filters", {})
    if probe_filters:
        lines.append("")
        lines.append("## Probe Filtering")
        for key in ["total_probes_raw", "after_detection", "after_snp", "after_autosome", "final_probes"]:
            if key in probe_filters:
                lines.append(f"- {key.replace('_', ' ').title()}: {probe_filters[key]}")

    covariates = summary.get("covariates", {})
    manual_covariates = summary.get("manual_covariates", {})
    lines.append("")
    lines.append("## Design")
    lines.append(f"- Covariates (numeric): {_format_mapping(covariates.get('numeric'))}")
    lines.append(f"- Covariates (factor): {_format_mapping(covariates.get('factor'))}")
    if manual_covariates:
        lines.append(f"- Manual covariates: {_format_mapping(manual_covariates)}")
    batch_methods = summary.get("batch_methods", {})
    if batch_methods:
        lines.append(f"- Batch correction: {_format_mapping(batch_methods)}")
    combat = summary.get("combat_models", {})
    if combat:
        lines.append(
            f"- ComBat attempts: {combat.get('attempted', 0)} (success {combat.get('succeeded', 0)})"
        )
        failures = combat.get("failure_messages") or []
        if failures:
            lines.append(f"  - Failure notes: {_format_mapping(failures)}")

    candidate_batches = summary.get("candidate_batches", {})
    if candidate_batches:
        lines.append(f"- Candidate batch columns: {_format_mapping(candidate_batches.get('final'))}")
        diag = candidate_batches.get("diagnostics")
        if diag:
            lines.append(f"  - Batch diagnostics: {_format_mapping(diag)}")

    significant = summary.get("significant_cpgs", {})
    if significant:
        lines.append("")
        lines.append("## Differential CpGs")
        lines.append(f"- minfi: {significant.get('minfi', 'n/a')}")
        lines.append(f"- sesame: {significant.get('sesame', 'n/a')}")
        lines.append(f"- intersection: {significant.get('intersection', 'n/a')}")

    output_path = summary_path.with_suffix(".md")
    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return output_path


def load_configure_dataframe(path: Path) -> pd.DataFrame:
    """Load configure.tsv, tolerating comment-prefixed header rows."""
    df = pd.read_csv(path, sep="\t", comment="#", encoding="utf-8-sig")
    if "dear_group" in df.columns:
        return df
    header_row = None
    with path.open("r", encoding="utf-8-sig") as handle:
        for idx, line in enumerate(handle):
            stripped = line.strip()
            if not stripped:
                continue
            normalized = stripped.lstrip("\ufeff")
            clean = normalized.lstrip("\"").strip()
            first_field = clean.split("\t", 1)[0].strip().strip("\"'").lower()
            if first_field == "dear_group":
                header_row = idx
                break
            if not normalized.startswith("#"):
                preview = stripped[:40]
                raise typer.BadParameter(
                    "configure.tsv must contain only comment lines (starting with '#') "
                    f"before the dear_group header; offending line {idx + 1}: {preview}"
                )
    if header_row is None:
        raise typer.BadParameter("Unable to locate 'dear_group' header in configure.tsv")
    df = pd.read_csv(path, sep="\t", comment="#", header=header_row, encoding="utf-8-sig")
    if "dear_group" not in df.columns:
        raise typer.BadParameter("configure.tsv missing 'dear_group' column")
    return df


@app.command()
def download(
    gse: str = typer.Option(..., "--gse", help="GEO series accession (e.g. GSE123456)."),
    min_samples: int = typer.Option(6, help="Abort if fewer samples with IDAT pairs remain."),
    cache_dir: Path = typer.Option(Path(".dearmeta_cache"), help="Cache directory for GEO matrices."),
) -> None:
    """Download GEO metadata and raw IDAT files, producing configure.tsv."""
    gse = validate_gse(gse)
    root = Path.cwd() / gse
    cache_dir = cache_dir.expanduser()
    if not cache_dir.is_absolute():
        cache_dir = (Path.cwd() / cache_dir).resolve()
    ensure_dir(cache_dir)
    paths = build_structure(root)
    pipeline_log = paths["runtime"] / "pipeline.log"
    ensure_log_file_handler(logger, pipeline_log)
    logger.info("Starting download workflow for %s", gse)

    client = GeoClient()
    try:
        series = client.fetch_series(gse)
    except ValueError as exc:
        raise typer.BadParameter(str(exc)) from exc
    logger.info(
        "Series title: %s | platform: %s | organism: %s",
        series.title,
        series.platform_version,
        series.organism,
    )
    matrix_path = client.download_series_matrix(gse, cache_dir=cache_dir)
    # Copy series matrix to metadata root for convenience
    matrix_dest = paths["metadata"] / Path(matrix_path).name
    if matrix_dest != matrix_path:
        shutil.copy2(matrix_path, matrix_dest)

    # Download platform annotation tables
    annotation_paths = []
    for platform_id in series.platform_ids:
        try:
            annot_path = client.download_platform_table(platform_id, paths["metadata"])
            annotation_paths.append(str(annot_path))
        except Exception as exc:  # pragma: no cover - external fetch best effort
            logger.warning("Failed to download platform table %s: %s", platform_id, exc)
    if annotation_paths:
        logger.info("Downloaded platform annotations: %s", ", ".join(annotation_paths))

    samples = client.fetch_samples(gse, cache_dir=cache_dir)
    logger.info("Cached series matrix at %s", matrix_dest)
    if not samples:
        raise typer.BadParameter("Unable to parse any sample metadata from the series matrix.")
    idat_pairs = client.discover_idat_pairs(gse, samples=samples, cache_dir=cache_dir)
    if not idat_pairs:
        raise typer.BadParameter("No supplementary IDAT pairs were found for this series.")

    usable_samples = {gsm: samples[gsm] for gsm in idat_pairs.keys() if gsm in samples}
    dropped = sorted(set(samples.keys()) - set(usable_samples.keys()))
    if len(usable_samples) < min_samples:
        raise typer.BadParameter(
            f"Only {len(usable_samples)} samples with complete IDAT pairs detected (<{min_samples})."
        )

    logger.info("Detected %s samples with IDAT pairs (dropped %s)", len(usable_samples), len(dropped))
    if dropped:
        logger.warning("Excluded samples without complete IDAT pairs: %s", ", ".join(sorted(dropped)))

    downloaded_pairs = {}
    for gsm, pair in idat_pairs.items():
        if gsm not in usable_samples:
            continue
        logger.debug("Downloading IDAT pair for %s", gsm)
        downloaded_pairs[gsm] = client.download_idat_pair(pair, paths["idat"])

    metadata_frame = assemble_metadata_frame(usable_samples, downloaded_pairs, series.platform_version, root)
    metadata_frame = harmonise_columns(metadata_frame)
    metadata_tsv = paths["metadata"] / "metadata.tsv"
    export_metadata_tsv(metadata_frame, metadata_tsv)

    metadata_frame, numeric_covariates, candidate_batches = infer_column_types(metadata_frame)
    metadata_summary = summarize_metadata(metadata_frame)
    configure_path = root / "configure.tsv"
    generate_configure_tsv(metadata_frame, configure_path, candidate_batches, numeric_covariates)

    run_config = {
        "gse": gse,
        "title": series.title,
        "platform_version": series.platform_version,
        "organism": series.organism,
        "sample_counts": {
            "total": series.sample_count,
            "with_metadata": len(samples),
            "with_idat_pairs": len(downloaded_pairs),
        },
        "dropped_samples": dropped,
        "candidate_batches": candidate_batches,
        "numeric_covariates": numeric_covariates,
        "metadata_summary": metadata_summary,
        "series_matrix": str(matrix_dest),
        "platform_annotations": annotation_paths,
    }
    write_json(paths["runtime"] / "run_config.json", run_config)

    readme_path = root / "README_GSE.txt"
    readme_content = f"""DearMeta workspace for {gse}

Structure:
- configure.tsv         : fill dear_group column before running analysis
- 01_download/idat/     : raw IDAT files (auto gunzipped)
- 01_download/metadata/ : GEO metadata artifacts
- 04_figures/           : static plots
- 05_interactive/       : HTML interactive reports
- runtime/              : logs, run_config.json

Next steps:
1. Review configure.tsv
2. Fill the dear_group column with biological group labels
3. Run `dearmeta analysis --gse {gse}`
"""
    write_text(readme_path, readme_content)

    console.print(f"[bold green]Download complete[/] for {gse}. Configure file at {configure_path}.")
    logger.info("Metadata summary: %s", json.dumps(metadata_summary, indent=2))


@app.command()
def analysis(
    gse: str = typer.Option(..., "--gse", help="GEO series accession (e.g. GSE123456)."),
    min_group_size: int = typer.Option(2, help="Minimum samples per group required for analysis."),
    group_ref: Optional[str] = typer.Option(
        None,
        "--group-ref",
        "--group_ref",
        help="Group label to use as the reference (baseline) in limma contrasts.",
    ),
    r_script: Optional[Path] = typer.Option(
        None,
        help="Path to the analysis R script. Defaults to the packaged resource.",
    ),
    top_n_cpgs: int = typer.Option(10000, help="Number of CpGs to retain (sorted by p-value) for plots and tables."),
    drop_sesame_failed: bool = typer.Option(
        False, help="Drop samples whose sesame pOOBAH failure rate exceeds the threshold."
    ),
    poobah_threshold: float = typer.Option(
        0.05,
        help="Threshold for sesame pOOBAH failure rate; used for flagging/dropping samples.",
    ),
    cell_comp_reference: str = typer.Option(
        "auto",
        help="Cell composition reference ('auto', 'blood', or 'none'). Auto enables inference for blood datasets only.",
    ),
) -> None:
    """Run the R-based analysis pipeline using prepared configure.tsv."""
    gse = validate_gse(gse)
    if not (0 < poobah_threshold < 1):
        raise typer.BadParameter("--poobah-threshold must be between 0 and 1")
    valid_cell_refs = {"auto", "blood", "none"}
    normalized_cell_ref = cell_comp_reference.strip().lower()
    if normalized_cell_ref not in valid_cell_refs:
        raise typer.BadParameter("--cell-comp-reference must be one of: auto, blood, none")
    root = Path.cwd() / gse
    paths = build_structure(root)
    pipeline_log = paths["runtime"] / "pipeline.log"
    ensure_log_file_handler(logger, pipeline_log)

    configure_path = root / "configure.tsv"
    legacy_config_path = paths["runtime"] / "configure.tsv"
    if not configure_path.exists():
        if legacy_config_path.exists():
            logger.warning(
                "Using legacy configure.tsv at %s (preferred location: %s)",
                legacy_config_path,
                configure_path,
            )
            configure_path = legacy_config_path
        else:
            raise typer.BadParameter(f"Missing configure.tsv at {configure_path}")

    with ExitStack() as stack:
        if r_script is None:
            script_path = stack.enter_context(analysis_script())
        else:
            if not r_script.exists():
                raise typer.BadParameter(f"Analysis script not found at {r_script}")
            script_path = r_script.resolve()

        logger.info("Loading configure TSV from %s", configure_path)
        config_df = load_configure_dataframe(configure_path)
        config_df["dear_group"] = config_df["dear_group"].astype("string").str.strip()
        retained = config_df[config_df["dear_group"].notna() & (config_df["dear_group"] != "")].copy()

        if retained.empty:
            raise typer.BadParameter("No samples with dear_group assigned. Edit configure.tsv first.")
        group_sizes = retained["dear_group"].value_counts()
        if group_sizes.shape[0] < 2:
            raise typer.BadParameter("At least two groups required in dear_group column.")
        if (group_sizes < min_group_size).any():
            problematic = group_sizes[group_sizes < min_group_size].to_dict()
            raise typer.BadParameter(f"Groups with insufficient sample count: {problematic}")
        normalized_group_ref: Optional[str] = None
        if group_ref is not None:
            candidate = group_ref.strip()
            if not candidate:
                raise typer.BadParameter("--group-ref cannot be empty")
            if candidate not in group_sizes.index:
                raise typer.BadParameter(
                    f"--group-ref '{candidate}' is not present in dear_group. Available groups: {sorted(group_sizes.index.tolist())}"
                )
            normalized_group_ref = candidate

        logger.info("Running R analysis for %s with %s samples", gse, retained.shape[0])
        logger.debug("Using R script at %s", script_path)
        extra_args = [
            "--top-n-cpgs",
            str(top_n_cpgs),
            "--poobah-threshold",
            str(poobah_threshold),
            "--cell-comp-reference",
            normalized_cell_ref,
        ]
        if normalized_group_ref:
            extra_args.extend(["--group-ref", normalized_group_ref])
        if drop_sesame_failed:
            extra_args.append("--drop-sesame-failed")
        env = {"DEARMETA_GROUPS_JSON": json.dumps(group_sizes.to_dict())}
        if normalized_group_ref:
            env["DEARMETA_GROUP_REFERENCE"] = normalized_group_ref
        run_r_analysis(
            gse=gse,
            project_root=root,
            config_path=configure_path,
            output_root=root,
            r_script=script_path,
            extra_args=extra_args,
            env=env,
        )
    summary_json = paths["runtime"] / "analysis_summary.json"
    summary_md: Optional[Path] = None
    if summary_json.exists():
        try:
            summary_md = generate_markdown_summary(summary_json)
        except Exception as exc:  # pragma: no cover - best effort summary rendering
            logger.warning("Failed to render Markdown summary for %s: %s", gse, exc)
    else:
        logger.warning("analysis_summary.json missing at %s", summary_json)

    console.print(
        f"[bold green]Analysis complete[/] for {gse}. See runtime/analysis_summary.json for details."
    )
    if summary_md is not None:
        try:
            rel_path = summary_md.relative_to(Path.cwd())
        except ValueError:
            rel_path = summary_md
        console.print(f"[dim]Markdown summary written to {rel_path}[/]")


def main() -> None:
    try:
        app()
    except Exception as exc:  # pragma: no cover - CLI fallback
        logger.exception("Application error: %s", exc)
        raise


if __name__ == "__main__":  # pragma: no cover
    main()
