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
    configure_path = paths["runtime"] / "configure.tsv"
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
- runtime/configure.tsv : fill dear_group column before running analysis
- 01_download/idat/     : raw IDAT files (auto gunzipped)
- 01_download/metadata/ : GEO metadata artifacts
- 04_figures/           : static plots
- 05_interactive/       : HTML interactive reports

Next steps:
1. Review runtime/configure.tsv
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
    r_script: Optional[Path] = typer.Option(
        None,
        help="Path to the analysis R script. Defaults to the packaged resource.",
    ),
    top_n_cpgs: int = typer.Option(10000, help="Number of CpGs to retain (sorted by p-value) for plots and tables."),
) -> None:
    """Run the R-based analysis pipeline using prepared configure.tsv."""
    gse = validate_gse(gse)
    root = Path.cwd() / gse
    paths = build_structure(root)
    pipeline_log = paths["runtime"] / "pipeline.log"
    ensure_log_file_handler(logger, pipeline_log)

    configure_path = paths["runtime"] / "configure.tsv"
    if not configure_path.exists():
        raise typer.BadParameter(f"Missing configure.tsv at {configure_path}")

    with ExitStack() as stack:
        if r_script is None:
            script_path = stack.enter_context(analysis_script())
        else:
            if not r_script.exists():
                raise typer.BadParameter(f"Analysis script not found at {r_script}")
            script_path = r_script.resolve()

        logger.info("Loading configure TSV from %s", configure_path)
        config_df = pd.read_csv(configure_path, sep="\t", comment="#")
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

        logger.info("Running R analysis for %s with %s samples", gse, retained.shape[0])
        logger.debug("Using R script at %s", script_path)
        run_r_analysis(
            gse=gse,
            project_root=root,
            config_path=configure_path,
            output_root=root,
            r_script=script_path,
            extra_args=["--top-n-cpgs", str(top_n_cpgs)],
            env={"DEARMETA_GROUPS_JSON": json.dumps(group_sizes.to_dict())},
        )
    console.print(f"[bold green]Analysis complete[/] for {gse}. See runtime/analysis_summary.json for details.")


def main() -> None:
    try:
        app()
    except Exception as exc:  # pragma: no cover - CLI fallback
        logger.exception("Application error: %s", exc)
        raise


if __name__ == "__main__":  # pragma: no cover
    main()
