"""Utilities for transforming GEO sample metadata into configure.tsv artifacts."""

from __future__ import annotations

import re
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import pandas as pd

from .files import ensure_dir, write_text
from .geo import IdatPair, SampleMetadata
from .logging_utils import get_logger

logger = get_logger(__name__)

CORE_COLUMNS = {
    "gsm_id",
    "sample_name",
    "idat_red",
    "idat_grn",
    "platform_version",
    "species",
}

CANDIDATE_BATCH_KEYS = {
    "slide",
    "barcode",
    "sentrix_id",
    "sentrix_position",
    "array",
    "plate",
    "batch",
    "processing_date",
    "center",
    "shipment_batch",
    "run_date",
    "chip",
}

MIN_BATCH_NONMISSING_RATIO = 0.5
MIN_BATCH_LEVEL_SIZE = 2
PROTECTED_BATCH_KEYWORDS = ("sex", "gender")
SENTRIX_PATTERN = re.compile(r"_(\d{8,12})_(R\d{2}C\d{2})", re.IGNORECASE)


def _normalize_name(name: str) -> str:
    return "".join(ch for ch in name.lower() if ch.isalnum())


def _is_protected_column(name: str) -> bool:
    normalized = _normalize_name(name)
    return any(keyword in normalized for keyword in PROTECTED_BATCH_KEYWORDS)


def _has_sufficient_batch_support(series: pd.Series, total_len: int) -> bool:
    series = series.dropna()
    if total_len <= 0 or series.empty:
        return False
    if len(series) < max(2 * MIN_BATCH_LEVEL_SIZE, 2):
        return False
    if len(series) / total_len < MIN_BATCH_NONMISSING_RATIO:
        return False
    counts = series.value_counts()
    return len(counts) >= 2 and (counts >= MIN_BATCH_LEVEL_SIZE).all()


def assemble_metadata_frame(
    samples: Dict[str, SampleMetadata],
    idat_pairs: Dict[str, IdatPair],
    platform_version: str,
    workspace_root: Path,
) -> pd.DataFrame:
    """Combine sample metadata with IDAT paths into a tidy DataFrame."""
    rows: List[dict] = []
    for gsm, sample in samples.items():
        pair = idat_pairs.get(gsm)
        if not pair:
            logger.debug("Skipping %s: no IDAT pair detected", gsm)
            continue
        row = {
            "gsm_id": gsm,
            "sample_name": sample.title or gsm,
            "platform_version": platform_version,
            "species": (sample.organism or "").strip() or "Homo sapiens",
            "idat_red": _relativise_path(pair.red, workspace_root),
            "idat_grn": _relativise_path(pair.green, workspace_root),
        }
        sentrix_id, sentrix_pos = _extract_sentrix_fields(pair.red or pair.green)
        if sentrix_id:
            row["sentrix_id"] = sentrix_id
        if sentrix_pos:
            row["sentrix_position"] = sentrix_pos
        for key, value in sample.characteristics.items():
            if value:
                row[key] = value
        rows.append(row)
    if not rows:
        raise ValueError("No samples with IDAT pairs available to build metadata frame.")
    frame = pd.DataFrame(rows)
    frame.sort_values("gsm_id", inplace=True)
    frame.reset_index(drop=True, inplace=True)
    return frame


def harmonise_columns(frame: pd.DataFrame) -> pd.DataFrame:
    """Standardise column names (lower snake_case) and trim whitespace."""
    def _normalise(name: str) -> str:
        name = name.strip().lower().replace(" ", "_")
        name = name.replace("-", "_")
        return name

    frame = frame.copy()
    frame.columns = [_normalise(col) for col in frame.columns]
    for column in frame.columns:
        if frame[column].dtype == object:
            frame[column] = frame[column].astype(str).str.strip()
            frame[column] = frame[column].replace({"nan": "", "none": ""})
    return frame


def infer_column_types(frame: pd.DataFrame) -> Tuple[pd.DataFrame, List[str], List[str]]:
    """Infer numeric vs categorical columns and identify candidate batch covariates."""
    frame = frame.copy()
    numeric_cols: List[str] = []
    categorical_cols: List[str] = []
    candidate_batches: List[str] = []

    for column in frame.columns:
        if column in CORE_COLUMNS:
            continue
        series = frame[column].fillna("")
        if not series.any():
            frame.drop(columns=[column], inplace=True)
            continue
        series_clean = series.replace("", pd.NA)
        non_missing = int(series_clean.notna().sum())
        unique_ratio = series_clean.nunique(dropna=True) / max(non_missing, 1)
        numeric_series = pd.to_numeric(series, errors="coerce")
        if numeric_series.notna().sum() >= len(series) * 0.6:
            frame[column] = numeric_series
            numeric_cols.append(column)
        else:
            frame[column] = series
            categorical_cols.append(column)
        if _is_protected_column(column):
            continue
        supports_batch = _has_sufficient_batch_support(series_clean.dropna(), len(series))
        if not supports_batch:
            continue
        keyword_match = any(key in column for key in CANDIDATE_BATCH_KEYS)
        if keyword_match and (unique_ratio < 0.9 or len(series) <= 12):
            candidate_batches.append(column)
        elif 0 < unique_ratio < 0.2 and series_clean.nunique(dropna=True) > 1:
            candidate_batches.append(column)

    # Remove redundant columns (duplicates)
    duplicates = _find_duplicate_columns(frame)
    if duplicates:
        frame.drop(columns=duplicates, inplace=True)
        numeric_cols = [col for col in numeric_cols if col not in duplicates]
        categorical_cols = [col for col in categorical_cols if col not in duplicates]
        candidate_batches = [col for col in candidate_batches if col not in duplicates]

    return frame, numeric_cols, sorted(set(candidate_batches))


def _find_duplicate_columns(frame: pd.DataFrame) -> List[str]:
    duplicates: List[str] = []
    seen = {}
    for column in frame.columns:
        signature = tuple(frame[column].fillna("").tolist())
        if signature in seen:
            previous = seen[signature]
            if column in CORE_COLUMNS or previous in CORE_COLUMNS:
                continue
            duplicates.append(column)
        else:
            seen[signature] = column
    return duplicates


def _relativise_path(location: str, workspace_root: Path) -> str:
    path = Path(location).resolve()
    try:
        return str(path.relative_to(workspace_root.resolve()))
    except ValueError:
        return str(path)


def _extract_sentrix_fields(path: str) -> tuple[Optional[str], Optional[str]]:
    name = Path(path).name
    match = SENTRIX_PATTERN.search(name)
    if not match:
        return None, None
    sentrix_id, position = match.groups()
    return sentrix_id, position.upper()


def generate_configure_tsv(
    frame: pd.DataFrame,
    path: Path,
    candidate_batches: Iterable[str],
    numeric_covariates: Iterable[str],
) -> None:
    """Write configure TSV with dear_group column first and header guidance."""
    ensure_dir(path.parent)
    candidate_batches = sorted(set(candidate_batches))
    numeric_covariates = sorted(set(numeric_covariates))

    guidance_lines = [
        "# dear_group: fill with biological group labels (e.g., case/control).",
        "# candidate_batch_columns: " + ",".join(candidate_batches) if candidate_batches else "# candidate_batch_columns: none detected",
        "# numeric_covariates: " + ",".join(numeric_covariates) if numeric_covariates else "# numeric_covariates: none detected",
        "# Only samples with complete dear_group values will be analysed.",
    ]

    frame = frame.copy()
    frame.insert(0, "dear_group", "")

    content = "\n".join(guidance_lines) + "\n"
    content += frame.to_csv(sep="\t", index=False, lineterminator="\n")
    write_text(path, content)


def export_metadata_tsv(frame: pd.DataFrame, path: Path) -> None:
    """Persist the raw metadata table for reference."""
    ensure_dir(path.parent)
    frame.to_csv(path, sep="\t", index=False)


def summarize_metadata(frame: pd.DataFrame) -> dict:
    """Return a compact JSON-serialisable summary of metadata coverage."""
    summary = {
        "n_samples": int(frame.shape[0]),
        "columns": {},
    }
    for column in frame.columns:
        series = frame[column]
        summary["columns"][column] = {
            "n_unique": int(series.nunique(dropna=True)),
            "pct_missing": float(series.isna().mean() * 100),
        }
    return summary
