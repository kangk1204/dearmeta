"""Helpers for invoking the R analysis pipeline."""

from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path
from typing import Dict, Iterable, Optional

from .logging_utils import get_logger

logger = get_logger(__name__)


class RRuntimeError(RuntimeError):
    """Raised when the R pipeline exits with a non-zero status."""


def _has_bitmap_capabilities(rscript: str) -> bool:
    """Return True if the given Rscript reports a bitmap device (cairo/png/jpeg)."""
    try:
        out = subprocess.check_output(
            [
                rscript,
                "-e",
                "caps<-capabilities(); cat(any(caps[c('cairo','png','jpeg')]))",
            ],
            text=True,
        )
        return out.strip().lower() in {"true", "t", "1"}
    except Exception as exc:  # pragma: no cover - defensive
        logger.debug("Failed to query capabilities from %s: %s", rscript, exc)
        return False


def _resolve_rscript_candidates() -> str:
    """Pick an Rscript, preferring one with bitmap support for interactive plots."""
    from shutil import which

    env_rscript = os.environ.get("DEARMETA_RSCRIPT") or os.environ.get("RSCRIPT")
    candidates = [env_rscript] if env_rscript else []
    candidates.extend(["Rscript", "/usr/bin/Rscript"])

    seen = set()
    resolved_candidates = []
    for candidate in candidates:
        if not candidate or candidate in seen:
            continue
        seen.add(candidate)
        if Path(candidate).is_absolute() or "/" in candidate:
            path = candidate if Path(candidate).exists() else None
        else:
            path = which(candidate)
        if path:
            resolved_candidates.append(path)

    if not resolved_candidates:
        raise FileNotFoundError(
            "Rscript executable not found. Ensure R (>=4.3) is installed or set DEARMETA_RSCRIPT=/path/to/Rscript."
        )

    bitmap_rscript = None
    for path in resolved_candidates:
        if _has_bitmap_capabilities(path):
            bitmap_rscript = path
            break

    chosen = bitmap_rscript or resolved_candidates[0]
    logger.info(
        "Using Rscript at %s (bitmap device: %s)",
        chosen,
        "yes" if _has_bitmap_capabilities(chosen) else "no",
    )
    return chosen


def run_r_analysis(
    gse: str,
    project_root: Path,
    config_path: Path,
    output_root: Path,
    r_script: Path,
    extra_args: Optional[Iterable[str]] = None,
    env: Optional[Dict[str, str]] = None,
) -> None:
    """Execute the R analysis script with the required arguments."""
    rscript_bin = _resolve_rscript_candidates()
    cmd = [
        rscript_bin,
        str(r_script),
        "--gse",
        gse,
        "--project-root",
        str(project_root),
        "--config",
        str(config_path),
        "--output-root",
        str(output_root),
    ]
    if extra_args:
        cmd.extend(list(extra_args))

    logger.info("Running R analysis: %s", " ".join(cmd))
    run_env = os.environ.copy()
    if env:
        run_env.update(env)
    repo_root = Path(__file__).resolve().parent.parent.parent
    run_env.setdefault("RENV_PROJECT", str(repo_root))

    process = subprocess.Popen(
        cmd,
        cwd=str(repo_root),
        env=run_env,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
    )
    captured_lines = []
    if process.stdout is None:
        process.terminate()
        raise RuntimeError("Failed to capture Rscript stdout; subprocess pipe was not created.")
    try:
        for raw_line in process.stdout:
            line = raw_line.rstrip("\n")
            print(line, flush=True)
            captured_lines.append(line)
    finally:
        process.wait()
    if process.returncode != 0:
        logger.error("R analysis failed with exit code %s", process.returncode)
        if captured_lines:
            logger.error("R output:\n%s", "\n".join(captured_lines))
        raise RRuntimeError("R analysis failed; see logs for details.")
    logger.info("R analysis completed successfully.")
