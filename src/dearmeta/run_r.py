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


def check_rscript_available(rscript: str = "Rscript") -> str:
    """Return the Rscript executable path if available, otherwise raise."""
    from shutil import which

    resolved = which(rscript)
    if not resolved:
        raise FileNotFoundError(
            "Rscript executable not found. Ensure R (>=4.3) is installed or use the Docker workflow."
        )
    return resolved


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
    rscript_bin = check_rscript_available()
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

    process = subprocess.Popen(
        cmd,
        cwd=str(project_root),
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
