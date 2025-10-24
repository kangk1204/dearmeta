"""Utility helpers for configuring structured logging across the project."""

from __future__ import annotations

import logging
import sys
from contextlib import contextmanager
from logging.handlers import RotatingFileHandler
from pathlib import Path
from typing import Iterator, Optional

from rich.console import Console
from rich.logging import RichHandler


def get_logger(name: str = "dearmeta", level: int = logging.INFO) -> logging.Logger:
    """Return a logger configured with Rich output suitable for CLI usage."""
    logger = logging.getLogger(name)
    if logger.handlers:
        return logger

    console = Console(stderr=True)
    rich_handler = RichHandler(
        console=console,
        rich_tracebacks=True,
        markup=False,
        show_path=False,
        omit_repeated_times=False,
    )
    formatter = logging.Formatter(fmt="%(asctime)s %(levelname)s %(message)s", datefmt="%Y-%m-%d %H:%M:%S")
    rich_handler.setFormatter(formatter)

    logger.addHandler(rich_handler)
    logger.setLevel(level)
    logger.propagate = False
    return logger


def ensure_log_file_handler(logger: logging.Logger, log_path: Path, max_bytes: int = 2 * 1024 * 1024) -> None:
    """Attach a log file handler writing to ``log_path`` if not already present."""
    log_path.parent.mkdir(parents=True, exist_ok=True)
    for handler in logger.handlers:
        if isinstance(handler, RotatingFileHandler) and Path(handler.baseFilename) == log_path:
            return

    file_handler = RotatingFileHandler(log_path, maxBytes=max_bytes, backupCount=3, encoding="utf-8")
    file_formatter = logging.Formatter(
        fmt="%(asctime)s %(levelname)s %(name)s %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    file_handler.setFormatter(file_formatter)
    logger.addHandler(file_handler)


@contextmanager
def temporary_log_level(logger: logging.Logger, level: int) -> Iterator[None]:
    """Context manager to temporarily set the log level."""
    original = logger.level
    logger.setLevel(level)
    try:
        yield
    finally:
        logger.setLevel(original)


def configure_stderr_logging(level: int = logging.INFO) -> logging.Logger:
    """Convenience function to configure and return the root project logger."""
    return get_logger(level=level)


def silence_external_loggers(level: int = logging.WARNING) -> None:
    """Reduce noise from third-party libraries."""
    for name in ("urllib3", "requests", "pandas", "fsspec"):
        logging.getLogger(name).setLevel(level)

