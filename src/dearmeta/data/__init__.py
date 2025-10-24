"""Package providing access to bundled R scripts."""

from __future__ import annotations

from contextlib import contextmanager
from importlib.resources import as_file, files
from pathlib import Path
from typing import Iterator


@contextmanager
def _resource_path(name: str) -> Iterator[Path]:
    """Yield a concrete filesystem path for a packaged resource."""
    resource = files(__package__).joinpath(name)
    with as_file(resource) as resolved:
        yield resolved


def analysis_script() -> Iterator[Path]:
    """Context manager yielding the packaged analysis R script path."""
    return _resource_path("analysis.R")


def install_script() -> Iterator[Path]:
    """Context manager yielding the packaged install R script path."""
    return _resource_path("install.R")
