"""Filesystem utilities supporting resilient downloads and artifact management."""

from __future__ import annotations

import gzip
import hashlib
import json
import shutil
import time
from pathlib import Path
from typing import Optional

import requests
from requests import Response
from tqdm import tqdm

CHUNK_SIZE = 1024 * 1024  # 1 MiB


class DownloadError(RuntimeError):
    """Raised when a download fails after retries."""


def ensure_dir(path: Path) -> Path:
    """Create ``path`` directory if needed and return it."""
    path.mkdir(parents=True, exist_ok=True)
    return path


def _stream_to_file(response: Response, destination: Path, show_progress: bool = True) -> str:
    """Stream an HTTP response to disk with optional progress bar; return md5 checksum."""
    total = int(response.headers.get("Content-Length", 0))
    checksum = hashlib.md5()

    if show_progress:
        progress = tqdm(total=total, unit="B", unit_scale=True, desc=f"Downloading {destination.name}")
    else:
        progress = None

    with destination.open("wb") as fh:
        for chunk in response.iter_content(chunk_size=CHUNK_SIZE):
            if not chunk:
                continue
            fh.write(chunk)
            checksum.update(chunk)
            if progress:
                progress.update(len(chunk))

    if progress:
        progress.close()

    return checksum.hexdigest()


def download_file(
    url: str,
    destination: Path,
    retries: int = 4,
    backoff_factor: float = 1.5,
    timeout: int = 60,
    show_progress: bool = True,
) -> str:
    """Download a file from ``url`` to ``destination``; return MD5 checksum."""
    ensure_dir(destination.parent)
    attempt = 0
    while True:
        try:
            response = requests.get(url, stream=True, timeout=timeout)
            response.raise_for_status()
            checksum = _stream_to_file(response, destination, show_progress=show_progress)
            return checksum
        except requests.RequestException as exc:  # pragma: no cover - error path
            attempt += 1
            if attempt > retries:
                raise DownloadError(f"Failed to download {url}: {exc}") from exc
            sleep_time = backoff_factor ** attempt
            time.sleep(sleep_time)


def write_json(path: Path, data: dict) -> None:
    """Write a dictionary to JSON with UTF-8 encoding."""
    ensure_dir(path.parent)
    tmp_path = path.with_suffix(path.suffix + ".tmp")
    with tmp_path.open("w", encoding="utf-8") as fh:
        json.dump(data, fh, indent=2, sort_keys=True)
    tmp_path.replace(path)


def write_text(path: Path, content: str) -> None:
    """Atomically persist text content."""
    ensure_dir(path.parent)
    tmp_path = path.with_suffix(path.suffix + ".tmp")
    with tmp_path.open("w", encoding="utf-8") as fh:
        fh.write(content)
    tmp_path.replace(path)


def gunzip_file(path: Path, remove_original: bool = False) -> Path:
    """Gunzip ``path`` if it ends with .gz; return the extracted file path."""
    if path.suffix != ".gz":
        return path
    target = path.with_suffix("")
    with gzip.open(path, "rb") as src, target.open("wb") as dst:
        shutil.copyfileobj(src, dst)
    if remove_original:
        path.unlink()
    return target


def compute_md5(path: Path) -> str:
    """Compute MD5 checksum for a file."""
    checksum = hashlib.md5()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(CHUNK_SIZE), b""):
            checksum.update(chunk)
    return checksum.hexdigest()


def relative_symlink(source: Path, link_name: Path) -> None:
    """Create a relative symlink from ``link_name`` to ``source``."""
    ensure_dir(link_name.parent)
    if link_name.exists() or link_name.is_symlink():
        link_name.unlink()
    link_name.symlink_to(source.resolve())

