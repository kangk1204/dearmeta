"""Tests for filesystem helper utilities."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from dearmeta import files


def test_ensure_dir_creates_and_returns_path(tmp_path: Path) -> None:
    target = tmp_path / "nested" / "dir"
    result = files.ensure_dir(target)
    assert target.exists()
    assert result == target


def test_write_json_and_text_are_atomic(tmp_path: Path) -> None:
    json_path = tmp_path / "data" / "sample.json"
    files.write_json(json_path, {"answer": 42})
    payload = json.loads(json_path.read_text())
    assert payload["answer"] == 42

    text_path = tmp_path / "logs" / "note.txt"
    files.write_text(text_path, "hello")
    assert text_path.read_text() == "hello"


def test_compute_md5_and_gunzip(tmp_path: Path) -> None:
    gz_path = tmp_path / "content.txt.gz"
    content = b"dearmeta-md5-check"
    files.ensure_dir(gz_path.parent)
    files.ensure_dir(tmp_path)

    import gzip

    with gzip.open(gz_path, "wb") as handle:
        handle.write(content)

    extracted = files.gunzip_file(gz_path)
    assert extracted.read_bytes() == content

    import hashlib

    checksum = files.compute_md5(extracted)
    assert checksum == hashlib.md5(content).hexdigest()


def test_relative_symlink(tmp_path: Path) -> None:
    source = tmp_path / "source.txt"
    source.write_text("example")
    link = tmp_path / "links" / "source.txt"
    files.relative_symlink(source, link)
    assert link.is_symlink()
    assert link.read_text() == "example"
