"""Tests for R invocation helpers."""

from __future__ import annotations

from pathlib import Path

from dearmeta.run_r import _detect_renv_project


def test_detects_renv_project_from_cwd(tmp_path, monkeypatch):
    marker = tmp_path / "renv.lock"
    marker.write_text("{}")
    monkeypatch.chdir(tmp_path)
    monkeypatch.delenv("RENV_PROJECT", raising=False)
    assert _detect_renv_project() == tmp_path
