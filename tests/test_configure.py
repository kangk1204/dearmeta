"""Tests for configure helpers using lightweight fixtures."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from dearmeta.configure import (
    assemble_metadata_frame,
    generate_configure_tsv,
    harmonise_columns,
    infer_column_types,
)
from dearmeta.geo import IdatPair, SampleMetadata


@pytest.fixture()
def sample_inputs(tmp_path: Path):
    workspace = tmp_path / "workspace"
    workspace.mkdir()
    samples = {}
    idats = {}
    for idx in range(1, 5):
        gsm = f"GSM{idx}"
        sentrix = f"12345678{idx:02d}"
        position = f"R0{idx}C01"
        samples[gsm] = SampleMetadata(
            gsm=gsm,
            title=f"Sample {idx}",
            organism="Homo sapiens",
            characteristics={"Slide": f"{120+idx}", "Batch": "A" if idx <= 2 else "B"},
            supplementary=[],
        )
        red_path = workspace / f"{gsm}_{sentrix}_{position}_Red.idat"
        green_path = workspace / f"{gsm}_{sentrix}_{position}_Grn.idat"
        idats[gsm] = IdatPair(
            gsm=gsm,
            red=str(red_path),
            green=str(green_path),
        )
    for pair in idats.values():
        Path(pair.red).write_text("red")
        Path(pair.green).write_text("green")
    return workspace, samples, idats


def test_assemble_metadata_frame_orders_and_relativises(sample_inputs):
    workspace, samples, idats = sample_inputs
    frame = assemble_metadata_frame(samples, idats, platform_version="EPICv2", workspace_root=workspace)
    assert list(frame["gsm_id"]) == ["GSM1", "GSM2", "GSM3", "GSM4"]
    assert "GSM1" in frame.loc[0, "idat_red"]
    assert frame["sentrix_id"].notna().all()
    assert frame["sentrix_position"].str.startswith("R0").all()


def test_infer_column_types_and_candidate_batches(sample_inputs):
    workspace, samples, idats = sample_inputs
    frame = assemble_metadata_frame(samples, idats, "EPICv2", workspace)
    frame = harmonise_columns(frame)
    frame["array"] = ["01", "02", "01", "02"]
    frame["age"] = [30, 41, 36, 44]
    clean, numeric, batches = infer_column_types(frame)
    assert "age" in numeric
    assert any(col.startswith("array") for col in batches)
    assert clean.shape[0] == 4


def test_generate_configure_tsv(tmp_path: Path):
    data = pd.DataFrame(
        {
            "gsm_id": ["GSM1"],
            "sample_name": ["Sample"],
            "platform_version": ["EPICv1"],
            "species": ["Homo sapiens"],
            "idat_red": ["path/red"],
            "idat_grn": ["path/grn"],
        }
    )
    target = tmp_path / "configure.tsv"
    generate_configure_tsv(data, target, candidate_batches=["array"], numeric_covariates=["age"])
    content = target.read_text()
    assert "candidate_batch_columns" in content
    assert content.endswith("\n")
