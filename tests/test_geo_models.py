"""Unit tests for GEO data models and helpers."""

from __future__ import annotations

import pytest

from dearmeta.geo import SeriesMetadata, _is_valid_cached_download, parse_series_matrix
from dearmeta.files import compute_md5


def test_series_metadata_infers_platform_and_normalises():
    metadata = SeriesMetadata.parse_obj(
        {
            "gse": "12345",
            "title": "Test",
            "taxon": "Homo sapiens",
            "gpl": ["GPL21145"],
            "gpltitle": ["Illumina EPIC"],
            "n_samples": 5,
        }
    )
    assert metadata.gse == "GSE12345"
    assert metadata.platform_version == "EPICv1"


def test_parse_series_matrix_extracts_samples(tmp_path):
    matrix = tmp_path / "series_matrix.txt"
    matrix.write_text(
        """!Sample_geo_accession\tGSM1\tGSM2\n
!Sample_title\tSample1\tSample2\n
!Sample_organism_ch1\tHomo sapiens\tHomo sapiens\n
!Sample_characteristics_ch1\tage: 30\tage: 40\n
!Sample_supplementary_file\tftp://example.com/GSM1_Red.idat.gz\tftp://example.com/GSM2_Red.idat.gz\n"""
    )
    samples = parse_series_matrix(matrix)
    assert len(samples) == 2
    assert samples[0].characteristics["age"] == "30"


def test_validates_cached_download(tmp_path):
    payload = b"abc123"
    path = tmp_path / "file.gz"
    path.write_bytes(payload)
    expected_md5 = compute_md5(path)
    assert _is_valid_cached_download(path, len(payload), expected_md5)
    assert not _is_valid_cached_download(path, len(payload) + 1, expected_md5)
    assert not _is_valid_cached_download(path, len(payload), "deadbeefdeadbeefdeadbeefdeadbeef")
