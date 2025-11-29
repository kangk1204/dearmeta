"""CLI-level behavioural tests."""

from __future__ import annotations

import pandas as pd
from pathlib import Path
from typer.testing import CliRunner

from dearmeta import cli


def test_analysis_forwards_min_group_size(tmp_path, monkeypatch):
    runner = CliRunner()
    gse = "GSE123456"
    workspace = tmp_path / gse
    workspace.mkdir()
    (workspace / "runtime").mkdir(parents=True, exist_ok=True)
    configure = workspace / "configure.tsv"
    rows = ["dear_group\tgsm_id\tidat_red\tidat_grn\n"]
    for idx in range(1, 6):
        rows.append(f"A\tGSMA{idx}\ta\tb\n")
    for idx in range(1, 6):
        rows.append(f"B\tGSMB{idx}\ta\tb\n")
    configure.write_text("".join(rows))

    captured = {}

    def fake_run_r_analysis(**kwargs):
        captured.update(kwargs)

    monkeypatch.setattr(cli, "run_r_analysis", fake_run_r_analysis)
    monkeypatch.chdir(tmp_path)
    result = runner.invoke(
        cli.app,
        [
            "analysis",
            "--gse",
            gse,
            "--min-group-size",
            "5",
            "--group-ref",
            "A",
        ],
    )
    assert result.exit_code == 0
    assert "--min-group-size" in captured.get("extra_args", [])
    assert "5" in captured.get("extra_args", [])


def test_analysis_forwards_intersection_choice(tmp_path, monkeypatch):
    runner = CliRunner()
    gse = "GSE654321"
    workspace = tmp_path / gse
    workspace.mkdir()
    (workspace / "runtime").mkdir(parents=True, exist_ok=True)
    configure = workspace / "configure.tsv"
    rows = ["dear_group\tgsm_id\tidat_red\tidat_grn\n", "A\tGSMA1\ta\tb\n", "B\tGSMB1\ta\tb\n"]
    configure.write_text("".join(rows))

    captured = {}

    def fake_run_r_analysis(**kwargs):
        captured.update(kwargs)

    monkeypatch.setattr(cli, "run_r_analysis", fake_run_r_analysis)
    monkeypatch.chdir(tmp_path)
    result = runner.invoke(
        cli.app,
        [
            "analysis",
            "--gse",
            gse,
            "--group-ref",
            "A",
            "--min-group-size",
            "1",
            "--workspace-root",
            str(tmp_path / gse),
            "--intersection-choice",
            "best",
        ],
    )
    assert result.exit_code == 0
    assert "--intersection-choice" in captured.get("extra_args", [])
    assert "best" in captured.get("extra_args", [])


def test_analysis_rejects_bad_intersection_choice(tmp_path):
    runner = CliRunner()
    gse = "GSE222222"
    workspace = tmp_path / gse
    workspace.mkdir()
    configure = workspace / "configure.tsv"
    configure.write_text("dear_group\tgsm_id\nA\tGSMA\nB\tGSMB\n")
    result = runner.invoke(
        cli.app,
        [
            "analysis",
            "--gse",
            gse,
            "--group-ref",
            "A",
            "--min-group-size",
            "1",
            "--workspace-root",
            str(tmp_path / gse),
            "--intersection-choice",
            "invalid",
        ],
    )
    assert result.exit_code != 0
    assert "intersection-choice" in result.stderr.lower() or "intersection-choice" in result.stdout.lower()


def test_subsample_configure_stratifies(tmp_path, monkeypatch):
    runner = CliRunner()
    gse = "GSE000111"
    workspace = tmp_path / gse
    workspace.mkdir()
    configure = workspace / "configure.tsv"

    rows = ["dear_group\tgsm_id\n"]
    for idx in range(6):
        rows.append(f"A\tGSMA{idx}\n")
    for idx in range(4):
        rows.append(f"B\tGSMB{idx}\n")
    configure.write_text("".join(rows))

    monkeypatch.chdir(tmp_path)
    result = runner.invoke(
        cli.app,
        ["subsample-configure", "--gse", gse, "--n", "4", "--seed", "123"],
    )
    assert result.exit_code == 0

    output_path = workspace / "configure_subsampled_4.tsv"
    assert output_path.exists()
    df = pd.read_csv(output_path, sep="\t", comment="#")
    counts = df["dear_group"].value_counts()
    assert counts.loc["A"] == 4
    assert counts.loc["B"] == 4


def test_subsample_configure_errors_when_insufficient(tmp_path, monkeypatch):
    runner = CliRunner()
    gse = "GSE000222"
    workspace = tmp_path / gse
    workspace.mkdir()
    configure = workspace / "configure.tsv"
    rows = ["dear_group\tgsm_id\n", "A\tGSMA1\n", "B\tGSMB1\n"]
    configure.write_text("".join(rows))
    monkeypatch.chdir(tmp_path)
    result = runner.invoke(
        cli.app, ["subsample-configure", "--gse", gse, "--n", "2", "--seed", "1"]
    )
    assert result.exit_code != 0
    assert "fewer samples" in result.output


def test_analysis_accepts_configure_and_output(tmp_path, monkeypatch):
    runner = CliRunner()
    gse = "GSE123000"
    workspace = tmp_path / gse
    workspace.mkdir()
    (workspace / "runtime").mkdir(parents=True, exist_ok=True)

    custom_config = tmp_path / "custom_config.tsv"
    custom_config_rel = tmp_path / "custom_config_rel.tsv"

    rows = [
        "dear_group\tgsm_id\tidat_red\tidat_grn\n",
        "A\tGSMA1\ta\tb\n",
        "A\tGSMA2\ta\tb\n",
        "B\tGSMB1\ta\tb\n",
        "B\tGSMB2\ta\tb\n",
    ]
    custom_config.write_text("".join(rows))
    custom_config_rel.write_text("".join(rows))

    captured = {}

    def fake_run_r_analysis(**kwargs):
        captured.update(kwargs)

    monkeypatch.setattr(cli, "run_r_analysis", fake_run_r_analysis)
    monkeypatch.chdir(tmp_path)
    output_root = tmp_path / "output_dir"
    result = runner.invoke(
        cli.app,
        [
            "analysis",
            "--gse",
            gse,
            "--group-ref",
            "A",
            "--configure",
            custom_config_rel.name,
            "--output",
            str(output_root),
        ],
    )
    assert result.exit_code == 0
    assert captured.get("config_path") == custom_config_rel.resolve()
    assert captured.get("output_root") == output_root.resolve()
