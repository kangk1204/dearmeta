"""CLI-level behavioural tests."""

from __future__ import annotations

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
