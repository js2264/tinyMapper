"""Tests for CLI entry point and argument parsing."""

from __future__ import annotations

import pytest
from click.testing import CliRunner

from tinymapper.cli import cli


@pytest.fixture
def runner():
    return CliRunner()


def test_help(runner):
    result = runner.invoke(cli, ["--help"])
    assert result.exit_code == 0
    assert "tinyMapper" in result.output
    assert "--mode" in result.output
    assert "--sample" in result.output
    assert "--genome" in result.output


def test_short_help(runner):
    result = runner.invoke(cli, ["-h"])
    assert result.exit_code == 0
    assert "--mode" in result.output


def test_version(runner):
    result = runner.invoke(cli, ["--version"])
    assert result.exit_code == 0
    assert "0.14.19" in result.output


def test_short_version(runner):
    result = runner.invoke(cli, ["-v"])
    assert result.exit_code == 0
    assert "0.14.19" in result.output


def test_missing_required_mode(runner):
    result = runner.invoke(cli, ["--sample", "/x", "--genome", "/y"])
    assert result.exit_code != 0


def test_missing_required_sample(runner):
    result = runner.invoke(cli, ["--mode", "chip", "--genome", "/y"])
    assert result.exit_code != 0


def test_missing_required_genome(runner):
    result = runner.invoke(cli, ["--mode", "chip", "--sample", "/x"])
    assert result.exit_code != 0


def test_invalid_mode(runner):
    result = runner.invoke(cli, ["--mode", "invalid", "--sample", "/x", "--genome", "/y"])
    assert result.exit_code != 0


def test_valid_modes(runner):
    """All six modes should be accepted (pipeline will fail later on missing files)."""
    for mode in [
        "chip",
        "rna",
        "atac",
        "mnase",
        "hic",
        "shotgun",
        "ChIP",
        "RNA",
        "ATAC",
        "MNase",
        "HiC",
    ]:
        result = runner.invoke(
            cli,
            ["--mode", mode, "--sample", "/nonexistent", "--genome", "/nonexistent"],
        )
        # We expect failure (missing files) but NOT an argument-parse failure
        # (exit code 2 is click argument error)
        assert result.exit_code != 2, f"Mode {mode!r} was rejected as invalid"


def test_calibration_without_input_cli_error(runner):
    """CLI should fail with a clear error when --calibration is given without --input."""
    result = runner.invoke(
        cli,
        [
            "--mode",
            "chip",
            "--sample",
            "/reads/JS001",
            "--genome",
            "/g/W303/W303",
            "--calibration",
            "/g/CBS138/CBS138",
        ],
    )
    assert result.exit_code != 0


def test_dry_run_flag_accepted(runner):
    """--dry-run should be accepted (pipeline still fails on missing files)."""
    result = runner.invoke(
        cli,
        [
            "--mode",
            "chip",
            "--sample",
            "/nonexistent",
            "--genome",
            "/nonexistent",
            "--dry-run",
        ],
    )
    # exit code 2 = click arg parse error (not acceptable); 0 or 1 ok
    assert result.exit_code != 2


def test_keepintermediate_flag_accepted(runner):
    result = runner.invoke(
        cli,
        [
            "--mode",
            "chip",
            "--sample",
            "/nonexistent",
            "--genome",
            "/nonexistent",
            "--keepIntermediate",
        ],
    )
    assert result.exit_code != 2


def test_duplicates_flag_accepted(runner):
    result = runner.invoke(
        cli,
        [
            "--mode",
            "chip",
            "--sample",
            "/nonexistent",
            "--genome",
            "/nonexistent",
            "--duplicates",
        ],
    )
    assert result.exit_code != 2


@pytest.mark.parametrize(
    "short,long,value",
    [
        ("-m", "--mode", "chip"),
        ("-s", "--sample", "/reads/JS001"),
        ("-g", "--genome", "/g/W303/W303"),
        ("-t", "--threads", "16"),
        ("-o", "--output", "/tmp/results"),
    ],
)
def test_short_flags_accepted(runner, short, long, value):
    """All short-form flags defined in tinyMapper.sh should work."""
    args_short = [short, value, "--mode", "chip", "--sample", "/x", "--genome", "/y"]
    args_long = [long, value, "--mode", "chip", "--sample", "/x", "--genome", "/y"]
    for args in [args_short, args_long]:
        result = runner.invoke(cli, args)
        assert result.exit_code != 2, f"Flag {args[0]!r} was rejected: {result.output}"
