"""Tests for TinyMapperRunner validation and FASTQ resolution."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

import pytest

from tinymapper.models import JobSpec
from tinymapper.runner import TinyMapperRunner

# ------------------------------------------------------------------ #
#  Fixtures                                                          #
# ------------------------------------------------------------------ #


@pytest.fixture
def tmp_outdir(tmp_path):
    return tmp_path / "results"


@pytest.fixture
def chip_spec(tmp_outdir):
    return JobSpec(
        mode="chip",
        sample="/reads/JS001",
        genome="/genome/W303/W303",
        output=tmp_outdir,
        hash="TSTRUN",
        dry_run=True,
    )


# ------------------------------------------------------------------ #
#  _setup_dirs                                                         #
# ------------------------------------------------------------------ #


def test_setup_dirs_creates_structure(chip_spec):
    runner = TinyMapperRunner(chip_spec)
    runner._setup_dirs()
    outdir = chip_spec.output
    assert (outdir / "bam/genome").is_dir()
    assert (outdir / "bam/spikein").is_dir()
    assert (outdir / "tracks").is_dir()
    assert (outdir / "logs").is_dir()
    assert (outdir / "stats").is_dir()


def test_setup_dirs_creates_sample_subdir(chip_spec):
    runner = TinyMapperRunner(chip_spec)
    runner._setup_dirs()
    assert (chip_spec.output / "bam/genome/JS001").is_dir()


# ------------------------------------------------------------------ #
#  _resolve_fastq                                                      #
# ------------------------------------------------------------------ #


def test_resolve_fastq_primary_pattern(tmp_path):
    r1 = tmp_path / "JS001_R1.fq.gz"
    r2 = tmp_path / "JS001_R2.fq.gz"
    r1.touch()
    r2.touch()
    found_r1, found_r2 = TinyMapperRunner._resolve_fastq(str(tmp_path / "JS001"), "sample")
    assert found_r1 == r1
    assert found_r2 == r2


def test_resolve_fastq_fastq_gz_pattern(tmp_path):
    r1 = tmp_path / "JS001_R1.fastq.gz"
    r2 = tmp_path / "JS001_R2.fastq.gz"
    r1.touch()
    r2.touch()
    found_r1, found_r2 = TinyMapperRunner._resolve_fastq(str(tmp_path / "JS001"), "sample")
    assert found_r1 == r1
    assert found_r2 == r2


def test_resolve_fastq_end_pattern(tmp_path):
    r1 = tmp_path / "JS001.end1.fq.gz"
    r2 = tmp_path / "JS001.end2.fq.gz"
    r1.touch()
    r2.touch()
    found_r1, found_r2 = TinyMapperRunner._resolve_fastq(str(tmp_path / "JS001"), "sample")
    assert found_r1 == r1
    assert found_r2 == r2


def test_resolve_fastq_not_found_raises(tmp_path):
    with pytest.raises(FileNotFoundError, match="Sample"):
        TinyMapperRunner._resolve_fastq(str(tmp_path / "JS001"), "sample")


def test_resolve_fastq_illumina_sXX_pattern(tmp_path):
    r1 = tmp_path / "JS001_S1_R1_001.fastq.gz"
    r2 = tmp_path / "JS001_S1_R2_001.fastq.gz"
    r1.touch()
    r2.touch()
    found_r1, found_r2 = TinyMapperRunner._resolve_fastq(str(tmp_path / "JS001"), "sample")
    assert found_r1 == r1
    assert found_r2 == r2


# ------------------------------------------------------------------ #
#  _validate_tools                                                     #
# ------------------------------------------------------------------ #


def test_validate_tools_missing_raises(chip_spec):
    runner = TinyMapperRunner(chip_spec)
    with patch("shutil.which", return_value=None) and pytest.raises(
        RuntimeError, match="not found in PATH"
    ):
        runner._validate_tools()


def test_validate_tools_passes_when_available(chip_spec):
    runner = TinyMapperRunner(chip_spec)
    with patch("shutil.which", return_value="/usr/bin/fake"):
        runner._validate_tools()  # should not raise


# ------------------------------------------------------------------ #
#  Runner.run() with dry_run=True                                      #
# ------------------------------------------------------------------ #


def test_run_dry_run_returns_zero(tmp_path):
    """In dry-run mode the pipeline should succeed even without real files."""
    r1 = tmp_path / "JS001_R1.fq.gz"
    r2 = tmp_path / "JS001_R2.fq.gz"
    r1.touch()
    r2.touch()

    # Provide a minimal "genome" prefix with fake bwa-mem2 + fa + chrom.sizes
    genome_dir = tmp_path / "genome"
    genome_dir.mkdir()
    genome_prefix = genome_dir / "W303"
    for ext in [
        ".fa",
        ".fa.fai",
        ".chrom.sizes",
        ".fa.0123",
        ".fa.bwt.2bit.64",
        ".fa.amb",
        ".fa.ann",
        ".fa.pac",
    ]:
        Path(str(genome_prefix) + ext).touch()

    outdir = tmp_path / "results"
    spec = JobSpec(
        mode="chip",
        sample=str(tmp_path / "JS001"),
        genome=str(genome_prefix),
        output=outdir,
        dry_run=True,
    )

    with patch("shutil.which", return_value="/usr/bin/fake"):
        runner = TinyMapperRunner(spec)
        exit_code = runner.run()

    assert exit_code == 0
    assert (outdir / "logs").is_dir()
