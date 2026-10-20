"""Tests for TinyMapperRunner validation and FASTQ resolution."""

from __future__ import annotations

import re
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


def _make_genome_dir(tmp_path: Path, name: str = "W303") -> Path:
    """Create a minimal fake genome index and return the prefix."""
    genome_dir = tmp_path / "genome"
    genome_dir.mkdir(exist_ok=True)
    prefix = genome_dir / name
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
        Path(str(prefix) + ext).touch()
    return prefix


# ------------------------------------------------------------------ #
#  Single-end / as-single-end validation                               #
# ------------------------------------------------------------------ #


def test_run_se_without_extend_reads_raises(tmp_path):
    """Single-end reads without --extend-reads should error for non-RNA modes."""
    r1 = tmp_path / "JS001_R1.fq.gz"
    r1.touch()
    genome_prefix = _make_genome_dir(tmp_path)
    spec = JobSpec(
        mode="chip",
        sample=str(tmp_path / "JS001"),
        genome=str(genome_prefix),
        output=tmp_path / "results",
        dry_run=True,
    )
    with patch("shutil.which", return_value="/usr/bin/fake"):
        runner = TinyMapperRunner(spec)
        exit_code = runner.run()
    assert exit_code == 1


def test_run_se_with_extend_reads_ok(tmp_path):
    """Single-end chip with --extend-reads set should not raise during validation."""
    r1 = tmp_path / "JS001_R1.fq.gz"
    r1.touch()
    genome_prefix = _make_genome_dir(tmp_path)
    spec = JobSpec(
        mode="chip",
        sample=str(tmp_path / "JS001"),
        genome=str(genome_prefix),
        output=tmp_path / "results",
        extend_reads_len=200,
        dry_run=True,
    )
    with patch("shutil.which", return_value="/usr/bin/fake"):
        runner = TinyMapperRunner(spec)
        exit_code = runner.run()
    assert exit_code == 0


def test_run_shotgun_as_single_end_without_extend_reads_raises(tmp_path):
    """--as-single-end shotgun without --extend-reads should error."""
    r1 = tmp_path / "JS001_R1.fq.gz"
    r2 = tmp_path / "JS001_R2.fq.gz"
    r1.touch()
    r2.touch()
    genome_prefix = _make_genome_dir(tmp_path)
    spec = JobSpec(
        mode="shotgun",
        sample=str(tmp_path / "JS001"),
        genome=str(genome_prefix),
        output=tmp_path / "results",
        as_single_end=True,
        dry_run=True,
    )
    with patch("shutil.which", return_value="/usr/bin/fake"):
        runner = TinyMapperRunner(spec)
        exit_code = runner.run()
    assert exit_code == 1


def test_run_shotgun_paired_end_default_ok(tmp_path):
    """Default paired-end shotgun (no --as-single-end) needs no --extend-reads."""
    r1 = tmp_path / "JS001_R1.fq.gz"
    r2 = tmp_path / "JS001_R2.fq.gz"
    r1.touch()
    r2.touch()
    genome_prefix = _make_genome_dir(tmp_path)
    spec = JobSpec(
        mode="shotgun",
        sample=str(tmp_path / "JS001"),
        genome=str(genome_prefix),
        output=tmp_path / "results",
        dry_run=True,
    )
    with patch("shutil.which", return_value="/usr/bin/fake"):
        runner = TinyMapperRunner(spec)
        exit_code = runner.run()
    assert exit_code == 0


# ------------------------------------------------------------------ #
#  Generated commands (dry-run)                                        #
# ------------------------------------------------------------------ #

_TRACK_TOOLS = ("bamCoverage ", "bamCompare ")

# Pipeline configurations whose tracks are built with deepTools.
_CONFIGS = {
    "chip-pe-input": {"mode": "chip", "with_input": True},
    "chip-pe-calibration": {"mode": "chip", "with_input": True, "with_calibration": True},
    "chip-se-input": {"mode": "chip", "paired": False, "with_input": True, "extend_reads_len": 200},
    "atac-pe": {"mode": "atac"},
    "atac-se": {"mode": "atac", "paired": False, "extend_reads_len": 200},
    "mnase-pe": {"mode": "mnase"},
    "mnase-se": {"mode": "mnase", "paired": False, "extend_reads_len": 200},
    "shotgun-pe": {"mode": "shotgun"},
    "shotgun-se": {"mode": "shotgun", "paired": False, "extend_reads_len": 200},
    "shotgun-as-se": {"mode": "shotgun", "as_single_end": True, "extend_reads_len": 200},
    "rna-pe": {"mode": "rna"},
    "rna-se": {"mode": "rna", "paired": False},
}
_SINGLE_END_CONFIGS = [
    "chip-se-input",
    "atac-se",
    "mnase-se",
    "shotgun-se",
    "shotgun-as-se",
    "rna-se",
]
# Configurations without upstream duplicate removal (duplicates are flagged instead).
_FLAGGING_CONFIGS = [*_SINGLE_END_CONFIGS, "rna-pe"]


def _dry_run_commands(
    tmp_path: Path,
    mode: str,
    *,
    paired: bool = True,
    with_input: bool = False,
    with_calibration: bool = False,
    **spec_kwargs,
) -> list[str]:
    """Dry-run *mode* on empty placeholder files and return the commands it logs."""
    reads = tmp_path / "reads"
    reads.mkdir(exist_ok=True)
    for name in ["JS001", "JS002"] if with_input else ["JS001"]:
        (reads / f"{name}_R1.fq.gz").touch()
        if paired:
            (reads / f"{name}_R2.fq.gz").touch()
    genome_prefix = _make_genome_dir(tmp_path)
    (genome_prefix / "STAR").mkdir(parents=True, exist_ok=True)  # RNA index directory
    spec = JobSpec(
        mode=mode,
        sample=str(reads / "JS001"),
        input=str(reads / "JS002") if with_input else "",
        calibration=str(_make_genome_dir(tmp_path, "CBS138")) if with_calibration else "",
        genome=str(genome_prefix),
        output=tmp_path / "results",
        hash="TSTRUN",
        dry_run=True,
        **spec_kwargs,
    )
    with patch("shutil.which", return_value="/usr/bin/fake"):
        runner = TinyMapperRunner(spec)
        assert runner.run() == 0
    log = runner.paths.logfile.read_text()
    return [line.removeprefix("[EXEC] ") for line in log.splitlines() if line.startswith("[EXEC] ")]


def _track_commands(cmds: list[str]) -> list[str]:
    return [c for c in cmds if c.startswith(_TRACK_TOOLS)]


@pytest.mark.parametrize("config", list(_CONFIGS))
def test_tracks_skip_duplicates_marked_upstream(tmp_path, config):
    """deepTools >= 4.0 has no --ignoreDuplicates: duplicates go through samtools markdup first."""
    cmds = _dry_run_commands(tmp_path, **_CONFIGS[config])
    tracks = _track_commands(cmds)
    assert tracks
    first_track = cmds.index(tracks[0])
    assert any("samtools markdup" in c for c in cmds[:first_track])
    for cmd in tracks:
        assert "--ignoreDuplicates" not in cmd
        assert "--samFlagExclude 1024" in cmd


@pytest.mark.parametrize("config", _FLAGGING_CONFIGS)
def test_single_end_and_rna_flag_duplicates_without_removing(tmp_path, config):
    """Without upstream removal, reads are kept in the BAM and only flagged as duplicates."""
    cmds = _dry_run_commands(tmp_path, **_CONFIGS[config])
    markdup = [c for c in cmds if "samtools markdup" in c]
    assert markdup
    for cmd in markdup:
        assert " -r " not in cmd.split("samtools markdup", 1)[1]


@pytest.mark.parametrize("config", _FLAGGING_CONFIGS)
def test_keep_duplicates_skips_flagging_and_exclusion(tmp_path, config):
    cmds = _dry_run_commands(tmp_path, **_CONFIGS[config], keep_duplicates=True)
    assert not any("samtools markdup" in c for c in cmds)
    assert not any("--samFlagExclude" in c for c in _track_commands(cmds))


@pytest.mark.parametrize("config", _SINGLE_END_CONFIGS)
def test_no_bare_extend_reads_on_single_end(tmp_path, config):
    """deepTools cannot infer a fragment length from single-end reads."""
    for cmd in _track_commands(_dry_run_commands(tmp_path, **_CONFIGS[config])):
        assert not re.search(r"--extendReads(?!\s+\d)", cmd), cmd


def test_single_end_chip_input_compare_tracks_use_fragment_length(tmp_path):
    cmds = _dry_run_commands(tmp_path, **_CONFIGS["chip-se-input"])
    compare = [c for c in cmds if c.startswith("bamCompare ")]
    assert len(compare) == 2
    for cmd in compare:
        assert "--extendReads 200 " in cmd


def test_single_end_rna_stranded_tracks(tmp_path):
    """Without --extend-reads, single-end RNA tracks are not extended at all."""
    tracks = _track_commands(_dry_run_commands(tmp_path, **_CONFIGS["rna-se"]))
    assert len(tracks) == 3
    assert not any("--extendReads" in c for c in tracks)


def test_single_end_rna_stranded_tracks_with_extend_reads(tmp_path):
    tracks = _track_commands(
        _dry_run_commands(tmp_path, **_CONFIGS["rna-se"], extend_reads_len=150)
    )
    stranded = [c for c in tracks if "--filterRNAstrand" in c]
    assert len(stranded) == 2
    for cmd in stranded:
        assert "--extendReads 150 " in cmd
