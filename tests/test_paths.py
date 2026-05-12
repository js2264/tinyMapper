"""Tests for RunPaths file path naming convention."""

from __future__ import annotations

from pathlib import Path

import pytest

from tinymapper._paths import RunPaths


@pytest.fixture
def paths():
    return RunPaths(
        outdir=Path("/results"),
        sample_base="JS001",
        genome="W303",
        hash_="ABCDEF",
        input_base="JS002",
        spikein="CBS138",
        mnase_min="130",
        mnase_max="200",
        mode="ChIP",
        hic_base_rez="1000",
    )


# ------------------------------------------------------------------ #
#  Naming convention: <sample>^<operation>^<hash>.<ext>               #
# ------------------------------------------------------------------ #


def test_sample_aligned_genome(paths):
    p = paths.sample_aligned_genome
    assert p.name == "JS001^mapped_W303^ABCDEF.sam"
    assert "bam/genome/JS001" in str(p)


def test_sample_aligned_genome_filtered(paths):
    p = paths.sample_aligned_genome_filtered
    assert p.name == "JS001^mapped_W303^filtered^ABCDEF.bam"


def test_sample_non_aligned_genome_prefix(paths):
    p = paths.sample_non_aligned_genome
    assert p.name == "JS001^unmapped_W303^ABCDEF"
    assert "fastq/genome/JS001" in str(p)


def test_input_aligned_genome_filtered(paths):
    p = paths.input_aligned_genome_filtered
    assert p.name == "JS002^mapped_W303^filtered^ABCDEF.bam"
    assert "bam/genome/JS002" in str(p)


def test_calibration_paths(paths):
    p = paths.sample_aligned_calibration
    assert "CBS138" in p.name
    assert "ABCDEF" in p.name
    assert "bam/spikein" in str(p)


def test_sample_non_aligned_calibration_aligned_genome_filtered(paths):
    p = paths.sample_non_aligned_calibration_aligned_genome_filtered
    assert "unmapped_CBS138" in p.name
    assert "mapped_W303" in p.name
    assert "filtered" in p.name
    assert p.suffix == ".bam"


def test_sample_raw_track_no_spikein():
    paths_no_spikein = RunPaths(
        outdir=Path("/results"),
        sample_base="JS001",
        genome="W303",
        hash_="ABCDEF",
        mode="ChIP",
    )
    p = paths_no_spikein.sample_raw_track
    assert p.name == "JS001^mapped_W303^ABCDEF.CPM.bw"
    assert "tracks/JS001" in str(p)


def test_sample_raw_track_with_spikein(paths):
    p = paths.sample_raw_track
    assert "unmapped_CBS138" in p.name
    assert "mapped_W303" in p.name
    assert p.suffix == ".bw"


def test_sample_raw_track_rna():
    paths_rna = RunPaths(
        outdir=Path("/results"),
        sample_base="AB4",
        genome="W303",
        hash_="ABCDEF",
        mode="RNA",
    )
    p = paths_rna.sample_raw_track
    assert "unstranded" in p.name


def test_spikeinscaled_track(paths):
    p = paths.sample_spikeinscaled_track
    assert "calibrated" in p.name
    assert "CBS138" in p.name


def test_mnase_readsize_track(paths):
    p = paths.sample_readsize_track
    assert "130-200" in p.name


def test_mnase_nucpos_track(paths):
    p = paths.sample_nucpos_track
    assert "nucpos" in p.name


def test_mnase_nuccov_track(paths):
    p = paths.sample_nuccov_track
    assert "nuccov" in p.name


def test_mnase_filtered_readsize40(paths):
    p = paths.sample_aligned_genome_filtered_readsize40
    assert "40bp" in p.name


def test_hic_paths(paths):
    assert "matrices" in str(paths.sample_cool)
    assert paths.sample_cool.suffix == ".cool"
    assert paths.sample_mcool.suffix == ".mcool"
    assert paths.sample_hic.suffix == ".hic"


def test_hic_pairs_paths(paths):
    assert "pairs" in str(paths.sample_pairs_valid_idx_pcrfree)
    assert "pcrfree" in paths.sample_pairs_valid_idx_pcrfree.name


def test_logfile_in_logs_dir(paths):
    assert "logs" in str(paths.logfile)
    assert paths.logfile.name.endswith(".log")
    assert paths.errfile.name.endswith(".err")
    assert "ABCDEF" in paths.logfile.name


def test_statfile(paths):
    p = paths.statfile
    assert "stats" in str(p)
    assert "JS001" in p.name
    assert "JS002" in p.name
    assert "W303" in p.name
    assert "CBS138" in p.name
    assert "ABCDEF" in p.name


def test_sample_input_track_with_spikein(paths):
    p = paths.sample_input_track
    assert "JS002" in p.name  # vs-{input_base}
    assert "CBS138" in p.name


def test_sample_input_log2_track_with_spikein(paths):
    p = paths.sample_input_log2_track
    assert "JS002" in p.name
    assert "CBS138" in p.name
    assert ".log2.bw" in p.name
