"""Shotgun sequencing mode: bwa-mem2 → samtools → bamCoverage.

Default (paired-end) pipeline:
  1. bwa-mem2 paired-end alignment (R1 + R2)
  2. samtools fixmate → markdup → filter → BAM
  3. bamCoverage CPM track (--extendReads inferred from fragment size)

Legacy single-end pipeline (--as-single-end):
  1. bwa-mem2 single-end: R1 and R2 concatenated and mapped independently
  2. samtools sort → filter → BAM (no fixmate/markdup)
  3. bamCoverage CPM track (--extendReads from --extend-reads)
"""

from __future__ import annotations

import logging
from pathlib import Path

from tinymapper._common import (
    align_bwamem2_paired,
    align_bwamem2_se,
    bam_coverage_cpm,
    filter_bam_paired,
    filter_bam_single,
)
from tinymapper._paths import RunPaths
from tinymapper._shell import run_cmd
from tinymapper.models import JobSpec

logger = logging.getLogger(__name__)

_SHOTGUN_FILTER = "-F 0x004 -q 10"


def run(
    spec: JobSpec,
    paths: RunPaths,
    sample_r1: Path,
    sample_r2: Path | None,
    input_r1: Path | None,
    input_r2: Path | None,
    log_file: Path,
) -> None:
    """Execute the shotgun pipeline (paired-end by default, SE with --as-single-end)."""
    if spec.as_single_end:
        logger.info(
            "Mapping sample reads to reference genome with bwa-mem2 "
            "(R1+R2 as single-end, legacy --as-single-end mode)"
        )
        _align_as_single(spec, paths, sample_r1, sample_r2, log_file)
        logger.info("Filtering sample BAM (single-end)")
        _filter_shotgun_se(spec, paths, log_file)
        logger.info("Generating CPM track for sample")
        bam_coverage_cpm(
            spec,
            paths.sample_aligned_genome_filtered,
            paths.sample_raw_track,
            log_file,
            extend_reads=spec.extend_reads_len,
        )
    else:
        if sample_r2 is not None:
            logger.info("Mapping sample reads to reference genome with bwa-mem2 (paired-end)")
            align_bwamem2_paired(
                spec,
                spec.genome,
                sample_r1,
                sample_r2,
                paths.sample_aligned_genome,
                paths.sample_non_aligned_genome,
                log_file,
            )
            logger.info("Filtering sample BAM (paired-end)")
            filter_bam_paired(
                spec,
                paths.sample_aligned_genome,
                paths.sample_aligned_genome_filtered,
                log_file,
            )
            logger.info("Generating CPM track for sample")
            bam_coverage_cpm(
                spec,
                paths.sample_aligned_genome_filtered,
                paths.sample_raw_track,
                log_file,
                extend_reads=True,
            )
        else:
            logger.info(
                "Mapping sample reads to reference genome with bwa-mem2 "
                "(single FASTQ file, single-end)"
            )
            align_bwamem2_se(
                spec,
                spec.genome,
                sample_r1,
                paths.sample_aligned_genome,
                log_file,
            )
            logger.info("Filtering sample BAM (single-end)")
            filter_bam_single(
                spec,
                paths.sample_aligned_genome,
                paths.sample_aligned_genome_filtered,
                log_file,
                filter_opts=_SHOTGUN_FILTER,
            )
            logger.info("Generating CPM track for sample")
            bam_coverage_cpm(
                spec,
                paths.sample_aligned_genome_filtered,
                paths.sample_raw_track,
                log_file,
                extend_reads=spec.extend_reads_len,
            )


# ------------------------------------------------------------------ #
#  Private helpers (legacy single-end path)                           #
# ------------------------------------------------------------------ #


def _align_as_single(
    spec: JobSpec,
    paths: RunPaths,
    r1: Path,
    r2: Path | None,
    log_file: Path,
) -> None:
    """bwa-mem2 in single-end mode: R1 (and optionally R2) treated as independent reads."""
    genome_fa = f"{spec.genome}.fa"
    aln = f" {spec.alignment}" if spec.alignment else ""
    if r2 is not None:
        cmd = (
            f"cat {r1} {r2} | "
            f"bwa-mem2 mem{aln} -v 1 -t {spec.threads} {genome_fa} - "
            f"> {paths.sample_aligned_genome}"
        )
    else:
        cmd = (
            f"bwa-mem2 mem{aln} -v 1 -t {spec.threads} {genome_fa} {r1} "
            f"> {paths.sample_aligned_genome}"
        )
    run_cmd(cmd, log_file, spec.dry_run)


def _filter_shotgun_se(spec: JobSpec, paths: RunPaths, log_file: Path) -> None:
    """Sort → view (shotgun SE filter flags) → sort → BAM, then index."""
    filter_bam_single(
        spec,
        paths.sample_aligned_genome,
        paths.sample_aligned_genome_filtered,
        log_file,
        filter_opts=_SHOTGUN_FILTER,
    )
