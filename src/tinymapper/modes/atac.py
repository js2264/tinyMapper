"""ATAC-seq mode: bowtie2 → samtools → bamCoverage → macs2.

ATAC-seq uses the same pipeline as ChIP-seq (paired-end bowtie2,
duplicate marking, CPM tracks, macs2 peaks).  Input and calibration
controls are not used in ATAC-seq; validation in the runner enforces
this.
"""

from __future__ import annotations

import logging
from pathlib import Path

from tinymapper._common import (
    align_bowtie2_paired,
    bam_coverage_cpm,
    call_peaks,
    filter_bam_paired,
)
from tinymapper._paths import RunPaths
from tinymapper.models import JobSpec

logger = logging.getLogger(__name__)


def run(
    spec: JobSpec,
    paths: RunPaths,
    sample_r1: Path,
    sample_r2: Path,
    input_r1: Path | None,
    input_r2: Path | None,
    log_file: Path,
) -> None:
    """Execute the ATAC-seq pipeline (ChIP-like, no input/calibration)."""
    logger.info("Mapping sample reads to reference genome with bowtie2")
    align_bowtie2_paired(
        spec,
        spec.genome,
        sample_r1,
        sample_r2,
        paths.sample_aligned_genome,
        paths.sample_non_aligned_genome,
        log_file,
    )

    logger.info("Filtering sample BAM")
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
    )

    logger.info("Calling peaks")
    out_dir = spec.output / "peaks" / spec.sample_base
    peak_name = f"{spec.sample_base}_genome-{spec.genome_base}_{spec.hash}"
    call_peaks(
        spec,
        paths.sample_aligned_genome_filtered,
        out_dir,
        peak_name,
        log_file,
    )
