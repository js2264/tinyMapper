"""ChIP-seq mode: bowtie2 → samtools fixmate/markdup/filter → bamCoverage → macs2.

Supports paired-end and single-end reads.
Paired-end configurations:
  1. Sample only (no input, no calibration)
  2. Sample + input (no calibration)
  3. Sample + input + calibration/spikein
Single-end configurations:
  1. Sample only
  2. Sample + input
  (spikein calibration is not supported for single-end)
"""

from __future__ import annotations

import logging
from pathlib import Path

from tinymapper._common import (
    align_bowtie2_paired,
    align_bowtie2_se,
    bam_compare,
    bam_coverage_cpm,
    bam_coverage_scaled,
    call_peaks,
    compute_scaling_factor,
    filter_bam_paired,
    filter_bam_single,
)
from tinymapper._paths import RunPaths
from tinymapper.models import JobSpec

logger = logging.getLogger(__name__)


def run(
    spec: JobSpec,
    paths: RunPaths,
    sample_r1: Path,
    sample_r2: Path | None,
    input_r1: Path | None,
    input_r2: Path | None,
    log_file: Path,
) -> None:
    """Execute the full ChIP-seq pipeline (paired-end or single-end)."""
    paired = sample_r2 is not None

    # 1. Alignment
    logger.info("Mapping sample reads to reference genome with bowtie2")
    if paired:
        align_bowtie2_paired(
            spec,
            spec.genome,
            sample_r1,
            sample_r2,  # type: ignore[arg-type]
            paths.sample_aligned_genome,
            paths.sample_non_aligned_genome,
            log_file,
        )
    else:
        align_bowtie2_se(
            spec,
            spec.genome,
            sample_r1,
            paths.sample_aligned_genome,
            log_file,
        )

    if spec.do_input:
        logger.info("Mapping input reads to reference genome with bowtie2")
        if paired:
            align_bowtie2_paired(
                spec,
                spec.genome,
                input_r1,
                input_r2,  # type: ignore[arg-type]
                paths.input_aligned_genome,
                paths.input_non_aligned_genome,
                log_file,
            )
        else:
            align_bowtie2_se(
                spec,
                spec.genome,
                input_r1,  # type: ignore[arg-type]
                paths.input_aligned_genome,
                log_file,
            )

    if spec.do_calibration:
        # do_calibration requires paired-end (enforced by runner)
        _align_spikein(spec, paths, sample_r1, sample_r2, input_r1, input_r2, log_file)
        _filter_with_calibration(spec, paths, log_file)
        scaling_factor = compute_scaling_factor(spec, paths, log_file)
        _tracks_with_calibration(spec, paths, scaling_factor, log_file)
    else:
        if paired:
            _filter_no_calibration(spec, paths, log_file)
        else:
            _filter_no_calibration_se(spec, paths, log_file)
        _tracks_no_calibration(spec, paths, log_file)

    # 2. Peak calling
    if spec.do_peaks:
        _peaks(spec, paths, log_file, paired=paired)


# ------------------------------------------------------------------ #
#  Private helpers                                                     #
# ------------------------------------------------------------------ #


def _align_spikein(
    spec: JobSpec,
    paths: RunPaths,
    sample_r1: Path,
    sample_r2: Path | None,
    input_r1: Path | None,
    input_r2: Path | None,
    log_file: Path,
) -> None:
    """Six additional alignments for spikein-based calibration."""
    spikein = spec.calibration
    genome = spec.genome

    logger.info("Mapping sample reads to spikein genome")
    align_bowtie2_paired(
        spec,
        spikein,
        sample_r1,
        sample_r2,
        paths.sample_aligned_calibration,
        paths.sample_non_aligned_calibration,
        log_file,
    )

    logger.info("Mapping sample reads (not on spikein) to reference genome")
    align_bowtie2_paired(
        spec,
        genome,
        Path(f"{paths.sample_non_aligned_calibration}.1.gz"),
        Path(f"{paths.sample_non_aligned_calibration}.2.gz"),
        paths.sample_non_aligned_calibration_aligned_genome,
        paths.sample_non_aligned_genome,  # unmapped-on-spikein prefix reuse
        log_file,
    )

    logger.info("Mapping sample reads (not on genome) to spikein genome")
    align_bowtie2_paired(
        spec,
        spikein,
        Path(f"{paths.sample_non_aligned_genome}.1.gz"),
        Path(f"{paths.sample_non_aligned_genome}.2.gz"),
        paths.sample_non_aligned_genome_aligned_calibration,
        paths.sample_non_aligned_genome,  # prefix reuse (unmapped discarded)
        log_file,
    )

    logger.info("Mapping input reads to spikein genome")
    align_bowtie2_paired(
        spec,
        spikein,
        input_r1,  # type: ignore[arg-type]
        input_r2,  # type: ignore[arg-type]
        paths.input_aligned_calibration,
        paths.input_non_aligned_calibration,
        log_file,
    )

    logger.info("Mapping input reads (not on spikein) to reference genome")
    align_bowtie2_paired(
        spec,
        genome,
        Path(f"{paths.input_non_aligned_calibration}.1.gz"),
        Path(f"{paths.input_non_aligned_calibration}.2.gz"),
        paths.input_non_aligned_calibration_aligned_genome,
        paths.input_non_aligned_genome,
        log_file,
    )

    logger.info("Mapping input reads (not on genome) to spikein genome")
    align_bowtie2_paired(
        spec,
        spikein,
        Path(f"{paths.input_non_aligned_genome}.1.gz"),
        Path(f"{paths.input_non_aligned_genome}.2.gz"),
        paths.input_non_aligned_genome_aligned_calibration,
        paths.input_non_aligned_genome,
        log_file,
    )


def _filter_no_calibration(spec: JobSpec, paths: RunPaths, log_file: Path) -> None:
    logger.info("Filtering sample BAM (no calibration)")
    filter_bam_paired(
        spec, paths.sample_aligned_genome, paths.sample_aligned_genome_filtered, log_file
    )
    if spec.do_input:
        logger.info("Filtering input BAM (no calibration)")
        filter_bam_paired(
            spec, paths.input_aligned_genome, paths.input_aligned_genome_filtered, log_file
        )


def _filter_no_calibration_se(spec: JobSpec, paths: RunPaths, log_file: Path) -> None:
    logger.info("Filtering sample BAM (single-end, no calibration)")
    filter_bam_single(
        spec, paths.sample_aligned_genome, paths.sample_aligned_genome_filtered, log_file
    )
    if spec.do_input:
        logger.info("Filtering input BAM (single-end, no calibration)")
        filter_bam_single(
            spec, paths.input_aligned_genome, paths.input_aligned_genome_filtered, log_file
        )


def _filter_with_calibration(spec: JobSpec, paths: RunPaths, log_file: Path) -> None:
    logger.info("Filtering calibration BAMs (4 files)")
    filter_bam_paired(
        spec,
        paths.sample_non_aligned_calibration_aligned_genome,
        paths.sample_non_aligned_calibration_aligned_genome_filtered,
        log_file,
    )
    filter_bam_paired(
        spec,
        paths.sample_non_aligned_genome_aligned_calibration,
        paths.sample_non_aligned_genome_aligned_calibration_filtered,
        log_file,
    )
    filter_bam_paired(
        spec,
        paths.input_non_aligned_calibration_aligned_genome,
        paths.input_non_aligned_calibration_aligned_genome_filtered,
        log_file,
    )
    filter_bam_paired(
        spec,
        paths.input_non_aligned_genome_aligned_calibration,
        paths.input_non_aligned_genome_aligned_calibration_filtered,
        log_file,
    )


def _tracks_no_calibration(spec: JobSpec, paths: RunPaths, log_file: Path) -> None:
    logger.info("Generating CPM track for sample")
    bam_coverage_cpm(spec, paths.sample_aligned_genome_filtered, paths.sample_raw_track, log_file)
    if spec.do_input:
        logger.info("Generating CPM track for input")
        bam_coverage_cpm(
            spec,
            paths.input_aligned_genome_filtered,
            paths.input_raw_track,
            log_file,
        )
        logger.info("Generating sample/input ratio track")
        bam_compare(
            spec,
            paths.sample_aligned_genome_filtered,
            paths.input_aligned_genome_filtered,
            paths.sample_input_track,
            log_file,
        )


def _tracks_with_calibration(
    spec: JobSpec,
    paths: RunPaths,
    scaling_factor: float,
    log_file: Path,
) -> None:
    sample_bam = paths.sample_non_aligned_calibration_aligned_genome_filtered
    input_bam = paths.input_non_aligned_calibration_aligned_genome_filtered

    logger.info("Generating CPM track for sample (spikein-calibrated run)")
    bam_coverage_cpm(spec, sample_bam, paths.sample_raw_track, log_file)

    logger.info("Generating CPM track for input (spikein-calibrated run)")
    bam_coverage_cpm(spec, input_bam, paths.input_raw_track, log_file)

    logger.info("Generating sample/input ratio track (spikein-calibrated run)")
    bam_compare(spec, sample_bam, input_bam, paths.sample_input_track, log_file)

    logger.info("Generating spikein-scaled CPM track")
    bam_coverage_scaled(
        spec, sample_bam, paths.sample_spikeinscaled_track, scaling_factor, log_file
    )


def _peaks(spec: JobSpec, paths: RunPaths, log_file: Path, *, paired: bool = True) -> None:
    if spec.do_calibration:
        sample_bam = paths.sample_non_aligned_calibration_aligned_genome_filtered
        input_bam: Path | None = paths.input_non_aligned_calibration_aligned_genome_filtered
    else:
        sample_bam = paths.sample_aligned_genome_filtered
        input_bam = paths.input_aligned_genome_filtered if spec.do_input else None

    if input_bam:
        peak_name = f"{spec.sample_base}_vs-{spec.input_base}_genome-{spec.genome_base}_{spec.hash}"
    else:
        peak_name = f"{spec.sample_base}_genome-{spec.genome_base}_{spec.hash}"

    out_dir = spec.output / "peaks" / spec.sample_base
    logger.info("Calling peaks: %s", peak_name)
    call_peaks(spec, sample_bam, out_dir, peak_name, log_file, input_bam=input_bam, paired=paired)
