"""MNase-seq mode: bwa-mem2 → samtools → [fragment-size filter] → bamCoverage.

Paired-end pipeline:
  1. bwa-mem2 paired-end alignment
  2. samtools fixmate → markdup → filter → BAM
  3. Filter BAM by insert size (MNASE_MIN–MNASE_MAX bp)
  4. Generate 40 bp-resized BAM (for nucleosome positioning)
  5. bamCoverage: size-filtered CPM, nucleosome-positioning (--MNase), coverage

Single-end pipeline (fragment size is unknown):
  1. bwa-mem2 single-end alignment
  2. samtools sort → filter → BAM
  3. bamCoverage: standard CPM track only
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


def run(
    spec: JobSpec,
    paths: RunPaths,
    sample_r1: Path,
    sample_r2: Path | None,
    input_r1: Path | None,
    input_r2: Path | None,
    log_file: Path,
) -> None:
    """Execute the MNase-seq pipeline (paired-end or single-end)."""
    paired = sample_r2 is not None

    logger.info("Mapping sample reads to reference genome with bwa-mem2")
    if paired:
        align_bwamem2_paired(
            spec,
            spec.genome,
            sample_r1,
            sample_r2,  # type: ignore[arg-type]
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
        _filter_by_size(spec, paths, log_file)
        _generate_40bp_bam(spec, paths, log_file)
        _mnase_tracks(spec, paths, log_file)
    else:
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
        )
        logger.info(
            "Single-end MNase: skipping fragment-size filter and nucleosome tracks; "
            "generating standard CPM track only"
        )
        bam_coverage_cpm(
            spec,
            paths.sample_aligned_genome_filtered,
            paths.sample_raw_track,
            log_file,
            extend_reads=spec.extend_reads_len,
        )


# ------------------------------------------------------------------ #
#  Private helpers                                                     #
# ------------------------------------------------------------------ #


def _filter_by_size(spec: JobSpec, paths: RunPaths, log_file: Path) -> None:
    """Keep only reads with insert size in [MNASE_MIN, MNASE_MAX]."""
    min_sz = spec.mnase_min_size
    max_sz = spec.mnase_max_size
    logger.info("Filtering BAM by fragment size (%s–%s bp)", min_sz, max_sz)

    in_bam = paths.sample_aligned_genome_filtered
    out_bam = paths.sample_aligned_genome_filtered_readsize

    cmd = (
        f"samtools view -@ {spec.threads} -h {in_bam} "
        f"| mawk '/^@/ || (sqrt(($9^2)) > {min_sz} && sqrt(($9^2)) < {max_sz})' "
        f"| samtools view -b - > {out_bam}"
    )
    run_cmd(cmd, log_file, spec.dry_run)
    run_cmd(f"samtools index -@ {spec.threads} {out_bam}", log_file, spec.dry_run)


def _generate_40bp_bam(spec: JobSpec, paths: RunPaths, log_file: Path) -> None:
    """Resize fragments to 40 bp centred on the fragment midpoint."""
    logger.info("Generating 40 bp-resized BAM")
    so = spec.samtools_thread_opts
    in_bam = paths.sample_aligned_genome_filtered_readsize
    out_bam = paths.sample_aligned_genome_filtered_readsize40

    # Shift read start and set TLEN to ±40 using mawk
    cmd = (
        f"samtools view -@ {spec.threads} -h {in_bam} "
        r"| mawk -F'\t' -v OFS='\t' '{ if ($1 ~ /^@/) { print $0 } else { "
        r"if ($9 > 0) { $4 = $4 + int(($9 - 40) / 2); $9 = 40 } "
        r"else { $4 = $4 - int(($9 + 40) / 2); $9 = -40 } print $0 } }' "
        f"| samtools view -b - "
        f"| samtools sort {so} --write-index -l 9 "
        f"-T {out_bam}_resized_sorting -o {out_bam}"
    )
    run_cmd(cmd, log_file, spec.dry_run)


def _mnase_tracks(spec: JobSpec, paths: RunPaths, log_file: Path) -> None:
    """Generate size-filtered CPM, nucleosome-positioning and coverage tracks."""
    bam = paths.sample_aligned_genome_filtered_readsize

    logger.info("Generating size-filtered CPM track")
    bam_coverage_cpm(spec, bam, paths.sample_readsize_track, log_file)

    logger.info("Generating nucleosome-positioning track (--MNase)")
    bam_coverage_cpm(
        spec,
        bam,
        paths.sample_nucpos_track,
        log_file,
        extra_flags="--smoothLength 10 --MNase",
    )

    logger.info("Generating nucleosome-coverage track (--centerReads)")
    bam_coverage_cpm(
        spec,
        bam,
        paths.sample_nuccov_track,
        log_file,
        extra_flags="--centerReads",
    )
