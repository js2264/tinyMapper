"""Shotgun sequencing mode: bowtie2 (single-end) → samtools → bamCoverage.

Shotgun mode overrides alignment options to ``--sensitive-local`` and
filter options to ``-F 0x004 -q 10`` (single-end flags), matching the
behaviour of ``tinyMapper.sh`` for this mode.
"""

from __future__ import annotations

import logging
from pathlib import Path

from tinymapper._common import bam_coverage_cpm
from tinymapper._paths import RunPaths
from tinymapper._shell import run_cmd
from tinymapper.models import JobSpec

logger = logging.getLogger(__name__)

_SHOTGUN_ALIGNMENT = "--sensitive-local"
_SHOTGUN_FILTER = "-F 0x004 -q 10"


def run(
    spec: JobSpec,
    paths: RunPaths,
    sample_r1: Path,
    sample_r2: Path,
    input_r1: Path | None,
    input_r2: Path | None,
    log_file: Path,
) -> None:
    """Execute the shotgun pipeline."""
    logger.info("Mapping sample reads to reference genome with bowtie2 (single-end)")
    _align_single(spec, paths, sample_r1, sample_r2, log_file)

    logger.info("Filtering sample BAM (single-end)")
    # Override filter options for single-end
    _filter_shotgun(spec, paths, log_file)

    logger.info("Generating CPM track for sample")
    bam_coverage_cpm(
        spec,
        paths.sample_aligned_genome_filtered,
        paths.sample_raw_track,
        log_file,
        extend_reads=False,
    )


# ------------------------------------------------------------------ #
#  Private helpers                                                     #
# ------------------------------------------------------------------ #


def _align_single(
    spec: JobSpec,
    paths: RunPaths,
    r1: Path,
    r2: Path,
    log_file: Path,
) -> None:
    """bowtie2 in single-end mode: R1 and R2 both passed as -U."""
    cmd = (
        f"bowtie2 {_SHOTGUN_ALIGNMENT} "
        f"--threads {spec.threads} "
        f"-x {spec.genome} "
        f"-U {r1},{r2} "
        f"> {paths.sample_aligned_genome}"
    )
    run_cmd(cmd, log_file, spec.dry_run)


def _filter_shotgun(spec: JobSpec, paths: RunPaths, log_file: Path) -> None:
    """Sort → view (shotgun filter flags) → sort → BAM, then index."""
    so = spec.samtools_thread_opts
    in_sam = paths.sample_aligned_genome
    out_bam = paths.sample_aligned_genome_filtered
    cmd = (
        f"samtools sort {so} -T {in_sam}_sorting {in_sam} "
        f"| samtools view {so} {_SHOTGUN_FILTER} -1 -b - "
        f"| samtools sort {so} -l 9 -T {in_sam}_sorting2 -o {out_bam}"
    )
    run_cmd(cmd, log_file, spec.dry_run)
    run_cmd(f"samtools index -@ {spec.threads} {out_bam}", log_file, spec.dry_run)
