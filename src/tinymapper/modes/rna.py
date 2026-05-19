"""RNA-seq mode: STAR alignment → samtools sort → bamCoverage (3 tracks).

Supports paired-end and single-end reads.
Produces:
  - unstranded CPM BigWig
  - forward-strand CPM BigWig
  - reverse-strand CPM BigWig
"""

from __future__ import annotations

import logging
from pathlib import Path

from tinymapper._common import bam_coverage_cpm
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
    """Execute the RNA-seq pipeline (paired-end or single-end)."""
    _align_star(spec, paths, sample_r1, sample_r2, log_file)
    _sort_bam(spec, paths, log_file)
    _rna_tracks(spec, paths, log_file)


# ------------------------------------------------------------------ #
#  Private helpers                                                     #
# ------------------------------------------------------------------ #


def _align_star(
    spec: JobSpec,
    paths: RunPaths,
    r1: Path,
    r2: Path | None,
    log_file: Path,
) -> None:
    """Run STAR; supports paired-end (r1 + r2) and single-end (r1 only)."""
    genome_dir = f"{spec.genome}/STAR/"
    out_bam = paths.sample_aligned_genome_bam
    out_bam.parent.mkdir(parents=True, exist_ok=True)
    read_files = f"{r1} {r2}" if r2 is not None else str(r1)
    logger.info("Mapping sample reads to reference genome with STAR")
    cmd = (
        f"STAR "
        f"--genomeDir {genome_dir} "
        f"--readFilesCommand zcat "
        f"--runThreadN {spec.threads} "
        f"--readFilesIn {read_files} "
        f"--outFileNamePrefix {out_bam}. "
        f"--outSAMtype BAM Unsorted "
        f"--outSAMunmapped None "
        f"--outSAMattributes Standard"
    )
    run_cmd(cmd, log_file, spec.dry_run)
    # STAR appends `.Aligned.out.bam`; rename to the canonical path
    run_cmd(f"mv {out_bam}.Aligned.out.bam {out_bam}", log_file, spec.dry_run)


def _sort_bam(spec: JobSpec, paths: RunPaths, log_file: Path) -> None:
    """Sort and index the STAR output BAM."""
    bam = paths.sample_aligned_genome_bam
    filtered = paths.sample_aligned_genome_filtered
    filtered.parent.mkdir(parents=True, exist_ok=True)
    so = spec.samtools_thread_opts
    logger.info("Sorting and indexing BAM")
    cmd = f"samtools sort {so} -l 9 -T {bam}_sorting -o {filtered} {bam}"
    run_cmd(cmd, log_file, spec.dry_run)
    run_cmd(f"samtools index -@ {spec.threads} {filtered}", log_file, spec.dry_run)


def _rna_tracks(spec: JobSpec, paths: RunPaths, log_file: Path) -> None:
    """Generate unstranded, forward-strand and reverse-strand BigWig tracks."""
    bam = paths.sample_aligned_genome_filtered

    logger.info("Generating unstranded CPM track")
    bam_coverage_cpm(spec, bam, paths.sample_raw_track, log_file, extend_reads=False)

    logger.info("Generating forward-strand CPM track")
    bam_coverage_cpm(
        spec,
        bam,
        paths.sample_track_fwd,
        log_file,
        extra_flags="--filterRNAstrand forward",
    )

    logger.info("Generating reverse-strand CPM track")
    bam_coverage_cpm(
        spec,
        bam,
        paths.sample_track_rev,
        log_file,
        extra_flags="--filterRNAstrand reverse",
    )
