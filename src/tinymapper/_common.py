"""Shared alignment, filtering and track-generation helpers.

Each mode module calls these helpers rather than building command strings
directly.  All functions accept a ``JobSpec``, a ``RunPaths`` instance,
optional FASTQ paths, and a log file path.
"""

from __future__ import annotations

import logging
import subprocess
from pathlib import Path

from tinymapper._paths import RunPaths
from tinymapper._shell import run_cmd
from tinymapper.models import JobSpec

logger = logging.getLogger(__name__)

# Filter options to use for single-end BAMs (no paired-in-sequencing flags).
_SE_FILTER_OPTS = "-F 0x004 -q 10"


# ------------------------------------------------------------------ #
#  Alignment                                                          #
# ------------------------------------------------------------------ #


def align_bowtie2_paired(
    spec: JobSpec,
    genome_prefix: str,
    r1: Path,
    r2: Path,
    out_sam: Path,
    unmapped_prefix: Path,
    log_file: Path,
) -> None:
    """Run bowtie2 in paired-end mode."""
    cmd = (
        f"bowtie2 {spec.alignment} "
        f"--threads {spec.threads} "
        f"-x {genome_prefix} "
        f"-1 {r1} -2 {r2} "
        f"--un-conc-gz {unmapped_prefix}.gz "
        f"> {out_sam}"
    )
    run_cmd(cmd, log_file, spec.dry_run)


def align_bowtie2_single(
    spec: JobSpec,
    genome_prefix: str,
    r1: Path,
    r2: Path,
    out_sam: Path,
    log_file: Path,
) -> None:
    """Run bowtie2 in single-end mode (shotgun: feeds both R1 and R2 as -U)."""
    cmd = (
        f"bowtie2 {spec.alignment} "
        f"--threads {spec.threads} "
        f"-x {genome_prefix} "
        f"-U {r1},{r2} "
        f"> {out_sam}"
    )
    run_cmd(cmd, log_file, spec.dry_run)


def align_bowtie2_se(
    spec: JobSpec,
    genome_prefix: str,
    r1: Path,
    out_sam: Path,
    log_file: Path,
) -> None:
    """Run bowtie2 with a single read file (true single-end data)."""
    cmd = (
        f"bowtie2 {spec.alignment} --threads {spec.threads} -x {genome_prefix} -U {r1} > {out_sam}"
    )
    run_cmd(cmd, log_file, spec.dry_run)


# ------------------------------------------------------------------ #
#  SAMtools filtering                                                  #
# ------------------------------------------------------------------ #


def filter_bam_paired(
    spec: JobSpec,
    in_sam: Path,
    out_bam: Path,
    log_file: Path,
) -> None:
    """fixmate → sort → markdup → view (filter) → sort → BAM, then index."""
    so = spec.samtools_thread_opts
    rd = spec.remove_duplicates_flag
    fo = spec.filter_opts
    cmd = (
        f"samtools fixmate {so} -m {in_sam} - "
        f"| samtools sort {so} -T {in_sam}_sorting - "
        f"| samtools markdup {so} {rd} - - "
        f"| samtools view {so} {fo} -1 -b - "
        f"| samtools sort {so} -l 9 -T {in_sam}_sorting2 -o {out_bam}"
    )
    run_cmd(cmd, log_file, spec.dry_run)
    run_cmd(f"samtools index -@ {spec.threads} {out_bam}", log_file, spec.dry_run)


def filter_bam_single(
    spec: JobSpec,
    in_sam: Path,
    out_bam: Path,
    log_file: Path,
    filter_opts: str = _SE_FILTER_OPTS,
) -> None:
    """sort → view (filter) → sort → BAM, then index (no fixmate/markdup for SE)."""
    so = spec.samtools_thread_opts
    fo = filter_opts
    cmd = (
        f"samtools sort {so} -T {in_sam}_sorting {in_sam} "
        f"| samtools view {so} {fo} -1 -b - "
        f"| samtools sort {so} -l 9 -T {in_sam}_sorting2 -o {out_bam}"
    )
    run_cmd(cmd, log_file, spec.dry_run)
    run_cmd(f"samtools index -@ {spec.threads} {out_bam}", log_file, spec.dry_run)


# ------------------------------------------------------------------ #
#  Track generation                                                    #
# ------------------------------------------------------------------ #


def bam_coverage_cpm(
    spec: JobSpec,
    bam: Path,
    out_bw: Path,
    log_file: Path,
    extend_reads: bool = True,
    extra_flags: str = "",
) -> None:
    """bamCoverage with CPM normalisation."""
    bl = spec.blacklist_options
    id_ = spec.ignore_duplicates_flag
    extend = "--extendReads" if extend_reads else ""
    cmd = (
        f"bamCoverage "
        f"--bam {bam} "
        f"--outFileName {out_bw} "
        f"--binSize 1 "
        f"--numberOfProcessors {spec.threads} "
        f"{bl} "
        f"--normalizeUsing CPM "
        f"--skipNonCoveredRegions "
        f"{extend} "
        f"{id_} "
        f"{extra_flags}"
    )
    run_cmd(cmd, log_file, spec.dry_run)


def bam_coverage_scaled(
    spec: JobSpec,
    bam: Path,
    out_bw: Path,
    scale_factor: float,
    log_file: Path,
) -> None:
    """bamCoverage with a fixed scale factor (spikein calibration)."""
    bl = spec.blacklist_options
    id_ = spec.ignore_duplicates_flag
    cmd = (
        f"bamCoverage "
        f"--bam {bam} "
        f"--outFileName {out_bw} "
        f"--binSize 1 "
        f"--numberOfProcessors {spec.threads} "
        f"{bl} "
        f"--normalizeUsing None "
        f"--scaleFactor {scale_factor} "
        f"--skipNonCoveredRegions "
        f"--extendReads "
        f"{id_}"
    )
    run_cmd(cmd, log_file, spec.dry_run)


def bam_compare(
    spec: JobSpec,
    sample_bam: Path,
    input_bam: Path,
    out_bw: Path,
    log_file: Path,
) -> None:
    """bamCompare (sample / input ratio track)."""
    bl = spec.blacklist_options
    id_ = spec.ignore_duplicates_flag
    cmd = (
        f"bamCompare "
        f"-b1 {sample_bam} -b2 {input_bam} "
        f"--outFileName {out_bw} "
        f"--scaleFactorsMethod readCount "
        f"--operation ratio "
        f"--skipZeroOverZero "
        f"--skipNAs "
        f"--numberOfProcessors {spec.threads} "
        f"{bl} "
        f"--binSize 5 "
        f"--skipNonCoveredRegions "
        f"--extendReads "
        f"{id_}"
    )
    run_cmd(cmd, log_file, spec.dry_run)


# ------------------------------------------------------------------ #
#  Peak calling                                                        #
# ------------------------------------------------------------------ #


def call_peaks(
    spec: JobSpec,
    sample_bam: Path,
    out_dir: Path,
    peak_name: str,
    log_file: Path,
    input_bam: Path | None = None,
    paired: bool = True,
) -> None:
    """macs2 callpeak.  Uses BAMPE for paired-end data, BAM for single-end."""
    ctrl = f"-c {input_bam} " if input_bam else ""
    fmt = "BAMPE" if paired else "BAM"
    cmd = (
        f"macs2 callpeak "
        f"-t {sample_bam} "
        f"{ctrl}"
        f"--format {fmt} "
        f"--gsize {spec.gsize} "
        f"--outdir {out_dir} "
        f"--name {peak_name}"
    )
    run_cmd(cmd, log_file, spec.dry_run)


# ------------------------------------------------------------------ #
#  Spikein calibration                                                 #
# ------------------------------------------------------------------ #


def _count_paired_reads(bam: Path) -> int:
    """Return the 'paired in sequencing' count from samtools flagstat."""
    result = subprocess.run(
        f"samtools flagstat {bam}",
        shell=True,
        capture_output=True,
        text=True,
        check=True,
    )
    for line in result.stdout.splitlines():
        if "paired in sequencing" in line:
            return int(line.split()[0])
    return 0


def compute_scaling_factor(
    spec: JobSpec,
    paths: RunPaths,
    log_file: Path,
) -> float:
    """Compute spikein-based ORI scaling factor and write stats file.

    Returns 1.0 in dry-run mode (no BAMs to read).
    """
    if spec.dry_run:
        return 1.0

    n_sample_genome = _count_paired_reads(
        paths.sample_non_aligned_calibration_aligned_genome_filtered
    )
    n_input_genome = _count_paired_reads(
        paths.input_non_aligned_calibration_aligned_genome_filtered
    )
    n_sample_spikein = _count_paired_reads(
        paths.sample_non_aligned_genome_aligned_calibration_filtered
    )
    n_input_spikein = _count_paired_reads(
        paths.input_non_aligned_genome_aligned_calibration_filtered
    )

    ori = (n_sample_genome / n_sample_spikein) / (n_input_genome / n_input_spikein)
    scaling_factor = ori / (n_sample_genome / 1_000_000)

    paths.statfile.parent.mkdir(parents=True, exist_ok=True)
    with open(paths.statfile, "w") as fh:
        fh.write(f"{spec.sample_base}_{spec.genome_base}_filtered\t{n_sample_genome}\n")
        fh.write(f"{spec.input_base}_{spec.genome_base}_filtered\t{n_input_genome}\n")
        fh.write(f"{spec.sample_base}_{spec.spikein_base}_filtered\t{n_sample_spikein}\n")
        fh.write(f"{spec.input_base}_{spec.spikein_base}_filtered\t{n_input_spikein}\n")
        fh.write(f"ORI\t{ori}\n")
        fh.write(f"Scaling\t{scaling_factor}\n")

    logger.info("Scaling factor for %s: %.6f (ORI=%.4f)", spec.sample_base, scaling_factor, ori)
    return scaling_factor
