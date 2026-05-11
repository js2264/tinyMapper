"""rich-click CLI — exact long-form flag parity with tinyMapper.sh.

Entry point: ``tinymapper`` (registered in pyproject.toml).
``data/tinyMapper.sh`` is a thin ``exec tinymapper "$@"`` wrapper so the
legacy call surface continues to work unchanged.
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path

import rich_click as click

from tinymapper import __version__
from tinymapper.models import JobSpec, TinyMapperMode
from tinymapper.runner import TinyMapperRunner

# Use a concise help style matching tinyMapper.sh output
click.rich_click.USE_MARKDOWN = False
click.rich_click.SHOW_ARGUMENTS = True
click.rich_click.OPTION_GROUPS = {
    "tinymapper": [
        {
            "name": "Required",
            "options": ["--mode", "--sample", "--genome"],
        },
        {
            "name": "Core optional",
            "options": ["--output", "--input", "--calibration", "--threads"],
        },
        {
            "name": "Alignment / filtering",
            "options": [
                "--alignment",
                "--filter",
                "--blacklist",
                "--gsize",
                "--duplicates",
            ],
        },
        {
            "name": "HiC",
            "options": [
                "--hicstuff",
                "--restriction",
                "--binning",
                "--balance",
            ],
        },
        {
            "name": "MNase",
            "options": ["--MNaseSizes"],
        },
        {
            "name": "Output",
            "options": [
                "--keepIntermediate",
                "--dry-run",
                "--help",
                "--version",
            ],
        },
    ]
}


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.version_option(__version__, "-v", "--version", prog_name="tinyMapper")
# ---- Required ----
@click.option(
    "-m",
    "--mode",
    required=True,
    type=click.Choice(
        [m.value for m in TinyMapperMode],
        case_sensitive=False,
    ),
    help="Mapping mode (ChIP, MNase, ATAC, RNA, HiC, shotgun).",
)
@click.option(
    "-s",
    "--sample",
    required=True,
    help=(
        "Path prefix to sample FASTQ files.  "
        "For ~/reads/JS001_R{1,2}.fq.gz use --sample ~/reads/JS001"
    ),
)
@click.option(
    "-g",
    "--genome",
    required=True,
    help=(
        "Path prefix to reference genome.  "
        "For ~/genome/W303/W303.fa use --genome ~/genome/W303/W303"
    ),
)
# ---- Core optional ----
@click.option(
    "-o",
    "--output",
    default="results",
    show_default=True,
    type=click.Path(path_type=Path),
    help="Directory to store results.",
)
@click.option(
    "-i",
    "--input",
    default="",
    help="(ChIP) Path prefix to input/control sample.",
)
@click.option(
    "-c",
    "--calibration",
    default="",
    help="(ChIP) Path prefix to spikein/calibration genome.",
)
@click.option(
    "-t",
    "--threads",
    default=8,
    show_default=True,
    type=int,
    help="Number of CPU threads.",
)
# ---- Alignment / filtering ----
@click.option(
    "-a",
    "--alignment",
    default="",
    show_default=False,
    help="Extra options passed to bwa-mem2 mem (use single quotes).",
)
@click.option(
    "-f",
    "--filter",
    "filter_opts",
    default="-f 0x001 -f 0x002 -F 0x004 -F 0x008 -q 10",
    show_default=True,
    help="Filtering options for samtools view (use single quotes).",
)
@click.option(
    "-bl",
    "--blacklist",
    default="",
    help="BED file of blacklist regions for bamCoverage.",
)
@click.option(
    "-gs",
    "--gsize",
    default="13000000",
    show_default=True,
    help="Effective genome size for macs3 peak calling.",
)
@click.option(
    "-d",
    "--duplicates",
    is_flag=True,
    default=False,
    help="Keep duplicate reads (default: remove duplicates).",
)
# ---- HiC ----
@click.option(
    "-hic",
    "--hicstuff",
    "hicstuff_opts",
    default="--mapping iterative --duplicates --filter --plot --no-cleanup",
    show_default=True,
    help="Extra arguments passed to hicstuff pipeline.",
)
@click.option(
    "-re",
    "--restriction",
    default="HpaII,HinfI",
    show_default=True,
    help="Restriction enzyme(s) for HiC (e.g. DpnII,HinfI).",
)
@click.option(
    "-b",
    "--binning",
    default="500",
    show_default=True,
    help="Minimum bin resolution for HiC matrix (bp); comma-separated for multi-res.",
)
@click.option(
    "-ba",
    "--balance",
    "balance_opts",
    default="--cis-only --min-nnz 3 --mad-max 7",
    show_default=True,
    help="Balancing options for cooler zoomify.",
)
# ---- MNase ----
@click.option(
    "-M",
    "--MNaseSizes",
    "mnase_sizes",
    default="130,200",
    show_default=True,
    help="Min,Max fragment size for MNase track.",
)
# ---- Output behaviour ----
@click.option(
    "-k",
    "--keepIntermediate",
    "keep_intermediate",
    is_flag=True,
    default=False,
    help="Keep intermediate SAM / unmapped FASTQ files.",
)
@click.option(
    "--dry-run",
    "dry_run",
    is_flag=True,
    default=False,
    help="Log commands without executing them.",
)
def cli(
    mode: str,
    sample: str,
    genome: str,
    output: Path,
    input: str,
    calibration: str,
    threads: int,
    alignment: str,
    filter_opts: str,
    blacklist: str,
    gsize: str,
    duplicates: bool,
    hicstuff_opts: str,
    restriction: str,
    binning: str,
    balance_opts: str,
    mnase_sizes: str,
    keep_intermediate: bool,
    dry_run: bool,
) -> None:
    """tinyMapper — map and process sequencing reads.

    \b
    Modes:
      ChIP    — ChIP-seq (bwa-mem2 → samtools → bamCoverage → macs3)
      RNA     — RNA-seq  (STAR → samtools → bamCoverage × 3)
      ATAC    — ATAC-seq (bwa-mem2 → samtools → bamCoverage → macs3)
      MNase   — MNase-seq (bwa-mem2 → samtools → size filter → 3 tracks)
      HiC     — Hi-C     (hicstuff pipeline → cooler → mcool)
      shotgun — Shotgun  (bwa-mem2 single-end → samtools → bamCoverage)

    \b
    Examples:
      tinymapper -m ChIP -s ~/HB44 -g ~/genomes/R64-1-1/R64-1-1 -o ~/results
      tinymapper -m RNA  -s ~/AB4  -g ~/genomes/W303/W303 -o ~/results
      tinymapper -m HiC  -s ~/CH266 -g ~/genomes/W303/W303 --binning 1000
    """
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        datefmt="%y-%m-%d %H:%M:%S",
    )

    try:
        spec = JobSpec(
            mode=mode,
            sample=sample,
            genome=genome,
            output=output,
            input=input,
            calibration=calibration,
            threads=threads,
            alignment=alignment,
            filter_opts=filter_opts,
            blacklist=blacklist,
            gsize=gsize,
            keep_duplicates=duplicates,
            hicstuff_opts=hicstuff_opts,
            restriction=restriction,
            binning=binning,
            balance_opts=balance_opts,
            mnase_sizes=mnase_sizes,
            keep_intermediate=keep_intermediate,
            dry_run=dry_run,
        )
    except Exception as exc:
        click.echo(f"Error: {exc}", err=True)
        sys.exit(1)

    runner = TinyMapperRunner(spec)
    exit_code = runner.run()
    sys.exit(exit_code)
