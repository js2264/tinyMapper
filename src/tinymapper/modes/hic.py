"""Hi-C mode: delegates to ``hicstuff pipeline``, then generates tracks.

Pipeline:
  1. Copy genome FASTA + bowtie2 index to tmp/ (hicstuff needs them co-located)
  2. Run ``hicstuff pipeline`` (alignment + pairing + matrix generation)
  3. Compute distance law file (``hicstuff distancelaw``)
  4. Merge fwd+rev BAMs → sort → markdup → bedtools genomecov → bedGraphToBigWig
  5. Generate .hic file with juicer_tools (optional, skipped if not installed)
  6. Rename / move output files to canonical paths
  7. Clean up tmp artefacts
"""

from __future__ import annotations

import logging
import shutil
import subprocess
from pathlib import Path

from tinymapper._paths import RunPaths
from tinymapper._shell import run_cmd
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
    """Execute the HiC pipeline via hicstuff."""
    _setup_tmp_genome(spec, paths, log_file)
    _run_hicstuff(spec, paths, sample_r1, sample_r2, log_file)
    _distance_law(spec, paths, log_file)
    _coverage_track(spec, paths, log_file)
    _hic_file(spec, paths, log_file)
    _move_outputs(spec, paths, log_file)
    _cleanup_tmp(spec, paths, log_file)


# ------------------------------------------------------------------ #
#  Private helpers                                                     #
# ------------------------------------------------------------------ #


def _setup_tmp_genome(spec: JobSpec, paths: RunPaths, log_file: Path) -> None:
    """Copy genome FASTA + bowtie2 index files to tmp/ for hicstuff."""
    logger.info("Copying genome index files for hicstuff into tmp/")
    paths.hic_tmp_dir.mkdir(parents=True, exist_ok=True)
    genome_base = Path(spec.genome)
    fasta_src = Path(str(genome_base) + ".fa")

    run_cmd(f"cp {fasta_src} {paths.hic_genome_fasta}", log_file, spec.dry_run)
    for ext in [".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"]:
        src = Path(str(genome_base) + ext)
        dst = Path(str(paths.hic_genome_fasta) + ext)
        run_cmd(f"cp {src} {dst}", log_file, spec.dry_run)


def _run_hicstuff(
    spec: JobSpec,
    paths: RunPaths,
    r1: Path,
    r2: Path,
    log_file: Path,
) -> None:
    """Run hicstuff pipeline."""
    prefix = f"{spec.sample_base}^{spec.hash}"
    logger.info("Running hicstuff pipeline")
    cmd = (
        f"hicstuff pipeline "
        f"--threads {spec.threads} "
        f"--enzyme {spec.restriction} "
        f"--outdir {spec.output} "
        f"--prefix {prefix} "
        f"{spec.hicstuff_opts} "
        f"--force "
        f"--binning {spec.binning} "
        f"--exclude Mito,chrM,MT "
        f"--genome {paths.hic_genome_fasta} "
        f"{r1} {r2}"
    )
    run_cmd(cmd, log_file, spec.dry_run)


def _distance_law(spec: JobSpec, paths: RunPaths, log_file: Path) -> None:
    """Generate distance-law (P(s)) plot + TSV."""
    logger.info("Generating distance law file")
    prefix = f"{spec.sample_base}^{spec.hash}"
    opts = spec.hicstuff_opts

    # Determine which pairs file hicstuff produced
    if "--duplicates" in opts or " -d" in opts:
        pairs = spec.output / f"{prefix}.valid_idx_pcrfree.pairs.gz"
    elif "--filter" in opts or " -f" in opts:
        pairs = spec.output / f"{prefix}.valid_idx_filtered.pairs.gz"
    else:
        pairs = spec.output / f"{prefix}.valid_idx.pairs.gz"

    ps_pdf = spec.output / f"{prefix}.ps.pdf"
    ps_tsv = spec.output / f"{prefix}.ps.tsv"
    cmd = f"hicstuff distancelaw --average -o {ps_pdf} -O {ps_tsv} --pairs {pairs}"
    # Non-fatal: continue if this step fails
    if spec.dry_run:
        run_cmd(cmd, log_file, spec.dry_run)
    else:
        try:
            subprocess.run(cmd, shell=True, check=True)
            with open(log_file, "a") as fh:
                fh.write(f"[EXEC] {cmd}\n")
        except subprocess.CalledProcessError:
            logger.warning("hicstuff distancelaw failed — continuing without it")


def _coverage_track(spec: JobSpec, paths: RunPaths, log_file: Path) -> None:
    """Merge fwd+rev BAMs → CPM BigWig track."""
    logger.info("Computing HiC coverage track")
    prefix = f"{spec.sample_base}^{spec.hash}"
    tmp = paths.hic_tmp_dir
    genome_sizes = Path(spec.genome + ".chrom.sizes")

    fwd_bam = spec.output / "tmp" / f"{prefix}.for.bam"
    rev_bam = spec.output / "tmp" / f"{prefix}.rev.bam"
    merged_bam = tmp / f"{prefix}.bam"
    sorted_bam = tmp / f"{prefix}.sorted.bam"
    bedgraph = tmp / f"{prefix}.bg"

    run_cmd(
        f"samtools merge -@ {spec.threads} {merged_bam} {fwd_bam} {rev_bam}",
        log_file,
        spec.dry_run,
    )
    run_cmd(
        f"samtools sort -@ {spec.threads} -T {tmp}/{prefix}_sorting {merged_bam} "
        f"| samtools markdup -@ {spec.threads} -r - - > {sorted_bam}",
        log_file,
        spec.dry_run,
    )

    if not spec.dry_run:
        result = subprocess.run(
            f"samtools stats {sorted_bam} | grep ^SN | grep 'reads mapped:' | cut -f 3",
            shell=True,
            capture_output=True,
            text=True,
            check=True,
        )
        mapped = int(result.stdout.strip())
        cov = 1_000_000 / mapped
    else:
        cov = 1.0

    run_cmd(
        f"bedtools genomecov -bg -scale {cov} -ibam {sorted_bam} | bedtools sort -i - > {bedgraph}",
        log_file,
        spec.dry_run,
    )
    run_cmd(
        f"bedGraphToBigWig {bedgraph} {genome_sizes} {paths.sample_raw_track}",
        log_file,
        spec.dry_run,
    )


def _hic_file(spec: JobSpec, paths: RunPaths, log_file: Path) -> None:
    """Generate .hic file with juicer_tools (optional)."""
    if not shutil.which("juicer_tools"):
        logger.warning("juicer_tools not found — skipping .hic file generation")
        return

    logger.info("Generating .hic file with juicer_tools")
    prefix = f"{spec.sample_base}^{spec.hash}"
    pcrfree = spec.output / f"{prefix}.valid_idx_pcrfree.pairs"
    chr_tsv = spec.output / f"{prefix}.chr.tsv"
    tmp1 = spec.output / "tmp" / f"{spec.hash}_tmp1"
    tmp2 = spec.output / "tmp" / f"{spec.hash}_tmp2"

    run_cmd(
        f"grep -v '^#' {pcrfree} "
        f"| sort -k2,2d -k4,4d "
        f"| sed -e '1 i\\## pairs format v1.0\\n#columns: readID chr1 position1 "
        f"chr2 position2 strand1 strand2' > {tmp1}",
        log_file,
        spec.dry_run,
    )
    run_cmd(
        f"sed '1d' {chr_tsv} | cut -f1,2 > {tmp2}",
        log_file,
        spec.dry_run,
    )
    run_cmd(
        f"juicer_tools pre -r {spec.binning} {tmp1} {paths.sample_hic} {tmp2}",
        log_file,
        spec.dry_run,
    )
    run_cmd(f"rm {tmp1} {tmp2}", log_file, spec.dry_run)


def _move_outputs(spec: JobSpec, paths: RunPaths, log_file: Path) -> None:
    """Rename hicstuff output files to canonical tinyMapper paths."""
    prefix = f"{spec.sample_base}^{spec.hash}"
    opts = spec.hicstuff_opts
    base_rez = spec.hic_base_rez

    def _mv(src: Path, dst: Path) -> None:
        if src.exists():
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.move(str(src), dst)

    # .cool and .mcool
    _mv(spec.output / f"{prefix}_{base_rez}.cool", paths.sample_cool)
    _mv(spec.output / f"{prefix}.mcool", paths.sample_mcool)

    # forward + reverse BAMs
    tmp_prefix = spec.output / "tmp" / prefix
    _mv(Path(f"{tmp_prefix}.for.bam"), paths.sample_aligned_genome_fwd)
    _mv(Path(f"{tmp_prefix}.rev.bam"), paths.sample_aligned_genome_rev)

    # pairs file (which variant hicstuff produced)
    if "--duplicates" in opts or " -d" in opts:
        _mv(
            spec.output / f"{prefix}.valid_idx_pcrfree.pairs.gz",
            Path(str(paths.sample_pairs_valid_idx_pcrfree) + ".gz"),
        )
    elif "--filter" in opts or " -f" in opts:
        _mv(
            spec.output / f"{prefix}.valid_idx_filtered.pairs.gz",
            Path(str(paths.sample_pairs_valid_idx_filtered) + ".gz"),
        )
    else:
        _mv(
            spec.output / f"{prefix}.valid_idx.pairs.gz",
            Path(str(paths.sample_pairs_valid_idx) + ".gz"),
        )

    # frags + plots
    _mv(spec.output / f"{prefix}.frags.tsv", paths.sample_pairs_frags)
    _mv(spec.output / "plots" / f"{prefix}_event_distance.pdf", paths.sample_pairs_dist)
    _mv(spec.output / "plots" / f"{prefix}_frags_hist.pdf", paths.sample_pairs_hist)
    _mv(
        spec.output / "plots" / f"{prefix}_event_distribution.pdf",
        paths.sample_pairs_distr,
    )
    _mv(spec.output / f"{prefix}.ps.pdf", paths.sample_pairs_law)
    _mv(spec.output / f"{prefix}.ps.tsv", paths.sample_pairs_law_tsv)


def _cleanup_tmp(spec: JobSpec, paths: RunPaths, log_file: Path) -> None:
    """Remove hicstuff intermediate files from the output root and tmp/."""
    logger.info("Cleaning up HiC temporary files")
    prefix = f"{spec.sample_base}^{spec.hash}"
    genome_fasta = paths.hic_genome_fasta

    for path in [
        genome_fasta,
        spec.output / f"{prefix}.hicstuff_0.log",
        spec.output / f"{prefix}.chr.tsv",
        spec.output / f"{prefix}.chr.tsv_filtered",
    ]:
        if path.exists():
            path.unlink(missing_ok=True)

    # bt2 index copies
    for ext in [".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"]:
        idx = Path(str(genome_fasta) + ext)
        if idx.exists():
            idx.unlink(missing_ok=True)

    # Sorted / merged BAMs in tmp/
    for fname in [f"{prefix}.bam", f"{prefix}.sorted.bam", f"{prefix}.bg"]:
        p = paths.hic_tmp_dir / fname
        if p.exists():
            p.unlink(missing_ok=True)
