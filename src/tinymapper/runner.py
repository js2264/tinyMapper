"""TinyMapperRunner — validates inputs and dispatches to mode modules."""

from __future__ import annotations

import contextlib
import glob
import logging
import re
import shutil
import subprocess
from datetime import datetime
from importlib.metadata import version
from pathlib import Path

from tinymapper._paths import RunPaths
from tinymapper.models import JobSpec, TinyMapperMode

VERSION = version("tinymapper")

logger = logging.getLogger(__name__)

_FASTQ_PATTERNS = [
    ("_R1.fq.gz", "_R2.fq.gz"),
    ("_R1.fastq.gz", "_R2.fastq.gz"),
    ("_nxq_R1.fq.gz", "_nxq_R2.fq.gz"),
    ("_nvq_R1.fq.gz", "_nvq_R2.fq.gz"),
    ("_R1.fq.gz", "_R2.fq.gz"),
    (".end1.fq.gz", ".end2.fq.gz"),
    (".end1.gz", ".end2.gz"),
]

# Single-end fallback patterns (R1 only).
_SINGLE_FASTQ_PATTERNS = [
    "_R1.fq.gz",
    "_R1.fastq.gz",
    "_nxq_R1.fq.gz",
    "_nvq_R1.fq.gz",
    ".fq.gz",
    ".fastq.gz",
]


class TinyMapperRunner:
    """Orchestrates validation, directory setup and mode dispatch."""

    def __init__(self, spec: JobSpec) -> None:
        self.spec = spec
        self.paths = RunPaths(
            outdir=spec.output,
            sample_base=spec.sample_base,
            genome=spec.genome_base,
            hash_=spec.hash,
            input_base=spec.input_base,
            spikein=spec.spikein_base,
            mnase_min=spec.mnase_min_size,
            mnase_max=spec.mnase_max_size,
            mode=spec.mode.value,
            hic_base_rez=spec.hic_base_rez,
        )

    # ------------------------------------------------------------------ #
    #  Public interface                                                    #
    # ------------------------------------------------------------------ #

    def run(self) -> int:
        """Execute the full pipeline.  Returns 0 on success, 1 on failure."""
        spec = self.spec
        paths = self.paths

        self._setup_dirs()
        paths.logfile.parent.mkdir(parents=True, exist_ok=True)
        paths.logfile.touch()
        paths.errfile.touch()
        paths.tmpfile.touch()

        from tinymapper._shell import configure_err_file

        configure_err_file(paths.errfile)

        try:
            self._validate_tools()
            self._validate_genome_index()
            if spec.do_calibration:
                self._validate_genome_index(is_spikein=True)

            sample_r1, sample_r2 = self._resolve_fastq(spec.sample, "sample")
            input_r1: Path | None = None
            input_r2: Path | None = None
            if spec.do_input:
                input_r1, input_r2 = self._resolve_fastq(spec.input, "input")

            paired = sample_r2 is not None
            if not paired:
                if spec.mode == TinyMapperMode.hic:
                    raise RuntimeError(
                        "HiC mode requires paired-end reads (R2 not found for sample prefix)."
                    )
                if spec.do_calibration:
                    raise RuntimeError(
                        "Spikein calibration requires paired-end reads "
                        "(R2 not found for sample prefix)."
                    )
                if spec.mode != TinyMapperMode.rna and spec.extend_reads_len is None:
                    raise RuntimeError(
                        "Single-end reads detected but --extend-reads was not set.  "
                        "Provide a fragment length (e.g. --extend-reads 200) so that "
                        "bamCoverage can extend single-end reads to the correct size."
                    )
                logger.warning(
                    "Single-end reads detected. "
                    "Fragment-size filtering (MNase) and spikein calibration "
                    "are not available in SE mode."
                )
            elif spec.as_single_end and spec.mode == TinyMapperMode.shotgun:
                if spec.extend_reads_len is None:
                    raise RuntimeError(
                        "--as-single-end requires --extend-reads to be set.  "
                        "Provide a fragment length (e.g. --extend-reads 200) so that "
                        "bamCoverage can extend single-end reads to the correct size."
                    )

            self._log_startup(sample_r1)

            from tinymapper.modes import atac, chip, hic, mnase, rna, shotgun

            _MODE_MAP = {
                TinyMapperMode.chip: chip.run,
                TinyMapperMode.rna: rna.run,
                TinyMapperMode.atac: atac.run,
                TinyMapperMode.mnase: mnase.run,
                TinyMapperMode.hic: hic.run,
                TinyMapperMode.shotgun: shotgun.run,
            }
            _MODE_MAP[spec.mode](
                spec=spec,
                paths=paths,
                sample_r1=sample_r1,
                sample_r2=sample_r2,
                input_r1=input_r1,
                input_r2=input_r2,
                log_file=paths.logfile,
            )

            self._log_finish()
            self._cleanup_intermediates()
            return 0

        except Exception as exc:
            logger.error("Pipeline failed: %s", exc)
            self._cleanup_on_failure()
            return 1
        finally:
            paths.tmpfile.unlink(missing_ok=True)
            self._remove_empty_dirs()

    # ------------------------------------------------------------------ #
    #  Directory setup                                                     #
    # ------------------------------------------------------------------ #

    def _setup_dirs(self) -> None:
        spec = self.spec
        outdir = spec.output
        subdirs = [
            "fastq/genome",
            "fastq/spikein",
            "bam/genome",
            "bam/spikein",
            "tracks",
            "peaks",
            "matrices",
            "pairs",
            "plots",
            "stats",
            "logs",
        ]
        for d in subdirs:
            (outdir / d).mkdir(parents=True, exist_ok=True)

        for name in [spec.sample_base, spec.input_base]:
            if not name:
                continue
            for d in [
                f"fastq/genome/{name}",
                f"fastq/spikein/{name}",
                f"bam/genome/{name}",
                f"bam/spikein/{name}",
                f"tracks/{name}",
                f"peaks/{name}",
                f"matrices/{name}",
                f"pairs/{name}",
                f"plots/{name}",
            ]:
                (outdir / d).mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------ #
    #  Validation                                                          #
    # ------------------------------------------------------------------ #

    def _validate_tools(self) -> None:
        spec = self.spec
        if spec.mode == TinyMapperMode.rna:
            required = ["STAR", "samtools", "bamCoverage"]
        elif spec.mode == TinyMapperMode.hic:
            required = ["hicstuff", "samtools", "bedtools", "bedGraphToBigWig"]
        elif spec.mode == TinyMapperMode.shotgun:
            required = ["bwa-mem2", "samtools", "bamCoverage"]
        else:
            required = ["bwa-mem2", "samtools", "bamCoverage"]

        for tool in required:
            if not shutil.which(tool):
                raise RuntimeError(
                    f"{tool} not found in PATH. "
                    f"Install it via: micromamba install -c bioconda {tool}"
                )
        if spec.do_peaks and not shutil.which("macs3"):
            raise RuntimeError(
                "macs3 not found in PATH. Install it via: micromamba install -c bioconda macs3"
            )

    def _validate_genome_index(self, is_spikein: bool = False) -> None:
        spec = self.spec
        genome_prefix = spec.calibration if is_spikein else spec.genome
        genome_base = Path(genome_prefix)
        genome_fa = Path(str(genome_base) + ".fa")
        label = "spikein" if is_spikein else "reference"

        if not genome_fa.exists():
            raise FileNotFoundError(f"{label.capitalize()} genome FASTA not found: {genome_fa}")

        if spec.mode == TinyMapperMode.rna and not is_spikein:
            star_dir = Path(str(genome_base) + "/STAR/")
            if not star_dir.is_dir():
                raise RuntimeError(
                    f"STAR genome index directory not found: {star_dir}\n"
                    f"Build it with: STAR --runMode genomeGenerate "
                    f"--genomeFastaFiles {genome_fa} --genomeDir {star_dir}"
                )
        elif spec.mode == TinyMapperMode.hic:
            # hicstuff calls bowtie2 internally — keep bowtie2 index for this mode.
            for ext in [".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"]:
                idx = Path(str(genome_base) + ext)
                if not idx.exists():
                    raise RuntimeError(
                        f"bowtie2 {label} index file missing: {idx}\n"
                        f"Build it with: bowtie2-build {genome_fa} {genome_base}"
                    )
        else:
            # bwa-mem2 index files are named after the FASTA file.
            for ext in [".0123", ".bwt.2bit.64", ".amb", ".ann", ".pac"]:
                idx = Path(str(genome_fa) + ext)
                if not idx.exists():
                    raise RuntimeError(
                        f"bwa-mem2 {label} index file missing: {idx}\n"
                        f"Build it with: bwa-mem2 index {genome_fa}"
                    )

        # Generate chrom.sizes if missing (non-fatal)
        sizes = Path(str(genome_base) + ".chrom.sizes")
        if not sizes.exists():
            logger.warning("%s chrom.sizes not found; generating from FASTA index", label)
            fai = Path(str(genome_fa) + ".fai")
            subprocess.run(
                f"samtools faidx {genome_fa} && cut -f1-2 {fai} > {sizes}",
                shell=True,
                check=True,
            )

    @staticmethod
    def _resolve_fastq(sample_prefix: str, label: str) -> tuple[Path, Path | None]:
        """Return (R1, R2) paths by trying known FASTQ naming patterns.

        Tries paired-end patterns first.  Falls back to single-end patterns
        (R2 = None) if no paired files are found.
        """
        prefix = Path(sample_prefix)
        for r1_suf, r2_suf in _FASTQ_PATTERNS:
            r1 = prefix.parent / f"{prefix.name}{r1_suf}"
            r2 = prefix.parent / f"{prefix.name}{r2_suf}"
            if r1.exists() and r2.exists():
                if (r1_suf, r2_suf) != _FASTQ_PATTERNS[0]:
                    logger.warning("Found %s FASTQ files: %s & %s", label, r1, r2)
                return r1, r2

        # Try Illumina S## naming pattern
        r1_candidates = [
            f
            for f in glob.glob(str(prefix.parent / f"{prefix.name}*_R1*.gz"))
            if re.search(r"_S\d{1,2}_R1", f)
        ]
        r2_candidates = [
            f
            for f in glob.glob(str(prefix.parent / f"{prefix.name}*_R2*.gz"))
            if re.search(r"_S\d{1,2}_R2", f)
        ]
        if len(r1_candidates) == 1 and len(r2_candidates) == 1:
            logger.warning(
                "Found %s FASTQ files via S## pattern: %s & %s",
                label,
                r1_candidates[0],
                r2_candidates[0],
            )
            return Path(r1_candidates[0]), Path(r2_candidates[0])

        # Fall back to single-end patterns
        for suffix in _SINGLE_FASTQ_PATTERNS:
            r1 = prefix.parent / f"{prefix.name}{suffix}"
            if r1.exists():
                logger.warning(
                    "Found single-end %s FASTQ file: %s (no paired R2 found)",
                    label,
                    r1,
                )
                return r1, None

        raise FileNotFoundError(
            f"{label.capitalize()} FASTQ files not found for prefix: {sample_prefix}\n"
            f"Files must be named e.g. {prefix.name}_R1.fq.gz & {prefix.name}_R2.fq.gz "
            f"(paired-end) or {prefix.name}_R1.fq.gz (single-end)"
        )

    # ------------------------------------------------------------------ #
    #  Logging                                                             #
    # ------------------------------------------------------------------ #

    def _log_startup(self, sample_r1: Path) -> None:
        spec = self.spec
        paths = self.paths
        now = datetime.now().strftime("%y-%m-%d %H:%M:%S")
        lines = [
            f"{now} | [INFO] tinyMapper v{VERSION}",
            f"{now} | [INFO] Hash        : {spec.hash}",
            f"{now} | [INFO] Log file    : {paths.logfile}",
            f"{now} | [INFO] Err file    : {paths.errfile}",
            "---",
            f"{now} | [INFO] MODE        : {spec.mode.value}",
            f"{now} | [INFO] SAMPLE      : {spec.sample}",
            f"{now} | [INFO] GENOME      : {spec.genome}",
        ]
        if spec.do_input:
            lines.append(f"{now} | [INFO] INPUT       : {spec.input}")
        if spec.do_calibration:
            lines.append(f"{now} | [INFO] SPIKEIN     : {spec.calibration}")
        lines += [
            f"{now} | [INFO] Threads     : {spec.threads}",
            f"{now} | [INFO] Output dir  : {spec.output}",
            "---",
        ]
        with open(paths.logfile, "a") as fh:
            fh.write("\n".join(lines) + "\n")

    def _log_finish(self) -> None:
        now = datetime.now().strftime("%y-%m-%d %H:%M:%S")
        with open(self.paths.logfile, "a") as fh:
            fh.write(f"---\n{now} | [INFO] Pipeline completed successfully\n")

    # ------------------------------------------------------------------ #
    #  Cleanup                                                             #
    # ------------------------------------------------------------------ #

    def _cleanup_intermediates(self) -> None:
        """Remove intermediate SAM/unmapped files when keep_intermediate is False."""
        if self.spec.keep_intermediate:
            return
        p = self.paths
        for path in [
            p.sample_aligned_genome,
            p.sample_non_aligned_calibration,
            p.sample_non_aligned_genome_aligned_calibration,
            p.sample_non_aligned_calibration_aligned_genome,
            p.input_aligned_genome,
            p.input_non_aligned_calibration,
            p.input_non_aligned_genome_aligned_calibration,
            p.input_non_aligned_calibration_aligned_genome,
        ]:
            if path.exists():
                path.unlink(missing_ok=True)

        # Unmapped FASTQ gz pairs
        for prefix in [
            p.sample_non_aligned_genome,
            p.sample_non_aligned_calibration,
            p.input_non_aligned_genome,
            p.input_non_aligned_calibration,
        ]:
            for gz in [Path(f"{prefix}.1.gz"), Path(f"{prefix}.2.gz")]:
                gz.unlink(missing_ok=True)

    def _cleanup_on_failure(self) -> None:
        """Remove any files matching this run's hash on failure (skipped with -k)."""
        if self.spec.keep_intermediate:
            return
        for path in self.spec.output.rglob(f"*{self.spec.hash}*"):
            if path.suffix not in (".log", ".err") and path.is_file():
                path.unlink(missing_ok=True)

    def _remove_empty_dirs(self) -> None:
        for d in sorted(self.spec.output.rglob("*"), reverse=True):
            if d.is_dir():
                with contextlib.suppress(OSError):
                    d.rmdir()
