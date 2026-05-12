"""Typed job specification for tinyMapper.

All parameters map 1-to-1 to tinyMapper.sh CLI flags.
"""

from __future__ import annotations

import random
import string
from enum import Enum
from pathlib import Path

from pydantic import BaseModel, field_validator, model_validator


class TinyMapperMode(str, Enum):
    """Supported mapping modes (matches tinyMapper.sh case-sensitive values)."""

    chip = "ChIP"
    rna = "RNA"
    atac = "ATAC"
    mnase = "MNase"
    hic = "HiC"
    shotgun = "shotgun"


_MODE_ALIASES: dict[str, str] = {
    "chip": "ChIP",
    "rna": "RNA",
    "atac": "ATAC",
    "mnase": "MNase",
    "hic": "HiC",
    "shotgun": "shotgun",
}


def _generate_hash() -> str:
    """Generate a 6-character uppercase-alphanumeric hash (matches tinyMapper.sh)."""
    return "".join(random.choices(string.ascii_uppercase + string.digits, k=6))


class JobSpec(BaseModel):
    """Complete specification for a single tinyMapper run.

    All field names and defaults mirror tinyMapper.sh variables.
    """

    model_config = {"populate_by_name": True}

    # ---- Required ----
    mode: TinyMapperMode
    sample: str
    genome: str

    # ---- Core optional ----
    output: Path = Path("results")
    input: str = ""  # --input / -i
    calibration: str = ""  # --calibration / -c  (spikein genome prefix)
    threads: int = 8

    # ---- Alignment / filtering ----
    alignment: str = ""
    filter_opts: str = "-f 0x001 -f 0x002 -F 0x004 -F 0x008 -q 10"
    blacklist: str = ""
    gsize: str = "13000000"
    keep_duplicates: bool = False  # --duplicates flag
    keep_intermediate: bool = False  # --keepIntermediate flag

    # ---- HiC-specific ----
    hicstuff_opts: str = "--mapping iterative --duplicates --filter --plot --no-cleanup"
    restriction: str = "HpaII,HinfI"
    binning: str = "500"
    balance_opts: str = "--cis-only --min-nnz 3 --mad-max 7"

    # ---- MNase-specific ----
    mnase_sizes: str = "130,200"

    # ---- Track generation ----
    extend_reads_len: int | None = None  # --extendReads value for single-end bamCoverage

    # ---- Runtime (not CLI flags) ----
    hash: str = ""  # populated in model_post_init
    dry_run: bool = False

    def model_post_init(self, __context: object) -> None:
        if not self.hash:
            object.__setattr__(self, "hash", _generate_hash())

    # ------------------------------------------------------------------ #
    #  Validators                                                          #
    # ------------------------------------------------------------------ #

    @field_validator("mode", mode="before")
    @classmethod
    def normalise_mode(cls, v: str) -> str:
        """Accept case-insensitive mode names (chip → ChIP, hic → HiC, etc.)."""
        return _MODE_ALIASES.get(str(v).lower(), v)

    @model_validator(mode="after")
    def _check_calibration_requires_input(self) -> JobSpec:
        if self.calibration and not self.input:
            raise ValueError("--calibration requires --input to also be provided")
        return self

    # ------------------------------------------------------------------ #
    #  Derived properties                                                  #
    # ------------------------------------------------------------------ #

    @property
    def sample_base(self) -> str:
        return Path(self.sample).name

    @property
    def input_base(self) -> str:
        return Path(self.input).name if self.input else ""

    @property
    def genome_base(self) -> str:
        return Path(self.genome).name

    @property
    def spikein_base(self) -> str:
        return Path(self.calibration).name if self.calibration else ""

    @property
    def do_input(self) -> bool:
        return bool(self.input)

    @property
    def do_calibration(self) -> bool:
        return bool(self.calibration)

    @property
    def do_peaks(self) -> bool:
        return self.mode in (TinyMapperMode.chip, TinyMapperMode.atac)

    @property
    def remove_duplicates_flag(self) -> str:
        """samtools markdup -r flag (empty string when keeping dups)."""
        return "" if self.keep_duplicates else "-r"

    @property
    def ignore_duplicates_flag(self) -> str:
        """bamCoverage --ignoreDuplicates flag."""
        return "" if self.keep_duplicates else "--ignoreDuplicates"

    @property
    def blacklist_options(self) -> str:
        """bamCoverage --blackListFileName option."""
        return f"--blackListFileName {self.blacklist}" if self.blacklist else ""

    @property
    def samtools_thread_opts(self) -> str:
        return f"-@ {self.threads} --output-fmt bam"

    @property
    def mnase_min_size(self) -> str:
        return self.mnase_sizes.split(",")[0]

    @property
    def mnase_max_size(self) -> str:
        return self.mnase_sizes.split(",")[1]

    @property
    def hic_base_rez(self) -> str:
        """Lowest resolution bin for HiC (first value in comma-separated list)."""
        return self.binning.split(",")[0]

    # ------------------------------------------------------------------ #
    #  Round-trip serialisation                                            #
    # ------------------------------------------------------------------ #

    def to_cli_args(self) -> list[str]:
        """Return the equivalent tinyMapper.sh long-form argument list."""
        args: list[str] = [
            "--mode",
            self.mode.value,
            "--sample",
            self.sample,
            "--genome",
            self.genome,
            "--output",
            str(self.output),
            "--threads",
            str(self.threads),
        ]
        if self.input:
            args += ["--input", self.input]
        if self.calibration:
            args += ["--calibration", self.calibration]
        if self.alignment:
            args += ["--alignment", self.alignment]
        if self.filter_opts != "-f 0x001 -f 0x002 -F 0x004 -F 0x008 -q 10":
            args += ["--filter", self.filter_opts]
        if self.blacklist:
            args += ["--blacklist", self.blacklist]
        if self.gsize != "13000000":
            args += ["--gsize", self.gsize]
        if self.keep_duplicates:
            args.append("--duplicates")
        if self.keep_intermediate:
            args.append("--keepIntermediate")
        if self.mode == TinyMapperMode.hic:
            args += ["--hicstuff", self.hicstuff_opts]
            args += ["--restriction", self.restriction]
            args += ["--binning", self.binning]
            args += ["--balance", self.balance_opts]
        if self.mode == TinyMapperMode.mnase:
            args += ["--MNaseSizes", self.mnase_sizes]
        return args
