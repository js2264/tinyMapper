"""Output file path computation following the tinyMapper naming convention.

All intermediate and final output files follow:
    <sample>^<operation>^<hash>.<ext>

This module centralises every path so mode modules never construct
path strings manually.
"""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime
from pathlib import Path


@dataclass
class RunPaths:
    """All computed file paths for a single tinyMapper run.

    Parameters
    ----------
    outdir:         root output directory (--output)
    sample_base:    basename of --sample prefix
    genome:         basename of --genome prefix
    hash_:          6-char run hash
    input_base:     basename of --input prefix (empty if no input)
    spikein:        basename of --calibration prefix (empty if no calibration)
    mnase_min:      MNase minimum fragment size
    mnase_max:      MNase maximum fragment size
    mode:           mode string (for track naming)
    hic_base_rez:   lowest HiC resolution bin
    """

    outdir: Path
    sample_base: str
    genome: str
    hash_: str
    input_base: str = ""
    spikein: str = ""
    mnase_min: str = "130"
    mnase_max: str = "200"
    mode: str = ""
    hic_base_rez: str = "500"

    # ------------------------------------------------------------------ #
    #  Log / tracking                                                      #
    # ------------------------------------------------------------------ #

    @property
    def _date(self) -> str:
        return datetime.now().strftime("%y%m%d")

    @property
    def logfile(self) -> Path:
        return self.outdir / "logs" / f"{self._date}-{self.sample_base}_{self.hash_}_log.txt"

    @property
    def cmdfile(self) -> Path:
        return self.outdir / "logs" / f"{self._date}-{self.sample_base}_{self.hash_}_commands.txt"

    @property
    def scriptfile(self) -> Path:
        return self.outdir / "logs" / f"{self._date}-{self.sample_base}_{self.hash_}_script.txt"

    @property
    def tmpfile(self) -> Path:
        return self.outdir / f"{self._date}-{self.sample_base}_{self.hash_}_INPROGRESS"

    @property
    def statfile(self) -> Path:
        return (
            self.outdir
            / "stats"
            / (
                f"sample-{self.sample_base}"
                f"_input-{self.input_base}"
                f"_genome-{self.genome}"
                f"_calibration-{self.spikein}"
                f"_{self.hash_}.counts.tsv"
            )
        )

    # ------------------------------------------------------------------ #
    #  Sample BAM files                                                    #
    # ------------------------------------------------------------------ #

    @property
    def sample_aligned_genome(self) -> Path:
        return (
            self.outdir
            / "bam/genome"
            / self.sample_base
            / f"{self.sample_base}^mapped_{self.genome}^{self.hash_}.sam"
        )

    @property
    def sample_aligned_genome_bam(self) -> Path:
        """Used by RNA mode (STAR produces BAM directly)."""
        return (
            self.outdir
            / "bam/genome"
            / self.sample_base
            / f"{self.sample_base}^mapped_{self.genome}^{self.hash_}.bam"
        )

    @property
    def sample_aligned_genome_filtered(self) -> Path:
        return (
            self.outdir
            / "bam/genome"
            / self.sample_base
            / f"{self.sample_base}^mapped_{self.genome}^filtered^{self.hash_}.bam"
        )

    @property
    def sample_non_aligned_genome(self) -> Path:
        """Prefix for unmapped reads (bowtie2 --un-conc-gz appends .1.gz/.2.gz)."""
        return (
            self.outdir
            / "fastq/genome"
            / self.sample_base
            / f"{self.sample_base}^unmapped_{self.genome}^{self.hash_}"
        )

    @property
    def sample_aligned_calibration(self) -> Path:
        return (
            self.outdir
            / "bam/spikein"
            / self.sample_base
            / f"{self.sample_base}^mapped_{self.spikein}^{self.hash_}.sam"
        )

    @property
    def sample_non_aligned_calibration(self) -> Path:
        return (
            self.outdir
            / "fastq/spikein"
            / self.sample_base
            / f"{self.sample_base}^unmapped_{self.spikein}^{self.hash_}"
        )

    @property
    def sample_non_aligned_genome_aligned_calibration(self) -> Path:
        return (
            self.outdir
            / "bam/spikein"
            / self.sample_base
            / (f"{self.sample_base}^unmapped_{self.genome}^mapped_{self.spikein}^{self.hash_}.sam")
        )

    @property
    def sample_non_aligned_calibration_aligned_genome(self) -> Path:
        return (
            self.outdir
            / "bam/genome"
            / self.sample_base
            / (f"{self.sample_base}^unmapped_{self.spikein}^mapped_{self.genome}^{self.hash_}.sam")
        )

    @property
    def sample_non_aligned_genome_aligned_calibration_filtered(self) -> Path:
        return (
            self.outdir
            / "bam/spikein"
            / self.sample_base
            / (
                f"{self.sample_base}"
                f"^unmapped_{self.genome}"
                f"^mapped_{self.spikein}"
                f"^filtered^{self.hash_}.bam"
            )
        )

    @property
    def sample_non_aligned_calibration_aligned_genome_filtered(self) -> Path:
        return (
            self.outdir
            / "bam/genome"
            / self.sample_base
            / (
                f"{self.sample_base}"
                f"^unmapped_{self.spikein}"
                f"^mapped_{self.genome}"
                f"^filtered^{self.hash_}.bam"
            )
        )

    # ------------------------------------------------------------------ #
    #  Input BAM files                                                     #
    # ------------------------------------------------------------------ #

    @property
    def input_aligned_genome(self) -> Path:
        return (
            self.outdir
            / "bam/genome"
            / self.input_base
            / f"{self.input_base}^mapped_{self.genome}^{self.hash_}.sam"
        )

    @property
    def input_aligned_genome_filtered(self) -> Path:
        return (
            self.outdir
            / "bam/genome"
            / self.input_base
            / f"{self.input_base}^mapped_{self.genome}^filtered^{self.hash_}.bam"
        )

    @property
    def input_non_aligned_genome(self) -> Path:
        return (
            self.outdir
            / "fastq/genome"
            / self.input_base
            / f"{self.input_base}^unmapped_{self.genome}^{self.hash_}"
        )

    @property
    def input_aligned_calibration(self) -> Path:
        return (
            self.outdir
            / "bam/spikein"
            / self.input_base
            / f"{self.input_base}^mapped_{self.spikein}^{self.hash_}.sam"
        )

    @property
    def input_non_aligned_calibration(self) -> Path:
        return (
            self.outdir
            / "fastq/spikein"
            / self.input_base
            / f"{self.input_base}^unmapped_{self.spikein}^{self.hash_}"
        )

    @property
    def input_non_aligned_genome_aligned_calibration(self) -> Path:
        return (
            self.outdir
            / "bam/spikein"
            / self.input_base
            / (f"{self.input_base}^unmapped_{self.genome}^mapped_{self.spikein}^{self.hash_}.sam")
        )

    @property
    def input_non_aligned_calibration_aligned_genome(self) -> Path:
        return (
            self.outdir
            / "bam/genome"
            / self.input_base
            / (f"{self.input_base}^unmapped_{self.spikein}^mapped_{self.genome}^{self.hash_}.sam")
        )

    @property
    def input_non_aligned_genome_aligned_calibration_filtered(self) -> Path:
        return (
            self.outdir
            / "bam/spikein"
            / self.input_base
            / (
                f"{self.input_base}"
                f"^unmapped_{self.genome}"
                f"^mapped_{self.spikein}"
                f"^filtered^{self.hash_}.bam"
            )
        )

    @property
    def input_non_aligned_calibration_aligned_genome_filtered(self) -> Path:
        return (
            self.outdir
            / "bam/genome"
            / self.input_base
            / (
                f"{self.input_base}"
                f"^unmapped_{self.spikein}"
                f"^mapped_{self.genome}"
                f"^filtered^{self.hash_}.bam"
            )
        )

    # ------------------------------------------------------------------ #
    #  Track files                                                         #
    # ------------------------------------------------------------------ #

    @property
    def sample_raw_track(self) -> Path:
        if self.spikein:
            suffix = f"^unmapped_{self.spikein}^mapped_{self.genome}^{self.hash_}.CPM.bw"
        elif self.mode == "RNA":
            suffix = f"^mapped_{self.genome}^{self.hash_}.unstranded.CPM.bw"
        else:
            suffix = f"^mapped_{self.genome}^{self.hash_}.CPM.bw"
        return self.outdir / "tracks" / self.sample_base / f"{self.sample_base}{suffix}"

    @property
    def input_raw_track(self) -> Path:
        if self.spikein:
            suffix = f"^unmapped_{self.spikein}^mapped_{self.genome}^{self.hash_}.CPM.bw"
        else:
            suffix = f"^mapped_{self.genome}^{self.hash_}.CPM.bw"
        return self.outdir / "tracks" / self.input_base / f"{self.input_base}{suffix}"

    @property
    def sample_input_track(self) -> Path:
        if self.spikein:
            suffix = (
                f"^unmapped_{self.spikein}^mapped_{self.genome}"
                f"^{self.hash_}.vs-{self.input_base}.bw"
            )
        else:
            suffix = f"^mapped_{self.genome}^{self.hash_}.vs-{self.input_base}.bw"
        return self.outdir / "tracks" / self.sample_base / f"{self.sample_base}{suffix}"

    @property
    def sample_spikeinscaled_track(self) -> Path:
        return (
            self.outdir
            / "tracks"
            / self.sample_base
            / (
                f"{self.sample_base}"
                f"^unmapped_{self.spikein}"
                f"^mapped_{self.genome}"
                f"^{self.hash_}.CPM.calibrated.bw"
            )
        )

    @property
    def sample_track_fwd(self) -> Path:
        return (
            self.outdir
            / "tracks"
            / self.sample_base
            / f"{self.sample_base}^mapped_{self.genome}^{self.hash_}.fwd.CPM.bw"
        )

    @property
    def sample_track_rev(self) -> Path:
        return (
            self.outdir
            / "tracks"
            / self.sample_base
            / f"{self.sample_base}^mapped_{self.genome}^{self.hash_}.rev.CPM.bw"
        )

    # ------------------------------------------------------------------ #
    #  MNase-specific                                                      #
    # ------------------------------------------------------------------ #

    @property
    def sample_aligned_genome_filtered_readsize(self) -> Path:
        return (
            self.outdir
            / "bam/genome"
            / self.sample_base
            / (
                f"{self.sample_base}"
                f"^mapped_{self.genome}"
                f"^filtered^{self.mnase_min}-{self.mnase_max}"
                f"^{self.hash_}.bam"
            )
        )

    @property
    def sample_aligned_genome_filtered_readsize40(self) -> Path:
        return (
            self.outdir
            / "bam/genome"
            / self.sample_base
            / f"{self.sample_base}^mapped_{self.genome}^filtered^40bp^{self.hash_}.bam"
        )

    @property
    def sample_readsize_track(self) -> Path:
        return (
            self.outdir
            / "tracks"
            / self.sample_base
            / (
                f"{self.sample_base}"
                f"^mapped_{self.genome}"
                f"^filtered^{self.mnase_min}-{self.mnase_max}"
                f"^{self.hash_}"
                f".{self.mnase_min}-{self.mnase_max}.CPM.bw"
            )
        )

    @property
    def sample_nucpos_track(self) -> Path:
        return (
            self.outdir
            / "tracks"
            / self.sample_base
            / (
                f"{self.sample_base}"
                f"^mapped_{self.genome}"
                f"^filtered^{self.mnase_min}-{self.mnase_max}"
                f"^{self.hash_}.nucpos.CPM.bw"
            )
        )

    @property
    def sample_nuccov_track(self) -> Path:
        return (
            self.outdir
            / "tracks"
            / self.sample_base
            / (
                f"{self.sample_base}"
                f"^mapped_{self.genome}"
                f"^filtered^{self.mnase_min}-{self.mnase_max}"
                f"^{self.hash_}.nuccov.CPM.bw"
            )
        )

    # ------------------------------------------------------------------ #
    #  HiC-specific                                                        #
    # ------------------------------------------------------------------ #

    @property
    def hic_tmp_dir(self) -> Path:
        return self.outdir / "tmp" / self.hash_

    @property
    def hic_genome_fasta(self) -> Path:
        return self.outdir / "tmp" / f"{self.sample_base}^{self.hash_}.genome.fasta"

    @property
    def sample_aligned_genome_fwd(self) -> Path:
        return (
            self.outdir
            / "bam/genome"
            / self.sample_base
            / f"{self.sample_base}^mapped_{self.genome}^{self.hash_}.fwd.bam"
        )

    @property
    def sample_aligned_genome_rev(self) -> Path:
        return (
            self.outdir
            / "bam/genome"
            / self.sample_base
            / f"{self.sample_base}^mapped_{self.genome}^{self.hash_}.rev.bam"
        )

    @property
    def sample_cool(self) -> Path:
        return (
            self.outdir
            / "matrices"
            / self.sample_base
            / f"{self.sample_base}^mapped_{self.genome}^{self.hash_}.cool"
        )

    @property
    def sample_mcool(self) -> Path:
        return (
            self.outdir
            / "matrices"
            / self.sample_base
            / f"{self.sample_base}^mapped_{self.genome}^{self.hash_}.mcool"
        )

    @property
    def sample_hic(self) -> Path:
        return (
            self.outdir
            / "matrices"
            / self.sample_base
            / f"{self.sample_base}^mapped_{self.genome}^{self.hash_}.hic"
        )

    @property
    def sample_pairs_valid_idx(self) -> Path:
        return (
            self.outdir
            / "pairs"
            / self.sample_base
            / f"{self.sample_base}^mapped_{self.genome}^{self.hash_}.valid_idx.pairs"
        )

    @property
    def sample_pairs_valid_idx_filtered(self) -> Path:
        return (
            self.outdir
            / "pairs"
            / self.sample_base
            / f"{self.sample_base}^mapped_{self.genome}^{self.hash_}.valid_idx_filtered.pairs"
        )

    @property
    def sample_pairs_valid_idx_pcrfree(self) -> Path:
        return (
            self.outdir
            / "pairs"
            / self.sample_base
            / f"{self.sample_base}^mapped_{self.genome}^{self.hash_}.valid_idx_pcrfree.pairs"
        )

    @property
    def sample_pairs_frags(self) -> Path:
        return (
            self.outdir
            / "pairs"
            / self.sample_base
            / f"{self.sample_base}^mapped_{self.genome}^{self.hash_}.frags.tsv"
        )

    @property
    def sample_pairs_dist(self) -> Path:
        return (
            self.outdir
            / "pairs"
            / self.sample_base
            / f"{self.sample_base}^mapped_{self.genome}^{self.hash_}.event_distance.pdf"
        )

    @property
    def sample_pairs_hist(self) -> Path:
        return (
            self.outdir
            / "plots"
            / self.sample_base
            / f"{self.sample_base}^mapped_{self.genome}^{self.hash_}.frags_hist.pdf"
        )

    @property
    def sample_pairs_distr(self) -> Path:
        return (
            self.outdir
            / "pairs"
            / self.sample_base
            / f"{self.sample_base}^mapped_{self.genome}^{self.hash_}.event_distribution.pdf"
        )

    @property
    def sample_pairs_law(self) -> Path:
        return (
            self.outdir
            / "plots"
            / self.sample_base
            / f"{self.sample_base}^mapped_{self.genome}^{self.hash_}.distance_law.pdf"
        )

    @property
    def sample_pairs_law_tsv(self) -> Path:
        return (
            self.outdir
            / "plots"
            / self.sample_base
            / f"{self.sample_base}^mapped_{self.genome}^{self.hash_}.distance_law.tsv"
        )
