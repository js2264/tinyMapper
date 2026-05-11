# tinyMapper

A minimalist yet versatile workflow to process ChIP-seq (with or without
input/spikein), RNA-seq, MNase-seq, ATAC-seq, Hi-C and shotgun sequencing data.
Hi-C mode delegates to [`hicstuff`](http://doi.org/10.5281/zenodo.4066363) and
`cooler`.

tinyMapper supports both **paired-end** and **single-end** reads. Hi-C mode
requires paired-end data. Spikein calibration (ChIP) also requires paired-end.
For single-end MNase, fragment-size filtering is skipped and only a standard
CPM track is produced.

> **Note:** tinyMapper is a Python package that orchestrates external CLI tools
> (bowtie2, STAR, samtools, deeptools, macs2, hicstuff). It does **not**
> re-implement alignment or peak-calling.

**DISCLAIMER:**

- This is by **no means** the "best" or "only" way to process sequencing data.
  Feedback and suggestions are welcome.
- This workflow does **NOT** include QC / validation. Run `fastqc` on raw reads
  at a minimum.

---

## Installation

tinyMapper is a Python package. The recommended install creates a micromamba
environment that bundles the Python package together with all bioinformatics
tools (bowtie2, STAR, samtools, deeptools, macs2, hicstuff, cooler, bedtools).

### Recommended — full install via micromamba

Requires [`micromamba`](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html).

```sh
git clone https://github.com/js2264/tinyMapper.git
micromamba env create -y -f tinyMapper/env/tinymapper.yaml
micromamba activate tinymapper
tinymapper --help
```

### Alternative — Python package only

If all bioinformatics tools are already available in your environment:

```sh
pip install git+https://github.com/js2264/tinyMapper.git
tinymapper --help
```

---

## Invocation

After activating the environment, there are two equivalent ways to call tinyMapper:

| Command | Description |
|---------|-------------|
| `tinymapper --mode ChIP ...` | Primary Python CLI (recommended) |
| `tinyMapper.sh --mode ChIP ...` | Legacy bash wrapper — forwards all arguments verbatim to `tinymapper` |

Both accept exactly the same flags. `tinyMapper.sh` is kept for compatibility
with existing Slurm scripts and `autotinymapper`.

---

## Usage

```
 Usage: tinymapper [OPTIONS]

 tinyMapper — map and process sequencing reads.
 Modes:
   ChIP    — ChIP-seq (bowtie2 → samtools → bamCoverage → macs2)
   RNA     — RNA-seq  (STAR → samtools → bamCoverage × 3)
   ATAC    — ATAC-seq (bowtie2 → samtools → bamCoverage → macs2)
   MNase   — MNase-seq (bowtie2 → samtools → size filter → 3 tracks)
   HiC     — Hi-C     (hicstuff pipeline → cooler → mcool)
   shotgun — Shotgun  (bowtie2 single-end → samtools → bamCoverage)


 Examples:
   tinymapper -m ChIP -s ~/HB44 -g ~/genomes/R64-1-1/R64-1-1 -o ~/results
   tinymapper -m RNA  -s ~/AB4  -g ~/genomes/W303/W303 -o ~/results
   tinymapper -m HiC  -s ~/CH266 -g ~/genomes/W303/W303 --binning 1000

╭─ Required ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╮
│ *  --mode    -m  [chip|rna|atac|mnase|hic|shotgun]  Mapping mode (ChIP, MNase, ATAC, RNA, HiC, shotgun). [required]                                            │
│ *  --sample  -s  TEXT                               Path prefix to sample FASTQ files.  For ~/reads/JS001_R{1,2}.fq.gz use --sample ~/reads/JS001 [required]   │
│ *  --genome  -g  TEXT                               Path prefix to reference genome.  For ~/genome/W303/W303.fa use --genome ~/genome/W303/W303 [required]     │
╰────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╯
╭─ Core optional ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╮
│ --output       -o  PATH     Directory to store results. [default: results]                                                                                     │
│ --input        -i  TEXT     (ChIP) Path prefix to input/control sample.                                                                                        │
│ --calibration  -c  TEXT     (ChIP) Path prefix to spikein/calibration genome.                                                                                  │
│ --threads      -t  INTEGER  Number of CPU threads. [default: 8]                                                                                                │
╰────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╯
╭─ Alignment / filtering ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╮
│ --alignment   -a   TEXT  Extra options passed to bowtie2 (use single quotes). [default: --maxins 1000]                                                         │
│ --filter      -f   TEXT  Filtering options for samtools view (use single quotes). [default: -f 0x001 -f 0x002 -F 0x004 -F 0x008 -q 10]                         │
│ --blacklist   -bl  TEXT  BED file of blacklist regions for bamCoverage.                                                                                        │
│ --gsize       -gs  TEXT  Effective genome size for macs2 peak calling. [default: 13000000]                                                                     │
│ --duplicates  -d         Keep duplicate reads (default: remove duplicates).                                                                                    │
╰────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╯
╭─ HiC ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╮
│ --hicstuff     -hic  TEXT  Extra arguments passed to hicstuff pipeline. [default: --mapping iterative --duplicates --filter --plot --no-cleanup]               │
│ --restriction  -re   TEXT  Restriction enzyme(s) for HiC (e.g. DpnII,HinfI). [default: HpaII,HinfI]                                                            │
│ --binning      -b    TEXT  Minimum bin resolution for HiC matrix (bp); comma-separated for multi-res. [default: 500]                                           │
│ --balance      -ba   TEXT  Balancing options for cooler zoomify. [default: --cis-only --min-nnz 3 --mad-max 7]                                                 │
╰────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╯
╭─ MNase ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╮
│ --MNaseSizes  -M  TEXT  Min,Max fragment size for MNase track. [default: 130,200]                                                                              │
╰────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╯
╭─ Output ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╮
│ --keepIntermediate  -k  Keep intermediate SAM / unmapped FASTQ files.                                                                                          │
│ --dry-run               Log commands without executing them.                                                                                                   │
│ --help              -h  Show this message and exit.                                                                                                            │
│ --version           -v  Show the version and exit.                                                                                                             │
╰────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╯
```

FASTQ files are detected automatically from the sample prefix. tinyMapper
tries paired-end patterns first, then falls back to single-end:

**Paired-end patterns** (both R1 and R2 must exist):
- `<SAMPLE>_R1.fq.gz` / `<SAMPLE>_R2.fq.gz` *(preferred)*
- `<SAMPLE>_R1.fastq.gz` / `<SAMPLE>_R2.fastq.gz`
- `<SAMPLE>_nxq_R1.fq.gz` / `<SAMPLE>_nxq_R2.fq.gz`
- `<SAMPLE>.end1.fq.gz` / `<SAMPLE>.end2.fq.gz`
- `<SAMPLE>.end1.gz` / `<SAMPLE>.end2.gz`
- Illumina `<SAMPLE>_S##_R1_*.gz` / `<SAMPLE>_S##_R2_*.gz`

**Single-end fallback** (R2 not found — only R1 required):
- `<SAMPLE>_R1.fq.gz`
- `<SAMPLE>_R1.fastq.gz`
- `<SAMPLE>_nxq_R1.fq.gz`
- `<SAMPLE>.fq.gz`
- `<SAMPLE>.fastq.gz`

| Mode | SE support | Notes |
|------|-----------|-------|
| ChIP | Yes | input control supported; spikein calibration requires PE |
| RNA | Yes | forward/reverse strand tracks still produced |
| ATAC | Yes | peaks called with `--format BAM` instead of `BAMPE` |
| MNase | Yes | fragment-size filter and nucleosome tracks skipped; CPM track only |
| shotgun | Yes | always single-end (R1+R2 concatenated as -U if both present) |
| HiC | **No** | paired-end required |

---

## Examples

### ChIP-seq

```sh
# Sample only (no input, no calibration)
tinymapper -m ChIP \
    -s ~/reads/JS001 \
    -g ~/genomes/R64-1-1/R64-1-1 \
    -o ~/results

# With input control
tinymapper -m ChIP \
    --sample ~/reads/JS001_IP \
    --input  ~/reads/JS001_input \
    --genome ~/genomes/R64-1-1/R64-1-1 \
    --output ~/results

# With input and spikein calibration
tinymapper -m ChIP \
    --sample      ~/reads/JS001_IP \
    --input       ~/reads/JS001_input \
    --genome      ~/genomes/R64-1-1/R64-1-1 \
    --calibration ~/genomes/Cglabrata/Cglabrata \
    --output      ~/results
```

### RNA-seq

```sh
tinymapper -m RNA -s ~/reads/JS001 -g ~/genomes/W303/W303 -o ~/results
```

### MNase-seq

```sh
tinymapper -m MNase -s ~/reads/JS001 -g ~/genomes/W303/W303 -o ~/results \
    --MNaseSizes 70,250
```

### ATAC-seq

```sh
tinymapper -m ATAC -s ~/reads/JS001 -g ~/genomes/W303/W303 -o ~/results
```

### Hi-C

```sh
tinymapper -m HiC \
    -s ~/reads/JS001 \
    -g ~/genomes/W303/W303 \
    -o ~/results \
    --binning 1000,2000,8000 \
    --restriction 'DpnII,HinfI'
```

### Shotgun

```sh
tinymapper -m shotgun -s ~/reads/JS001 -g ~/genomes/W303/W303 -o ~/results
```

---

## Output layout

Results are written under `--output` with the following structure:

```
<output>/
  bam/genome/          filtered BAM files (genome)
  bam/spikein/         filtered BAM files (spikein, ChIP only)
  tracks/              BigWig coverage tracks (CPM, calibrated, fwd/rev for RNA)
  peaks/               MACS2 peak files (ChIP, ATAC)
  pairs/               contact pairs (Hi-C only)
  matrices/            .cool matrices (Hi-C only)
  logs/                per-run log and command files
  tmp/                 temporary files (removed on success unless --keepIntermediate)
```

Files follow the naming convention `<sample>^<operation>^<hash>.<ext>` where
`<hash>` is a 6-character alphanumeric string unique to each run.

---

## Running on a Slurm cluster (e.g. Maestro)

Activate the environment and submit with `sbatch`:

```sh
micromamba activate tinymapper

# Generic
sbatch --mem 40G -c 10 --wrap \
    "tinymapper --mode ChIP --sample <SAMPLE> --genome <GENOME> --output <OUTPUT> --threads 8"

# ChIP examples
sbatch --mem 40G -c 10 --wrap \
    "tinymapper -m ChIP -s ~/reads/JS001_IP -g ~/genomes/S288c/S288c --threads 8"
sbatch --mem 40G -c 10 --wrap \
    "tinymapper -m ChIP -s ~/reads/JS001_IP -i ~/reads/JS001_input -g ~/genomes/S288c/S288c --threads 8"

# RNA
sbatch --mem 40G -c 10 --wrap \
    "tinymapper -m RNA -s ~/reads/JS001 -g ~/genomes/S288c/S288c --threads 8"

# Hi-C
sbatch --mem 40G -c 10 --wrap \
    "tinymapper -m HiC -s ~/reads/JS001 -g ~/genomes/S288c/S288c --threads 8"
```

`tinyMapper.sh` can be used as a drop-in replacement for the legacy command
surface (e.g. from `autotinymapper` Slurm scripts):

```sh
sbatch --mem 40G -c 10 --wrap \
    "tinyMapper.sh -m ChIP -s ~/reads/JS001_IP -g ~/genomes/S288c/S288c --threads 8"
```

---

## Acknowledgments

- A. Cournac, A. Bignaud & F. Girard for tests.
- H. Bordelet for sharing her mapping scripts and configuration.
- L. Meneu for suggestions of improvements in documentation and raising bugs.
