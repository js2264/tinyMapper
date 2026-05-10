# tinyMapper

A minimalist yet versatile workflow to process ChIP-seq (with or without
input/spikein), RNA-seq, MNase-seq, ATAC-seq, Hi-C and shotgun sequencing data.
Hi-C mode delegates to [`hicstuff`](http://doi.org/10.5281/zenodo.4066363) and
`cooler`. All modes require **paired-end** data.

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
micromamba env create -n tinymapper -f tinyMapper/env/tinymapper.yaml
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
tinymapper [OPTIONS]

Required:
  -m, --mode      ChIP | MNase | ATAC | RNA | HiC | shotgun
  -s, --sample    Path prefix to sample FASTQ files
                  (e.g. ~/reads/JS001 for ~/reads/JS001_R{1,2}.fq.gz)
  -g, --genome    Path prefix to reference genome index
                  (e.g. ~/genomes/W303/W303 for bowtie2/STAR index)

Core optional:
  -o, --output         Output directory (default: results/)
  -i, --input          (ChIP) Path prefix to input/control sample
  -c, --calibration    (ChIP) Path prefix to spikein/calibration genome index
  -t, --threads        Number of CPU threads (default: 8)

Alignment / filtering:
  -a, --alignment      Extra bowtie2 options (default: '--maxins 1000')
  -f, --filter         samtools view filter options (default: '-f 2 -q 10')
  -bl, --blacklist     BED file of blacklist regions
  --gsize              Effective genome size for macs2 (default: 1.2e7)
  -d, --duplicates     Keep duplicate reads (flag, default: remove duplicates)

HiC:
  --hicstuff           Extra arguments for hicstuff (default: '--iterative --duplicates --filter --plot')
  --restriction        Restriction enzyme(s) (default: 'DpnII,HinfI')
  --binning            Matrix resolutions, comma-separated (default: '10000,20000,40000,160000,1280000')
  --balance            Balance Hi-C matrix (flag)

MNase:
  -M, --MNaseSizes     Fragment size range MIN,MAX (default: '70,250')

Output:
  -k, --keepIntermediate    Keep intermediate files
  --dry-run                 Print commands without executing
  -v, --version             Show version and exit
  -h, --help                Show this message and exit
```

FASTQ files must be paired-end. The sample prefix is matched against these
filename patterns (in order):

- `<SAMPLE>_R1.fq.gz` / `<SAMPLE>_R2.fq.gz` *(preferred)*
- `<SAMPLE>_R1.fastq.gz` / `<SAMPLE>_R2.fastq.gz`
- `<SAMPLE>_nxq_R1.fq.gz` / `<SAMPLE>_nxq_R2.fq.gz`
- `<SAMPLE>.end1.fq.gz` / `<SAMPLE>.end2.fq.gz`
- `<SAMPLE>.end1.gz` / `<SAMPLE>.end2.gz`

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

## Fetching reads

```sh
rsync -vhrP user@sftpcampus:/pasteur/gaia/projets/p02/Rsg_reads/<FOLDER>/<ID>*.fq.gz .
```

---

## Processing details

| Step | Detail |
|------|--------|
| Alignment | bowtie2 (all modes except RNA → STAR) |
| BAM filtering | fixmate, markdup, `samtools view -f 2 -q 10` |
| Duplicate removal | `samtools markdup` (skip with `--duplicates`) |
| Coverage tracks | `bamCoverage` — CPM; calibrated (ChIP+spikein); fwd/rev (RNA) |
| Peak calling | `macs2 callpeak` (ChIP, ATAC) |
| MNase fragment filter | keep fragments 70–250 bp (configurable with `--MNaseSizes`) |
| Hi-C | `hicstuff pipeline` → `cooler zoomify` → optionally balanced |
| Logs | `<output>/logs/<sample>^<op>^<hash>-log.txt` and `-commands.txt` |

---

## Acknowledgments

- A. Cournac, A. Bignaud & F. Girard for tests.
- H. Bordelet for sharing her mapping scripts and configuration.
- L. Meneu for suggestions of improvements in documentation and raising bugs.
