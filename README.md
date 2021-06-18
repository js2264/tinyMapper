J. Serizay, with contribution from H. Bordelet  
GPL-3.0

# tinyMapper 

The goal of `tinyMapper.sh` is to provide a minimal -but working!- workflow to process ChIP-seq (with or without input/spikein), RNA-seq, MNase-seq and ATAC-seq data. Currently, this workflow only works for **paired-end** data. 

The main steps are: 

- Mapping with bowtie2 (against spikein ref. as well if needed)
- Filtering `bam` files: 
    - Fixing mates
    - Removing duplicates
    - Removing reads with mapping quality < Q10
    - Removing unpaired reads
    - For Mase: extra filtering to keep only fragments between 70-250 bp
- Generating tracks: 
    - CPM (counts per million) tracks
    - Input and spikein-based calibrated tracks 
    - For Mase: extra track for nucleosome positions
    - For RNA-seq: directed tracks (`fwd` and `rev` transcription)
- Extract some *very* succint stats on mapping results
- Keep everything tidy, organized, documented and reproducible.

**DISCLAIMER:** 

- This is by **no means** the "best" or "only" way to process sequencing data. Do not hesitate to give suggestions / feedbacks!
- This workflow does **NOT** include any proper QC / validation of the data. At the very least, do run `fastqc` on the sequencing data. Further QC checks are highly recommended, and will vary depending on which assay is performed. 

### Usage 

Just download `tinyMapper.sh` script and use it!

```
Usage: tinyMapper.sh -m <MODE> -s <SAMPLE> -g <GENOME> -o <OUTPUT> [ -i <INPUT> | -c <CALIBRATION> | -t <THREADS> | -M <MEMORY> | -k <1, 0> ]

      -m | --mode                 Mapping mode ('ChIP', 'MNase', 'ATAC', 'RNA') (Default: ChIP)
      -s | --sample               Path prefix to sample <SAMPLE>_R*.fastq.gz (e.g. for ~/reads/JS001_R*.fastq.gz files: --sample ~/reads/JS001)
      -g | --genome               Path prefix to reference genome (e.g. for ~/genome/W303/W303.fa fasta file: --genome ~/genome/W303/W303)
      -o | --output               Path to store results
      -i | --input                (Optional) Path prefix to input <INPUT>_R*.fastq.gz
      -c | --calibration          (Optional) Path prefix to genome used for calibration
      -t | --threads              (Optional) Number of threads (Default: 8)
      -M | --memory               (Optional) Memory in bits (Default: 12294967296, which is 12Gb)
      -k | --keepIntermediate     (Optional) Keep intermediate mapping files (Default: 1 (i.e. 'false'))
      -h | --help                 Print this message
```

Note that fastq files *MUST* be named following this convention:
   
- **Read 1:** \*_R1.fastq.gz
- **Read 2:** \*_R2.fastq.gz

### Examples

* **ChIP-seq mode**:

    - Without input:               

        ```
        ./tinyMapper.sh --mode ChIP -s ~/testIP -g ~/genomes/R64-1-1/R64-1-1 -o ~/results
        ```
    
    - With input:

        ```
        ./tinyMapper.sh --mode ChIP \
            --sample ~/testIP \
            --input ~/testInput \
            --genome ~/genomes/R64-1-1/R64-1-1 \
            --output ~/results`
        ```
    
    - With input and calibration:

        ```
        ./tinyMapper.sh --mode ChIP \
            --sample ~/testIP \
            --input ~/testInput \
            --genome ~/genomes/R64-1-1/R64-1-1 \
            --calibration ~/genomes/Cglabrata/Cglabrata \
            --output ~/results
        ```
    
* **RNA-seq mode**:

    ```
    ./tinyMapper.sh --mode RNA -s ~/testRNA -g ~/genomes/W303/W303 -o ~/results
    ```

* **MNase-seq mode**:

    ```
    ./tinyMapper.sh --mode MNase -s ~/testMNase -g ~/genomes/W303/W303 -o ~/results
    ```

### Required utilities

The dependencies can be installed as follows (provided that you are working in a dedicated `conda` env.): 

```
conda create -n tinyMapper
conda install -c conda-forge -c bioconda \
    bowtie2 \
    samtools \
    deeptools \
    macs2
```

### To do

- Allow to specify filtering options in `samtools view`
- Allow to specify whether duplicates should be removed or not
- Keep track of each command that was run to save them in the log file.
- Allow different filtering of MNase fragment sizes
