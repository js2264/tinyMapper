J. Serizay, with contribution from H. Bordelet  
GPL-3.0

The goal of `tinyMapper.sh` is to provide a minimal -but working!- workflow to process ChIP-seq (with or without input/spikein), RNA-seq, MNase-seq and ATAC-seq data. 

Currently, this workflow only works for **paired-end** data. 

The main steps are: 

- Mapping with bowtie2 (against spikein ref. as well if needed)
- Filtering `bam` files: 
    - Fixing mates
    - Removing duplicates
    - Removing unpaired reads
    - Filtering reads based on their mapping quality (Q10)
    - For Mase: extra filtering to fragments between 70-250 bp
- Generating tracks: 
    - CPM (counts per million) tracks
    - Input and spikein-based calibrated tracks 
    - For Mase: extra track for nucleosome positions
    - For RNA-seq: directed tracks (`fwd` and `rev` transcription)
- Extract some *very* succint stats on mapping results
- Keep everything tidy, organized and documented.

**DISCLAIMER:** This is by **no means** the "best" or "only" way to process sequencing data. Do not hesitate to give suggestions / feedbacks!

### Usage 

```
Usage: tinyMapper.sh -m <MODE> -s <SAMPLE> -g <GENOME> -o <OUTPUT> [ -i <INPUT> | -c <CALIBRATION> | -t <THREADS> | -M <MEMORY> | -k <1, 0> ]

      -m | --mode                 Mapping mode ('ChIP', 'MNase', 'ATAC', 'RNA') (Default: ChIP)
      -s | --sample               Path prefix to sample *_R*.fastq.gz (e.g. for ~/reads/JS001_R*.fastq.gz files: '~/reads/JS001')
      -g | --genome               Path prefix to reference genome (e.g. for ~/genome/W303/W303.fa fasta file: '~/genome/W303/W303')
      -o | --output               Path to store results
      -i | --input                (Optional) Path prefix to input *_R*.fastq.gz
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

        `./tinyMapper.sh -m ChIP -s ~/HB44 -g ~/genomes/R64-1-1/R64-1-1 -o ~/results`
    
    - With input:

        `./tinyMapper.sh -m ChIP -s ~/HB44 -i ~/HB42 -g ~/genomes/R64-1-1/R64-1-1 -o ~/results`
    
    - With input and calibration:

        `./tinyMapper.sh -m ChIP -s ~/HB44 -i ~/HB42 -g ~/genomes/R64-1-1/R64-1-1 -c ~/genomes/Cglabrata/Cglabrata -o ~/results`
    
* **RNA-seq mode**:

    `./tinyMapper.sh -m RNA -s ~/AB4 -g ~/genomes/W303/W303 -o ~/results`

* **MNase-seq mode**:

    `./tinyMapper.sh -m MNase -s ~/CH266 -g ~/genomes/W303/W303 -o ~/results`

### Required utilities

`bowtie2`  
`samtools`  
`deeptools`  
`macs2`  
