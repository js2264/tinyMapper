# tinyMapper 

A minimalist yet versatile workflow to process ChIP-seq (with or without input/spikein), RNA-seq, MNase-seq and ATAC-seq data.  
`tinyMapper` can also operate as a thin wrapper to process HiC data with `hicstuff` ([Cyril Matthey-Doret et al.](http://doi.org/10.5281/zenodo.4066363)) and `cooler`.  
Currently, this workflow only works for **paired-end** data. 

**DISCLAIMER:** 

- This is by **no means** the "best" or "only" way to process sequencing data. Do not hesitate to give suggestions / feedbacks to improve this workflow.
- This workflow does **NOT** include any proper QC / validation of the data. At the very least, do run `fastqc` on the sequencing data. Further QC checks are highly recommended, and will vary depending on which assay is performed. 

### Installation

`tinyMapper.sh` can be directly downloaded from `GitHub`. It is a single script 
file (`.sh`) which can be executed with: `./tinyMapper.sh`. 

```sh
wget https://raw.githubusercontent.com/js2264/tinyMapper/master/tinyMapper.sh
chmod +x tinyMapper.sh
```

Dependencies can be installed via conda. You can create and activate a conda 
environment using the yaml file we provide as follows:

```sh
conda env create -f https://raw.githubusercontent.com/js2264/tinyMapper/master/tinymapper.yaml
conda activate tinymapper
```

If HiC is going to be mapped with `hicstuff`, don't forget to install it as well! 

```sh
pip install hicstuff
```

### Usage 

Just download `tinyMapper.sh` script and use it!

```sh
wget https://raw.githubusercontent.com/js2264/tinyMapper/master/tinyMapper.sh
chmod +x tinyMapper.sh
./tinyMapper.sh -h
```

```
Usage: ./tinyMapper.sh --mode <MODE> --sample <SAMPLE> --genome <GENOME> --output <OUTPUT> [ additional arguments ]

---------------------- BASIC ARGUMENTS -----------------------------------------

   -m|--mode <MODE>                 Mapping mode (ChIP, MNase, ATAC, RNA, HiC) (Default: ChIP)
   -s|--sample <SAMPLE>             Path prefix to sample \`<SAMPLE>_R{1,2}.fastq.gz\` (e.g. for \`~/reads/JS001_R{1,2}.fastq.gz\` files, use \`--sample ~/reads/JS001\`)
   -g|--genome <GENOME>             Path prefix to reference genome (e.g. for \`~/genome/W303/W303.fa\` fasta file, use \`--genome ~/genome/W303/W303\`)
   -h|--help                        Print help ('--help' for examples)


---------------------- ADVANCED ARGUMENTS --------------------------------------

   -i|--input <INPUT>               (Optional) Path prefix to input \`<INPUT>_R{1,2}.fastq.gz\`
   -c|--calibration <CALIBRATION>   (Optional) Path prefix to genome used for calibration
   -bl|--blacklist <BED>            Bed file of blacklist regions
   -a|--alignment <ALIGN.>          Alignment options for \`bowtie2\` (between single quotes)
                                    Default: '' (no specific options)
   -f|--filter <FILTER>             Filtering options for \`samtools view\` (between single quotes)
                                    Default: '-f 2 -q 10' ('-f 2' to only keep concordant mapped and paired reads, '-q 10' to filter out reads with mapping quality score < 10)
   -d|--duplicates                  Keep duplicate reads
   -hic|--hicstuff <OPT>            Additional arguments passed to hicstuff (default: \`--iterative --duplicates --filter --plot\`)
   -r|--resolutions <#>             Resolution of final matrix file (default: '10000,20000,40000,160000,1280000')
   -re|--restriction <RE>           Restriction enzyme(s) used for HiC (default: Arima \`--restriction DpnII,HinfI\`)
   -M|--MNaseSizes <MIN,MAX>        Minimum and maximum fragment size for MNase track (default: \`--MNaseSizes 70,250\`)


---------------------- OUTPUT ARGUMENTS ----------------------------------------

   -t|--threads <THREADS>           (Optional) Number of threads (Default: 8)
   -o|--output <OUTPUT>             Path to store results (Default: \`./results/\`)
   -k|--keepIntermediate            (Optional) Keep intermediate mapping files
```

Note that fastq files *MUST* be named following this convention:
   
- **Read 1:** \*_R1.fastq.gz
- **Read 2:** \*_R2.fastq.gz

### Examples

* **ChIP-seq mode**:

    - Without input:               

        ```
        ./tinyMapper.sh \
            --mode ChIP \
            -s ~/testIP \
            -g ~/genomes/R64-1-1/R64-1-1 \
            -o ~/results
        ```
    
    - With input:

        ```
        ./tinyMapper.sh --mode ChIP \
            --sample ~/testIP \
            --input ~/testInput \
            --genome ~/genomes/R64-1-1/R64-1-1 \
            --output ~/results
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
    ./tinyMapper.sh --mode RNA -s ./testRNA -g ~/genomes/W303/W303 -o ~/results
    ```

* **MNase-seq mode**:

    ```
    ./tinyMapper.sh --mode MNase -s ./testMNase -g ~/genomes/W303/W303 -o ~/results
    ```

* **HiC mode**:

    ```
    ./tinyMapper.sh --mode HiC -s ./testHiC -g ~/genomes/W303/W303 -o ~/results --resolutions 1000,2000,8000 --restriction 'DpnII,HinfI'
    ```

### Processing details

The default steps are: 

- Mapping with bowtie2 (against spikein ref. as well if needed)
- Filtering `bam` files: 
    - Fixing mates
    - Removing duplicates (can be skipped with `--duplicates`)
    - Removing reads with mapping quality < Q10 (can be adjusted / skipped with `--filter <SAMTOOLS VIEW OPTIONS>`)
    - Removing unpaired reads (can be adjusted / skipped with `--filter <SAMTOOLS VIEW OPTIONS>`)
    - For Mase: extra filtering to keep only fragments between 70-250 bp
    - For Hi-C: process `fastq` files with `hicstuff` and binnify/balance/zoomify with `cooler`
- Generating tracks: 
    - CPM (counts per million) tracks
    - Input and spikein-based calibrated tracks 
    - For Mase: extra track for nucleosome positions
    - For RNA-seq: directed tracks (`fwd` and `rev` transcription)
- Extracting some *very* succint stats on mapping results
- Keeping everything tidy, organized, documented and reproducible. Notably, when running `tinyMapper.sh`, three files are generated: 
    - `*-log.txt`: A detailed `log` file
    - `*-commands.txt`: A list of the actual commands that were executed in the pipeline
    - `*-script.txt`: A backup copy of the `tinyMapper.sh` entire script as it was at the time of the execution

### Acknowledgments

Many thanks to H. Bordelet for sharing her mapping scripts and configuration. 
