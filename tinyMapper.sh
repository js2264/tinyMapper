#!/bin/bash

VERSION=0.9.111

INVOC=$(printf %q "$BASH_SOURCE")$((($#)) && printf ' %q' "$@")
HASH=`LC_CTYPE=C tr -dc 'A-Z0-9' < /dev/urandom | head -c 6`

## ------------------------------------------------------------------
## -------- HELPER FUNCTIONS ----------------------------------------
## ------------------------------------------------------------------

function usage() {
    echo -e ""
    echo -e "tinyMapper v${VERSION}"
    echo -e "GPL-3.0"
    echo -e ""
    echo -e "Usage: tinyMapper.sh --mode <MODE> --sample <SAMPLE> --genome <GENOME> --output <OUTPUT> [ additional arguments ]"
    echo -e ""
    echo -e "================================================================================"
    echo -e ""
    echo -e "---------------------- BASIC ARGUMENTS -----------------------------------------"
    echo -e ""
    echo -e "   -m|--mode <MODE>                 Mapping mode (ChIP, MNase, ATAC, RNA, HiC) (Default: ChIP)"
    echo -e "   -s|--sample <SAMPLE>             Path prefix to sample \`<SAMPLE>_R{1,2}.fastq.gz\` (e.g. for \`~/reads/JS001_R{1,2}.fastq.gz\` files, use \`--sample ~/reads/JS001\`)"
    echo -e "   -g|--genome <GENOME>             Path prefix to reference genome (e.g. for \`~/genome/W303/W303.fa\` fasta file, use \`--genome ~/genome/W303/W303\`)"
    echo -e "   -h|--help                        Print help ('--help' for examples)"
    echo -e ""
    echo -e ""
    echo -e "---------------------- ADVANCED ARGUMENTS --------------------------------------"
    echo -e ""
    echo -e "   -i|--input <INPUT>               (Optional) Path prefix to input \`<INPUT>_R{1,2}.fastq.gz\`"
    echo -e "   -c|--calibration <CALIBRATION>   (Optional) Path prefix to genome used for calibration"
    echo -e "   -bl|--blacklist <BED>            Bed file of blacklist regions"
    echo -e "   -a|--alignment <ALIGN.>          Alignment options for \`bowtie2\` (between single quotes)"
    echo -e "                                    Default: '' (no specific options)"
    echo -e "   -f|--filter <FILTER>             Filtering options for \`samtools view\` (between single quotes)"
    echo -e "                                    Default: '-f 2 -q 10' ('-f 2' to only keep concordant mapped and paired reads, '-q 10' to filter out reads with mapping quality score < 10)"
    echo -e "   -d|--duplicates                  Keep duplicate reads"
    echo -e ""
    echo -e "   -hic|--hicstuff <OPT>            Additional arguments passed to hicstuff (default: \`--iterative --duplicates --filter --plot\`)"
    echo -e "   -r|--resolutions <#>             Resolution of final matrix file (default: '10000,20000,40000,160000,1280000')"
    echo -e "   -re|--restriction <RE>           Restriction enzyme(s) used for HiC (default: Arima \`--restriction DpnII,HinfI\`)"
    echo -e ""
    echo -e "   -M|--MNaseSizes <MIN,MAX>        Minimum and maximum fragment size for MNase track (default: \`--MNaseSizes 70,250\`)"
    echo -e ""
    echo -e ""
    echo -e "---------------------- OUTPUT ARGUMENTS ----------------------------------------"
    echo -e ""
    echo -e "   -t|--threads <THREADS>           (Optional) Number of threads (Default: 8)"
    echo -e "   -o|--output <OUTPUT>             Path to store results (Default: \`./results/\`)"
    echo -e "   -k|--keepIntermediate            (Optional) Keep intermediate mapping files"
    echo -e ""
}

function usage_extended() {
    usage
    echo -e "================================================================================"
    echo -e ""
    echo -e "---------------------- EXAMPLES ------------------------------------------------"
    echo -e ""
    echo -e "   ChIP-seq mode:"
    echo -e ""
    echo -e "      Without input:               ./tinyMapper.sh -m ChIP -s ~/HB44 -g ~/genomes/R64-1-1/R64-1-1 -o ~/results"
    echo -e "      With input:                  ./tinyMapper.sh -m ChIP -s ~/HB44 -i ~/HB42 -g ~/genomes/R64-1-1/R64-1-1 -o ~/results"
    echo -e "      With input and calibration:  ./tinyMapper.sh -m ChIP -s ~/HB44 -i ~/HB42 -g ~/genomes/R64-1-1/R64-1-1 -c ~/genomes/Cglabrata/Cglabrata -o ~/results"
    echo -e ""
    echo -e "   RNA-seq mode:"
    echo -e ""
    echo -e "      ./tinyMapper.sh -m RNA -s ~/AB4 -g ~/genomes/W303/W303 -o ~/results"
    echo -e ""
    echo -e "   MNase-seq mode:"
    echo -e ""
    echo -e "      ./tinyMapper.sh -m MNase -s ~/CH266 -g ~/genomes/W303/W303 -o ~/results"
    echo -e ""
    echo -e "   HiC mode (through hicstuff):"
    echo -e ""
    echo -e "      ./tinyMapper.sh -m HiC -s ~/CH266 -g ~/genomes/W303/W303 -o ~/results --resolutions 1000,2000,4000 --restriction 'DpnII,HinfI'"
    echo -e ""
    echo -e "================================================================================"
    echo -e ""
    echo -e "---------------------- REQUIRED UTILITIES --------------------------------------"
    echo -e ""
    echo -e "   bowtie2"
    echo -e "   samtools"
    echo -e "   deeptools"
    echo -e "   macs2 (for ChIP/ATAC)"
    echo -e "   hicstuff (for HiC)"
    echo -e "   cooler (for HiC)"
    echo -e ""
    echo -e "================================================================================"
    echo -e ""
    echo -e "tinyMapper v${VERSION}"
    echo -e ""
}

function is_set() { 
    test -n "${1}" ; echo $?
}

function fn_log {
    date=`date "+%y-%m-%d %H:%M:%S"`
    BOLD="\e[1m"
    GREEN="\e[32m"
    RED="\e[31m"
    BLUE="\e[96m"
    YELLOW="\e[33m"
    RESET="\e[0m"
    echo -e "${BOLD}${BLUE}${date} | ${GREEN}[INFO]${RESET} $@"
}

function fn_error {
    date=`date "+%y-%m-%d %H:%M:%S"`
    BOLD="\e[1m"
    GREEN="\e[32m"
    RED="\e[31m"
    BLUE="\e[96m"
    YELLOW="\e[33m"
    RESET="\e[0m"
    echo -e "${BOLD}${BLUE}${date} | ${RED}[ERR.]${RESET} $@"
}

function fn_warning {
    date=`date "+%y-%m-%d %H:%M:%S"`
    BOLD="\e[1m"
    GREEN="\e[32m"
    RED="\e[31m"
    MAGENTA="\e[35m"
    BLUE="\e[96m"
    YELLOW="\e[33m"
    RESET="\e[0m"
    echo -e "${BOLD}${BLUE}${date} | ${MAGENTA}[WAR.]${RESET} $@"
}

function fn_exec {
    date=`date "+%y-%m-%d %H:%M:%S"`
    BOLD="\e[1m"
    GREEN="\e[32m"
    RED="\e[31m"
    BLUE="\e[96m"
    YELLOW="\e[33m"
    RESET="\e[0m"
    cmd=`echo $1 | tr -s '' | sed 's,2>>.*,,' `
    if test `is_set $2` == 0 ; then 
        echo -e "${BOLD}${BLUE}${date} | ${YELLOW}[EXEC]${RESET} ${cmd}" >> $2
    fi
    eval ${cmd}
}

function fastqfastcnt {

    fix_base_count() {
        local counts=($(cat))
        echo "${counts[0]}"
    }

    gzip -dc $1 \
        | awk 'NR % 4 == 2' \
        | wc -l \
        | fix_base_count

}

## ------------------------------------------------------------------
## -------- PARSING ARGUMENTS ---------------------------------------
## ------------------------------------------------------------------

# Default values of arguments

MODE='ChIP'
SAMPLE=''
INPUT=''
GENOME_=''
SPIKEIN_=''
OUTDIR='results'
BOWTIEOPTIONS=''
FILTEROPTIONS='-f 2 -q 10'
HICSTUFFOPTIONS=' --mapping iterative --duplicates --filter --plot --no-cleanup'
HICREZ='10000,20000,40000,160000,1280000'
MNASESIZES='70,250'
RE=' DpnII,HinfI '
BLACKLISTBEDFILE=''
KEEPDUPLICATES=1
KEEPFILES=1
CPU=8

# Custom values for test

# MODE=ChIP
# SAMPLE=sftpcampus:tmp/HB44
# INPUT=sftpcampus:tmp/HB42
# GENOME_=~/genomes/R64-1-1/R64-1-1
# SPIKEIN_=''
# OUTDIR=results
# FILTEROPTIONS='-F 4 -q 5'
# KEEPDUPLICATES=0
# CPU=16
# KEEPFILES=0

# MODE=HiC
# SAMPLE=test
# GENOME_=~/genomes/W303/W303
# HICREZ='1000,2000,4000,8000'
# KEEPFILES=0

if test `is_set "${1}"` == 1 ; then
    usage && exit 0
fi

for arg in "$@"
do
    case $arg in
        #####
        ##### BASIC ARGUMENTS
        #####
        -m|--mode)
        MODE="${2}"
        shift 
        shift 
        ;;
        -s|--sample)
        SAMPLE="${2}"
        shift 
        shift 
        ;;
        -i|--input)
        INPUT="${2}"
        shift 
        shift 
        ;;
        -g|--genome)
        GENOME_="${2}"
        shift 
        shift 
        ;;
        -c|--calibration)
        SPIKEIN_="${2}"
        shift 
        shift 
        ;;
        -o|--output)
        OUTDIR="${2}"
        shift 
        shift 
        ;;
        #####
        ##### ADVANCED ARGUMENTS
        #####
        -a|--alignment)
        BOWTIEOPTIONS=${2}
        shift 
        shift 
        ;;
        -f|--filter)
        FILTEROPTIONS=${2}
        shift 
        shift 
        ;;
        -bl|--blacklist)
        BLACKLISTBEDFILE=${2}
        shift 
        shift 
        ;;
        -hic|--hicstuff)
        HICSTUFFOPTIONS=${2}
        shift 
        shift 
        ;;
        -re|--restriction)
        RE=${2}
        shift 
        shift 
        ;;
        -r|--resolutions)
        HICREZ=${2}
        shift 
        shift 
        ;;
        -M|--MNaseSizes)
        MNASESIZES=${2}
        shift 
        shift 
        ;;
        -d|--duplicates)
        KEEPDUPLICATES=0
        shift 
        ;;
        #####
        ##### OPTIONAL ARGUMENTS
        #####
        -t|--threads)
        CPU="${2}"
        shift 
        shift 
        ;;
        -k|--keepIntermediate)
        KEEPFILES=0
        shift 
        ;;
        -h)
        usage && exit 0
        ;;
        --help)
        usage_extended && exit 0
        ;;
        -v|--version)
        echo -e "tinyMapper v${VERSION}" && exit 0
        ;;
    esac
done

## ------------------------------------------------------------------
## -------- PREPARING VARIABLES -------------------------------------
## ------------------------------------------------------------------

# - Genome(s) variable
GENOME_DIR=`dirname "${GENOME_}"`
GENOME=`basename "${GENOME_}"`
GENOME_BASE="${GENOME_DIR}"/"${GENOME}"
GENOME_FA="${GENOME_BASE}.fa"
SPIKEIN_DIR=`dirname "${SPIKEIN_}"`
SPIKEIN=`basename "${SPIKEIN_}"`
SPIKEIN_BASE="${SPIKEIN_DIR}"/"${SPIKEIN}"
SPIKEIN_FA="${SPIKEIN_BASE}.fa"

# - Sample (input) variable
SAMPLE_DIR=`dirname "${SAMPLE}"`
SAMPLE_BASE=`basename "${SAMPLE}"`
SAMPLE_R1="${SAMPLE}_R1.fastq.gz"
SAMPLE_R2="${SAMPLE}_R2.fastq.gz"
INPUT_DIR=`dirname "${INPUT}"`
INPUT_BASE=`basename "${INPUT}"`
INPUT_R1="${INPUT}_R1.fastq.gz"
INPUT_R2="${INPUT}_R2.fastq.gz"

# - Specific files
STATFILE="${OUTDIR}/stats/sample-${SAMPLE_BASE}_input-${INPUT_BASE}_genome-${GENOME}_calibration-${SPIKEIN}_${HASH}".counts.tsv
LOGFILE="${OUTDIR}/logs/`date "+%y%m%d"`-${HASH}-log.txt"
CMDFILE="${OUTDIR}/logs/`date "+%y%m%d"`-${HASH}-commands.txt"
SCRIPTFILE="${OUTDIR}/logs/`date "+%y%m%d"`-${HASH}-script.txt"
TMPFILE="${OUTDIR}/`date "+%y%m%d"`-${HASH}-INPROGRESS"

# - Advanced options
SAMTOOLS_OPTIONS=" -@ ${CPU} --output-fmt bam "
DO_INPUT=`is_set "${INPUT}"`
DO_CALIBRATION=`is_set "${SPIKEIN}"`
DO_PEAKS=`if test "${MODE}" == 'ChIP' || test "${MODE}" == 'ATAC'; then echo 0; else echo 1; fi`
REMOVE_DUPLICATES=`if test "${KEEPDUPLICATES}" == 1 ; then echo " -r " ; else echo " " ; fi`
BLACKLIST_OPTIONS=`if test $(is_set "${BLACKLISTBEDFILE}") == 0 ; then echo " --blackListFileName ${BLACKLISTBEDFILE} " ; else echo " " ; fi`
FIRSTREZ=`echo "${HICREZ}" | sed 's/,.*//' | sed 's,000$,kb,'`
MNASE_MINSIZE=`echo "${MNASESIZES}" | sed 's/,.*//'`
MNASE_MAXSIZE=`echo "${MNASESIZES}" | sed 's/.*,//'`

# - Bam file names
SAMPLE_ALIGNED_GENOME="${OUTDIR}"/bam/genome/"${SAMPLE_BASE}"/"${SAMPLE_BASE}"^mapped_"${GENOME}"^"${HASH}".sam
SAMPLE_ALIGNED_GENOME_FILTERED="${OUTDIR}"/bam/genome/"${SAMPLE_BASE}"/"${SAMPLE_BASE}"^mapped_"${GENOME}"^filtered^"${HASH}".bam
SAMPLE_NON_ALIGNED_GENOME="${OUTDIR}"/fastq/genome/"${SAMPLE_BASE}"/"${SAMPLE_BASE}"^unmapped_"${GENOME}"^"${HASH}"
SAMPLE_ALIGNED_CALIBRATION="${OUTDIR}"/bam/spikein/"${SAMPLE_BASE}"/"${SAMPLE_BASE}"^mapped_"${SPIKEIN}"^"${HASH}".sam
SAMPLE_NON_ALIGNED_CALIBRATION="${OUTDIR}"/fastq/spikein/"${SAMPLE_BASE}"/"${SAMPLE_BASE}"^unmapped_"${SPIKEIN}"^"${HASH}"
SAMPLE_NON_ALIGNED_GENOME_ALIGNED_CALIBRATION="${OUTDIR}"/bam/spikein/"${SAMPLE_BASE}"/"${SAMPLE_BASE}"^unmapped_"${GENOME}"^mapped_"${SPIKEIN}"^"${HASH}".sam
SAMPLE_NON_ALIGNED_CALIBRATION_ALIGNED_GENOME="${OUTDIR}"/bam/genome/"${SAMPLE_BASE}"/"${SAMPLE_BASE}"^unmapped_"${SPIKEIN}"^mapped_"${GENOME}"^"${HASH}".sam
SAMPLE_NON_ALIGNED_GENOME_ALIGNED_CALIBRATION_FILTERED="${OUTDIR}"/bam/spikein/"${SAMPLE_BASE}"/"${SAMPLE_BASE}"^unmapped_"${GENOME}"^mapped_"${SPIKEIN}"^filtered^"${HASH}".bam
SAMPLE_NON_ALIGNED_CALIBRATION_ALIGNED_GENOME_FILTERED="${OUTDIR}"/bam/genome/"${SAMPLE_BASE}"/"${SAMPLE_BASE}"^unmapped_"${SPIKEIN}"^mapped_"${GENOME}"^filtered^"${HASH}".bam
INPUT_ALIGNED_GENOME="${OUTDIR}"/bam/genome/"${INPUT_BASE}"/"${INPUT_BASE}"^mapped_"${GENOME}"^"${HASH}".sam
INPUT_ALIGNED_GENOME_FILTERED="${OUTDIR}"/bam/genome/"${INPUT_BASE}"/"${INPUT_BASE}"^mapped_"${GENOME}"^filtered^"${HASH}".bam
INPUT_NON_ALIGNED_GENOME="${OUTDIR}"/fastq/genome/"${INPUT_BASE}"/"${INPUT_BASE}"^unmapped_"${GENOME}"^"${HASH}"
INPUT_ALIGNED_CALIBRATION="${OUTDIR}"/bam/spikein/"${INPUT_BASE}"/"${INPUT_BASE}"^mapped_"${SPIKEIN}"^"${HASH}".sam
INPUT_NON_ALIGNED_CALIBRATION="${OUTDIR}"/fastq/spikein/"${INPUT_BASE}"/"${INPUT_BASE}"^unmapped_"${SPIKEIN}"^"${HASH}"
INPUT_NON_ALIGNED_GENOME_ALIGNED_CALIBRATION="${OUTDIR}"/bam/spikein/"${INPUT_BASE}"/"${INPUT_BASE}"^unmapped_"${GENOME}"^mapped_"${SPIKEIN}"^"${HASH}".sam
INPUT_NON_ALIGNED_CALIBRATION_ALIGNED_GENOME="${OUTDIR}"/bam/genome/"${INPUT_BASE}"/"${INPUT_BASE}"^unmapped_"${SPIKEIN}"^mapped_"${GENOME}"^"${HASH}".sam
INPUT_NON_ALIGNED_GENOME_ALIGNED_CALIBRATION_FILTERED="${OUTDIR}"/bam/spikein/"${INPUT_BASE}"/"${INPUT_BASE}"^unmapped_"${GENOME}"^mapped_"${SPIKEIN}"^filtered^"${HASH}".bam
INPUT_NON_ALIGNED_CALIBRATION_ALIGNED_GENOME_FILTERED="${OUTDIR}"/bam/genome/"${INPUT_BASE}"/"${INPUT_BASE}"^unmapped_"${SPIKEIN}"^mapped_"${GENOME}"^filtered^"${HASH}".bam

# - Track file names
SAMPLE_RAW_TRACK="${OUTDIR}"/tracks/"${SAMPLE_BASE}"/"${SAMPLE_BASE}"^mapped_"${GENOME}"^"${HASH}".CPM.bw
INPUT_RAW_TRACK="${OUTDIR}"/tracks/"${INPUT_BASE}"/"${INPUT_BASE}"^mapped_"${GENOME}"^"${HASH}".CPM.bw
SAMPLE_INPUT_TRACK="${OUTDIR}"/tracks/"${SAMPLE_BASE}"/"${SAMPLE_BASE}"^mapped_"${GENOME}"^"${HASH}".vs-"${INPUT_BASE}".bw

if test "${DO_CALIBRATION}" == 0 ; then
    SAMPLE_RAW_TRACK="${OUTDIR}"/tracks/"${SAMPLE_BASE}"/"${SAMPLE_BASE}"^unmapped_"${SPIKEIN}"^mapped_"${GENOME}"^"${HASH}".CPM.bw
    INPUT_RAW_TRACK="${OUTDIR}"/tracks/"${INPUT_BASE}"/"${INPUT_BASE}"^unmapped_"${SPIKEIN}"^mapped_"${GENOME}"^"${HASH}".CPM.bw
    SAMPLE_INPUT_TRACK="${OUTDIR}"/tracks/"${SAMPLE_BASE}"/"${SAMPLE_BASE}"^unmapped_"${SPIKEIN}"^mapped_"${GENOME}"^"${HASH}".vs-"${INPUT_BASE}".bw
    SAMPLE_SPIKEINSCALED_TRACK="${OUTDIR}"/tracks/"${SAMPLE_BASE}"/"${SAMPLE_BASE}"^unmapped_"${SPIKEIN}"^mapped_"${GENOME}"^"${HASH}".CPM.calibrated.bw
fi

if test "${MODE}" == MNase ; then 
    SAMPLE_ALIGNED_GENOME_FILTERED_READSIZE="${OUTDIR}"/bam/genome/"${SAMPLE_BASE}"/"${SAMPLE_BASE}"^mapped_"${GENOME}"^filtered^"${MNASE_MINSIZE}"-"${MNASE_MAXSIZE}"^"${HASH}".bam
    SAMPLE_READSIZE_TRACK="${OUTDIR}"/tracks/"${SAMPLE_BASE}"/"${SAMPLE_BASE}"^mapped_"${GENOME}"^filtered^"${MNASE_MINSIZE}"-"${MNASE_MAXSIZE}"^"${HASH}"."${MNASE_MINSIZE}"-"${MNASE_MAXSIZE}".CPM.bw
    SAMPLE_NUCPOS_TRACK="${OUTDIR}"/tracks/"${SAMPLE_BASE}"/"${SAMPLE_BASE}"^mapped_"${GENOME}"^filtered^"${MNASE_MINSIZE}"-"${MNASE_MAXSIZE}"^"${HASH}".nucpos.CPM.bw
fi

if test "${MODE}" == RNA ; then 
    SAMPLE_RAW_TRACK="${OUTDIR}"/tracks/"${SAMPLE_BASE}"/"${SAMPLE_BASE}"^mapped_"${GENOME}"^"${HASH}".unstranded.CPM.bw
    SAMPLE_TRACK_FWD="${OUTDIR}"/tracks/"${SAMPLE_BASE}"/"${SAMPLE_BASE}"^mapped_"${GENOME}"^"${HASH}".fwd.CPM.bw
    SAMPLE_TRACK_REV="${OUTDIR}"/tracks/"${SAMPLE_BASE}"/"${SAMPLE_BASE}"^mapped_"${GENOME}"^"${HASH}".rev.CPM.bw
fi

# - HiC-related file names
if test "${MODE}" == HiC ; then 
    SAMPLE_ALIGNED_GENOME_FWD="${OUTDIR}"/bam/genome/"${SAMPLE_BASE}"/"${SAMPLE_BASE}"^mapped_"${GENOME}"^"${HASH}".fwd.bam
    SAMPLE_ALIGNED_GENOME_REV="${OUTDIR}"/bam/genome/"${SAMPLE_BASE}"/"${SAMPLE_BASE}"^mapped_"${GENOME}"^"${HASH}".rev.bam
    SAMPLE_COOL="${OUTDIR}"/matrices/"${SAMPLE_BASE}"/"${SAMPLE_BASE}"^"${HASH}".cool
    SAMPLE_MCOOL="${OUTDIR}"/matrices/"${SAMPLE_BASE}"/"${SAMPLE_BASE}"^"${HASH}".mcool
fi

mkdir -p "${OUTDIR}"/logs
touch "${LOGFILE}"
touch "${TMPFILE}"

## ------------------------------------------------------------------
## -------- CHECKING THAT ALL REQUIRED FILES EXIST ------------------
## ------------------------------------------------------------------

# Abort if trying to calibrate without input
if test "${DO_INPUT}" == 1 && test "${DO_CALIBRATION}" == 0 ; then
    fn_error "Calibration can only be done if an input is provided." 2>&1 | tee -a "${LOGFILE}"
    fn_error "Aborting now." 2>&1 | tee -a "${LOGFILE}"
    rm --force "${LOGFILE}"
    exit 1
fi

# Check if a sample is provided
if test `is_set "${SAMPLE_BASE}"` == 1 ; then
    fn_error "Please provide a sample." 2>&1 | tee -a "${LOGFILE}"
    fn_error "Aborting now." 2>&1 | tee -a "${LOGFILE}"
    rm --force "${LOGFILE}"
    exit 1
fi

# If sample files are accessed through ssh, download them first
if [[ "${SAMPLE_DIR}" == *:* ]] ; then 
    echo -e "Fetching sample files from remote to \`./data/reads\`."
    p data/reads
    scp "${SAMPLE}"* data/reads/
    SAMPLE="data/reads/`basename ${SAMPLE}`"
    SAMPLE_DIR=`dirname "${SAMPLE}"`
    SAMPLE_BASE=`basename "${SAMPLE}"`
    SAMPLE_R1="${SAMPLE}_R1.fastq.gz"
    SAMPLE_R2="${SAMPLE}_R2.fastq.gz"
fi

# If input files are accessed through ssh, download them first
if [[ "${INPUT_DIR}" == *:* ]] ; then 
    echo -e "Fetching input files from remote \`./data/reads\`."
    mkdir -p data/reads
    scp "${INPUT}"* data/reads/
    INPUT="data/reads/`basename ${INPUT}`"
    INPUT_DIR=`dirname "${INPUT}"`
    INPUT_BASE=`basename "${INPUT}"`
    INPUT_R1="${INPUT}_R1.fastq.gz"
    INPUT_R2="${INPUT}_R2.fastq.gz"
fi

# If several sample files are found, zcat all of them
filestomerge=`ls "${SAMPLE}"*1*gz 2>/dev/null | grep -v "${SAMPLE}_R1.fastq.gz" | grep -v '.*_R2_.*' | grep -v '.*_2_.*' | grep -v '.*end2.*'`
cnt=`echo ${filestomerge} | tr ' ' '\n' | wc -l`
if [ "$cnt" -gt "1" ] ; then 
    echo -e "Merging several sample files together (R1):  `echo -n ${filestomerge}`"
    zcat `ls "${SAMPLE}"*1*gz 2>/dev/null | grep -v "${SAMPLE}_R1.fastq.gz"` | gzip > "${SAMPLE}_R1.fastq.gz"
    SAMPLE_R1="${SAMPLE}_R1.fastq.gz"
fi
filestomerge=`ls "${SAMPLE}"*2*gz 2>/dev/null | grep -v "${SAMPLE}_R2.fastq.gz" | grep -v '.*_R1_.*' | grep -v '.*_1_.*' | grep -v '.*end1.*'`
cnt=`echo ${filestomerge} | tr ' ' '\n' | wc -l`
if [ "$cnt" -gt "1" ] ; then 
    echo -e "Merging several sample files together (R2):  `echo -n ${filestomerge}`"
    zcat `ls "${SAMPLE}"*2*gz 2>/dev/null | grep -v "${SAMPLE}_R2.fastq.gz"` | gzip > "${SAMPLE}_R2.fastq.gz"
    SAMPLE_R2="${SAMPLE}_R2.fastq.gz"
fi

# If several input files are found, zcat all of them
filestomerge=`ls "${INPUT}"*1*gz 2>/dev/null | grep -v "${INPUT}_R1.fastq.gz" | grep -v '.*_R2_.*' | grep -v '.*_2_.*' | grep -v '.*end2.*'`
cnt=`echo ${filestomerge} | tr ' ' '\n' | wc -l`
if test "$cnt" -gt "1" && test "${DO_INPUT}" ==  0 ; then 
    echo -e "Merging several input files together (R1):  `echo -n ${filestomerge}`"
    zcat "${INPUT}"*R1*gz | gzip > "${INPUT}_R1.fastq.gz"
    INPUT_R1="${INPUT}_R1.fastq.gz"
fi
filestomerge=`ls "${INPUT}"*2*gz 2>/dev/null | grep -v "${INPUT}_R2.fastq.gz" | grep -v '.*_R1_.*' | grep -v '.*_1_.*' | grep -v '.*end1.*'`
cnt=`echo ${filestomerge} | tr ' ' '\n' | wc -l`
if test "$cnt" -gt "1" && test "${DO_INPUT}" ==  0 ; then 
    echo -e "Merging several input files together (R2):  `echo -n ${filestomerge}`"
    zcat "${INPUT}"*R2*gz | gzip > "${INPUT}_R2.fastq.gz"
    INPUT_R2="${INPUT}_R2.fastq.gz"
fi

# Check that sample files exist
if test ! -f "${SAMPLE_R1}" || test ! -f "${SAMPLE_R2}" ; then
    SAMPLE_R1="${SAMPLE}.end1.fastq.gz"
    SAMPLE_R2="${SAMPLE}.end2.fastq.gz"
    if test -f "${SAMPLE_R1}" && test -f "${SAMPLE_R2}" ; then
        fn_warning "Sample files found here: ${SAMPLE_R1} & ${SAMPLE_R2}" 2>&1 | tee -a "${LOGFILE}" 
        fn_warning "Renaming '\${SAMPLE_R1}' & '\${SAMPLE_R2}' variables" 2>&1 | tee -a "${LOGFILE}" 
    else
        fn_error "Sample files are missing. Check sample directory: ${SAMPLE_DIR}/." 2>&1 | tee -a "${LOGFILE}"
        fn_error "Files *must* be named as follows: ${SAMPLE_BASE}_R1.fastq.gz & ${SAMPLE_BASE}_R2.fastq.gz" 2>&1 | tee -a "${LOGFILE}"
        fn_error "Aborting now." 2>&1 | tee -a "${LOGFILE}"
        rm --force "${LOGFILE}"
        exit 1
    fi
fi

# Check if a genome is provided
if test `is_set "${GENOME}"` == 1 ; then
    fn_error "Please provide a genome." 2>&1 | tee -a "${LOGFILE}"
    fn_error "Aborting now." 2>&1 | tee -a "${LOGFILE}"
    rm --force "${LOGFILE}"
    exit 1
fi

# Check that the genome fasta file exists
if test ! -f "${GENOME_FA}" ; then
    fn_error ""${GENOME_FA}" does not exist." 2>&1 | tee -a "${LOGFILE}"
    fn_error "Aborting now." 2>&1 | tee -a "${LOGFILE}"
    rm --force "${LOGFILE}"
    exit 1
fi

# If providing input, check that the input files exist
if test "${DO_INPUT}" == 0 ; then
    if test ! -f "${INPUT_R1}" || test ! -f "${INPUT_R2}" ; then
        INPUT_R1="${INPUT}.end1.fastq.gz"
        INPUT_R2="${INPUT}.end2.fastq.gz"
        if test -f "${INPUT_R1}" && test -f "${INPUT_R2}" ; then
            fn_warning "Sample files found here: ${INPUT_R1} & ${INPUT_R2}" 2>&1 | tee -a "${LOGFILE}" 
            fn_warning "Renaming '\${INPUT_R1}' & '\${INPUT_R2}' variables" 2>&1 | tee -a "${LOGFILE}" 
        else
            fn_error "Input files are missing. Check input directory: ${INPUT_DIR}/." 2>&1 | tee -a "${LOGFILE}"
            fn_error "Files *must* be named as follows: ${INPUT_BASE}_R1.fastq.gz & ${INPUT_BASE}_R2.fastq.gz" 2>&1 | tee -a "${LOGFILE}"
            fn_error "Aborting now." 2>&1 | tee -a "${LOGFILE}"
            rm --force "${LOGFILE}"
            exit 1
        fi
    fi
fi

# If providing calibration, check that the calibration genome exists
if test "${DO_CALIBRATION}" == 0 ; then
    if test ! -f "${SPIKEIN_FA}" ; then
        fn_error ""${SPIKEIN_FA}" does not exist." 2>&1 | tee -a "${LOGFILE}"
        fn_error "Aborting now." 2>&1 | tee -a "${LOGFILE}"
        rm --force "${LOGFILE}"
        exit 1
    fi
fi

# If provided reference genome is not indexed abort
if test ! -f "${GENOME_BASE}".1.bt2 || test ! -f "${GENOME_BASE}".2.bt2 || test ! -f "${GENOME_BASE}".3.bt2 || test ! -f "${GENOME_BASE}".4.bt2 || test ! -f "${GENOME_BASE}".rev.1.bt2 || test ! -f "${GENOME_BASE}".rev.2.bt2 ; then
    fn_error "Genome bowtie2 index files are missing. Please run the following command first:" 2>&1 | tee -a "${LOGFILE}"
    echo -e "bowtie2-build ${GENOME_FA} ${GENOME_BASE}" 2>&1 | tee -a "${LOGFILE}"
    fn_error "Aborting now." 2>&1 | tee -a "${LOGFILE}"
    rm --force "${LOGFILE}"
    exit 1
fi

# If provided calibration genome is not indexed abort
if test "${DO_CALIBRATION}" == 0 ; then
    if test ! -f "${SPIKEIN_BASE}".1.bt2 || test ! -f "${SPIKEIN_BASE}".2.bt2 || test ! -f "${SPIKEIN_BASE}".3.bt2 || test ! -f "${SPIKEIN_BASE}".4.bt2 || test ! -f "${SPIKEIN_BASE}".rev.1.bt2 || test ! -f "${SPIKEIN_BASE}".rev.2.bt2 ; then
        fn_error "Calibration genome bowtie2 index files are missing. Please run the following command first:" 2>&1 | tee -a "${LOGFILE}"
        echo -e "bowtie2-build ${SPIKEIN_FA} ${SPIKEIN_BASE}" 2>&1 | tee -a "${LOGFILE}"
        fn_error "Aborting now." 2>&1 | tee -a "${LOGFILE}"
        rm --force "${LOGFILE}"
        exit 1
    fi
fi

# Check that blacklist file exists, if provided
if test `is_set "${BLACKLISTBEDFILE}"` == 0 && test ! -f "${BLACKLISTBEDFILE}" ; then
    fn_error "The provided blaklist file does not exist." 2>&1 | tee -a "${LOGFILE}"
    fn_error "Aborting now." 2>&1 | tee -a "${LOGFILE}"
    rm --force "${LOGFILE}"
    exit 1
fi

# Check that utils are available
for util in bowtie2 samtools deeptools
do
    if test -z `command -v "${util}"` ; then
        fn_error "${util} does not seem to be installed or loaded. Most likely, it can be installed as follows:" 2>&1 | tee -a "${LOGFILE}"
        echo -e "conda install -c bioconda ${util}" 2>&1 | tee -a "${LOGFILE}"
        fn_error "Aborting now." 2>&1 | tee -a "${LOGFILE}"
        rm --force "${LOGFILE}"
        exit 1
    fi
done

if test "${DO_PEAKS}" == 0 ; then
    util=macs2
    if test -z `command -v "${util}"` ; then
        fn_error "${util} does not seem to be installed or loaded. Most likely, it can be installed as follows:" 2>&1 | tee -a "${LOGFILE}"
        echo -e "conda install -c bioconda ${util}" 2>&1 | tee -a "${LOGFILE}"
        fn_error "Aborting now." 2>&1 | tee -a "${LOGFILE}"
        rm --force "${LOGFILE}"
        exit 1
    fi
fi

if test "${MODE}" == HiC ; then
    util=hicstuff
    if test -z `command -v "${util}"` ; then
        fn_error "${util} does not seem to be installed or loaded. Install it as follows:" 2>&1 | tee -a "${LOGFILE}"
        echo -e "pip install -U ${util}" 2>&1 | tee -a "${LOGFILE}"
        fn_error "Aborting now." 2>&1 | tee -a "${LOGFILE}"
        rm --force "${LOGFILE}"
        exit 1
    fi
fi

## ------------------------------------------------------------------
## -------- PREPARING RESULT DIRECTORIES ----------------------------
## ------------------------------------------------------------------

mkdir -p "${OUTDIR}"/fastq/
mkdir -p "${OUTDIR}"/fastq/genome/
mkdir -p "${OUTDIR}"/fastq/genome/"${SAMPLE_BASE}"/
mkdir -p "${OUTDIR}"/fastq/genome/"${INPUT_BASE}"/
mkdir -p "${OUTDIR}"/fastq/spikein/
mkdir -p "${OUTDIR}"/fastq/spikein/"${SAMPLE_BASE}"/
mkdir -p "${OUTDIR}"/fastq/spikein/"${INPUT_BASE}"/
mkdir -p "${OUTDIR}"/bam/
mkdir -p "${OUTDIR}"/bam/genome/
mkdir -p "${OUTDIR}"/bam/genome/"${SAMPLE_BASE}"/
mkdir -p "${OUTDIR}"/bam/genome/"${INPUT_BASE}"/
mkdir -p "${OUTDIR}"/bam/spikein/
mkdir -p "${OUTDIR}"/bam/spikein/"${SAMPLE_BASE}"/
mkdir -p "${OUTDIR}"/bam/spikein/"${INPUT_BASE}"/
mkdir -p "${OUTDIR}"/tracks/
mkdir -p "${OUTDIR}"/tracks/"${SAMPLE_BASE}"
mkdir -p "${OUTDIR}"/tracks/"${INPUT_BASE}"
mkdir -p "${OUTDIR}"/peaks/
mkdir -p "${OUTDIR}"/peaks/"${SAMPLE_BASE}"
mkdir -p "${OUTDIR}"/matrices/
mkdir -p "${OUTDIR}"/matrices/"${SAMPLE_BASE}"/
mkdir -p "${OUTDIR}"/pairs/
mkdir -p "${OUTDIR}"/pairs/"${SAMPLE_BASE}"/
mkdir -p "${OUTDIR}"/stats/
mkdir -p "${OUTDIR}"/logs/

## ------------------------------------------------------------------
## -------- PRINT STARTUP INFO --------------------------------------
## ------------------------------------------------------------------

fn_log "Pipeline started on `date`" 2>&1 | tee "${LOGFILE}"
fn_log "Command     : ${INVOC}" 2>&1 | tee -a "${LOGFILE}"
fn_log "Hash        : ${HASH}" 2>&1 | tee -a "${LOGFILE}"
fn_log "Log file    : ${LOGFILE}" 2>&1 | tee -a "${LOGFILE}"
fn_log "Version     : ${VERSION}" 2>&1 | tee -a "${LOGFILE}"
echo -e "---" 2>&1 | tee -a "${LOGFILE}"
fn_log "MODE        : ${MODE}" 2>&1 | tee -a "${LOGFILE}"
fn_log "SAMPLE      : ${SAMPLE}" 2>&1 | tee -a "${LOGFILE}"
fn_log "GENOME      : ${GENOME}" 2>&1 | tee -a "${LOGFILE}"
if test "${DO_INPUT}" == 0 && test "${MODE}" == ChIP ; then
    fn_log "INPUT       : ${INPUT}" 2>&1 | tee -a "${LOGFILE}"
else 
    fn_warning "Input reads not provided. Processing without input." 2>&1 | tee -a "${LOGFILE}"
fi
if test "${DO_CALIBRATION}" == 0 && test "${MODE}" == ChIP ; then
    fn_log "SPIKEIN     : ${SPIKEIN}" 2>&1 | tee -a "${LOGFILE}"
else
    fn_warning "Spikein genome not provided. Processing without calibration." 2>&1 | tee -a "${LOGFILE}"
fi
fn_log "Keep dups.  : `if test ${KEEPDUPLICATES} == 0 ; then echo yes ; else echo no ; fi`" 2>&1 | tee -a "${LOGFILE}"
fn_log "Align. opt. : ${BOWTIEOPTIONS}" 2>&1 | tee -a "${LOGFILE}"
fn_log "Filt. opt.  : ${FILTEROPTIONS}" 2>&1 | tee -a "${LOGFILE}"
fn_log "CPU         : ${CPU}" 2>&1 | tee -a "${LOGFILE}"
fn_log "OUTDIR      : ${OUTDIR}" 2>&1 | tee -a "${LOGFILE}"
echo -e "---" 2>&1 | tee -a "${LOGFILE}"
fn_log "bowtie2     : `type -P bowtie2` (version: `bowtie2 --version | head -n1 | sed 's,.* ,,g'`)" 2>&1 | tee -a "${LOGFILE}"
fn_log "samtools    : `type -P samtools` (version: `samtools --version | head -n1 | sed 's,.* ,,'`)" 2>&1 | tee -a "${LOGFILE}"
fn_log "deeptools   : `type -P deeptools` (version: `deeptools --version | head -n1 | sed 's,.* ,,g'`)" 2>&1 | tee -a "${LOGFILE}"
fn_log "macs2       : `type -P macs2` (version: `macs2 --version | head -n1 | sed 's,.* ,,g'`)" 2>&1 | tee -a "${LOGFILE}"
fn_log "hicstuff    : `type -P hicstuff` (version: `hicstuff --version | head -n1 | sed 's,.* ,,g'`)" 2>&1 | tee -a "${LOGFILE}"
fn_log "cooler      : `type -P cooler` (version: `cooler --version | head -n1 | sed 's,.* ,,g'`)" 2>&1 | tee -a "${LOGFILE}"
echo -e "---" 2>&1 | tee -a "${LOGFILE}"

## ------------------------------------------------------------------
## ------------------- RUN HICSTUFF ---------------------------------
## ------------------------------------------------------------------

if test "${MODE}" == HiC ; then

    mkdir -p "${OUTDIR}"/tmp/
    cp "${GENOME_BASE}".fa "${OUTDIR}"/tmp/${SAMPLE_BASE}.genome.fasta
    cp "${GENOME_BASE}".1.bt2 "${OUTDIR}"/tmp/${SAMPLE_BASE}.genome.fasta.1.bt2
    cp "${GENOME_BASE}".2.bt2 "${OUTDIR}"/tmp/${SAMPLE_BASE}.genome.fasta.2.bt2
    cp "${GENOME_BASE}".3.bt2 "${OUTDIR}"/tmp/${SAMPLE_BASE}.genome.fasta.3.bt2
    cp "${GENOME_BASE}".4.bt2 "${OUTDIR}"/tmp/${SAMPLE_BASE}.genome.fasta.4.bt2
    cp "${GENOME_BASE}".rev.1.bt2 "${OUTDIR}"/tmp/${SAMPLE_BASE}.genome.fasta.rev.1.bt2
    cp "${GENOME_BASE}".rev.2.bt2 "${OUTDIR}"/tmp/${SAMPLE_BASE}.genome.fasta.rev.2.bt2

    fn_log "Processing sample reads with hicstuff" 2>&1 | tee -a "${LOGFILE}"
    cmd="hicstuff pipeline \
        --threads "${CPU}" \
        --enzyme "${RE}" \
        --outdir "${OUTDIR}" \
        --prefix "${SAMPLE_BASE}" \
        "${HICSTUFFOPTIONS}" --force \
        --matfmt cool \
        --genome "${OUTDIR}"/tmp/${SAMPLE_BASE}.genome.fasta \
        "${SAMPLE_R1}" "${SAMPLE_R2}""
    fn_exec "${cmd}" "${LOGFILE}" 2>> "${LOGFILE}"

    fn_log "Binning cool file to ${FIRSTREZ} bp" 2>&1 | tee -a "${LOGFILE}"
    cmd="hicstuff rebin \
        --binning "${FIRSTREZ}" \
        --frags "${OUTDIR}"/"${SAMPLE_BASE}".frags.tsv \
        --chroms "${OUTDIR}"/"${SAMPLE_BASE}".chr.tsv \
        --force \
        "${OUTDIR}"/"${SAMPLE_BASE}".cool \
        "${OUTDIR}"/"${SAMPLE_BASE}"_"${FIRSTREZ}""
    fn_exec "${cmd}" "${LOGFILE}" 2>> "${LOGFILE}"

    fn_log "Generating .mcool file" 2>&1 | tee -a "${LOGFILE}"
    cmd="cooler zoomify \
        --nproc "${CPU}" \
        --resolutions "${HICREZ}" \
        --balance \
        --balance-args \"--cis-only --min-nnz 3 --mad-max 7\" \
        --out "${SAMPLE_MCOOL}" \
        "${OUTDIR}"/"${SAMPLE_BASE}"_"${FIRSTREZ}".cool"
    fn_exec "${cmd}" "${LOGFILE}" 2>> "${LOGFILE}"

## ------------------------------------------------------------------
## ------------------- MAPPING --------------------------------------
## ------------------------------------------------------------------

else 

fn_log "Mapping sample reads to reference genome" 2>&1 | tee -a "${LOGFILE}"
cmd="bowtie2 ${BOWTIEOPTIONS} \
    --threads "${CPU}" \
    -x "${GENOME_BASE}" \
    -1 "${SAMPLE_R1}" \
    -2 "${SAMPLE_R2}" \
    --un-conc-gz "${SAMPLE_NON_ALIGNED_GENOME}".gz \
    > "${SAMPLE_ALIGNED_GENOME}"" 
fn_exec "${cmd}" "${LOGFILE}" 2>> "${LOGFILE}"

if test "${DO_CALIBRATION}" == 0 ; then
    fn_log "Mapping sample reads to spikein genome" 2>&1 | tee -a "${LOGFILE}"
    cmd="bowtie2 ${BOWTIEOPTIONS} \
        --threads "${CPU}" \
        -x "${SPIKEIN_BASE}" \
        -1 "${SAMPLE_R1}" \
        -2 "${SAMPLE_R2}" \
        --un-conc-gz "${SAMPLE_NON_ALIGNED_CALIBRATION}".gz \
        > "${SAMPLE_ALIGNED_CALIBRATION}""
        fn_exec "${cmd}" "${LOGFILE}" 2>> "${LOGFILE}"
    
    fn_log "Mapping sample reads non-mapped on spikein genome to reference genome" 2>&1 | tee -a "${LOGFILE}"
    cmd="bowtie2 ${BOWTIEOPTIONS} \
        --threads "${CPU}" \
        -x "${GENOME_BASE}" \
        -1 "${SAMPLE_NON_ALIGNED_CALIBRATION}".1.gz \
        -2 "${SAMPLE_NON_ALIGNED_CALIBRATION}".2.gz \
        > "${SAMPLE_NON_ALIGNED_CALIBRATION_ALIGNED_GENOME}""
    fn_exec "${cmd}" "${LOGFILE}" 2>> "${LOGFILE}"

    fn_log "Mapping sample reads non-mapped on reference genome to spikein genome" 2>&1 | tee -a "${LOGFILE}"
    cmd="bowtie2 ${BOWTIEOPTIONS} \
        --threads "${CPU}" \
        -x "${SPIKEIN_BASE}" \
        -1 "${SAMPLE_NON_ALIGNED_GENOME}".1.gz \
        -2 "${SAMPLE_NON_ALIGNED_GENOME}".2.gz \
        > "${SAMPLE_NON_ALIGNED_GENOME_ALIGNED_CALIBRATION}""
    fn_exec "${cmd}" "${LOGFILE}" 2>> "${LOGFILE}"
fi

if test "${DO_INPUT}" == 0 ; then
    fn_log "Mapping input reads to reference genome" 2>&1 | tee -a "${LOGFILE}"
    cmd="bowtie2 ${BOWTIEOPTIONS} \
        --threads "${CPU}" \
        -x "${GENOME_BASE}" \
        -1 "${INPUT_R1}" \
        -2 "${INPUT_R2}" \
        --un-conc-gz "${INPUT_NON_ALIGNED_GENOME}".gz \
        > "${INPUT_ALIGNED_GENOME}""
    fn_exec "${cmd}" "${LOGFILE}" 2>> "${LOGFILE}"

    if test "${DO_CALIBRATION}" == 0 ; then
        fn_log "Mapping input reads to spikein genome" 2>&1 | tee -a "${LOGFILE}"
        cmd="bowtie2 ${BOWTIEOPTIONS} \
            --threads "${CPU}" \
            -x "${SPIKEIN_BASE}" \
            -1 "${INPUT_R1}" \
            -2 "${INPUT_R2}" \
            --un-conc-gz "${INPUT_NON_ALIGNED_CALIBRATION}".gz \
            > "${INPUT_ALIGNED_CALIBRATION}""
        fn_exec "${cmd}" "${LOGFILE}" 2>> "${LOGFILE}"

        fn_log "Mapping input reads non-mapped on spikein genome to reference genome" 2>&1 | tee -a "${LOGFILE}"
        cmd="bowtie2 ${BOWTIEOPTIONS} \
            --threads "${CPU}" \
            -x "${GENOME_BASE}" \
            -1 "${INPUT_NON_ALIGNED_CALIBRATION}".1.gz \
            -2 "${INPUT_NON_ALIGNED_CALIBRATION}".2.gz \
            > "${INPUT_NON_ALIGNED_CALIBRATION_ALIGNED_GENOME}""
            fn_exec "${cmd}" "${LOGFILE}" 2>> "${LOGFILE}"

        fn_log "Mapping input reads non-mapped on reference genome to spikein genome" 2>&1 | tee -a "${LOGFILE}"
        cmd="bowtie2 ${BOWTIEOPTIONS} \
            --threads "${CPU}" \
            -x "${SPIKEIN_BASE}" \
            -1 "${INPUT_NON_ALIGNED_GENOME}".1.gz \
            -2 "${INPUT_NON_ALIGNED_GENOME}".2.gz \
            > "${INPUT_NON_ALIGNED_GENOME_ALIGNED_CALIBRATION}""
        fn_exec "${cmd}" "${LOGFILE}" 2>> "${LOGFILE}"
    fi

fi

## ------------------------------------------------------------------
## ------------------- FILTERING & INDEXING -------------------------
## ------------------------------------------------------------------

if test "${DO_CALIBRATION}" == 1 && test "${MODE}" != HiC ; then

    fn_log "Filtering sample bam file of reads mapped to reference genome" 2>&1 | tee -a "${LOGFILE}"
    cmd="samtools fixmate ${SAMTOOLS_OPTIONS} -m "${SAMPLE_ALIGNED_GENOME}" - \
        | samtools sort ${SAMTOOLS_OPTIONS} -T "${SAMPLE_ALIGNED_GENOME}"_sorting - \
        | samtools markdup ${SAMTOOLS_OPTIONS} ${REMOVE_DUPLICATES} -T "${SAMPLE_ALIGNED_GENOME}"_markdup - - \
        | samtools view ${SAMTOOLS_OPTIONS} ${FILTEROPTIONS} -1 -b - \
        | samtools sort ${SAMTOOLS_OPTIONS} -l 9 -T "${SAMPLE_ALIGNED_GENOME}"_sorting2 \
        -o "${SAMPLE_ALIGNED_GENOME_FILTERED}""
    fn_exec "${cmd}" "${LOGFILE}"
    cmd="samtools index -@ "${CPU}" "${SAMPLE_ALIGNED_GENOME_FILTERED}""
    fn_exec "${cmd}" "${LOGFILE}" 2>> "${LOGFILE}"

    if test "${DO_INPUT}" == 0 ; then

        fn_log "Filtering input bam file of reads mapped to reference genome" 2>&1 | tee -a "${LOGFILE}"
        cmd="samtools fixmate ${SAMTOOLS_OPTIONS} -m "${INPUT_ALIGNED_GENOME}" - \
            | samtools sort ${SAMTOOLS_OPTIONS} -T "${INPUT_ALIGNED_GENOME}"_sorting - \
            | samtools markdup ${SAMTOOLS_OPTIONS} ${REMOVE_DUPLICATES} -T "${INPUT_ALIGNED_GENOME}"_markdup - - \
            | samtools view ${SAMTOOLS_OPTIONS} ${FILTEROPTIONS} -1 -b - \
            | samtools sort ${SAMTOOLS_OPTIONS} -l 9 -T "${INPUT_ALIGNED_GENOME}"_sorting2 \
            -o "${INPUT_ALIGNED_GENOME_FILTERED}""
        fn_exec "${cmd}" "${LOGFILE}"
        cmd="samtools index -@ "${CPU}" "${INPUT_ALIGNED_GENOME_FILTERED}""
        fn_exec "${cmd}" "${LOGFILE}" 2>> "${LOGFILE}"

    fi

fi

if test "${DO_CALIBRATION}" == 0 ; then

    fn_log "Filtering sample bam file of reads unmapped on spikein genome and mapped to reference genome" 2>&1 | tee -a "${LOGFILE}"
    cmd="samtools fixmate ${SAMTOOLS_OPTIONS} -m "${SAMPLE_NON_ALIGNED_CALIBRATION_ALIGNED_GENOME}" - \
        | samtools sort ${SAMTOOLS_OPTIONS} -T "${SAMPLE_NON_ALIGNED_CALIBRATION_ALIGNED_GENOME}"_sorting - \
        | samtools markdup ${SAMTOOLS_OPTIONS} ${REMOVE_DUPLICATES} -T "${SAMPLE_NON_ALIGNED_CALIBRATION_ALIGNED_GENOME}"_markdup - - \
        | samtools view ${SAMTOOLS_OPTIONS} ${FILTEROPTIONS} -1 -b - \
        | samtools sort ${SAMTOOLS_OPTIONS} -l 9 -T "${SAMPLE_NON_ALIGNED_CALIBRATION_ALIGNED_GENOME}"_sorting2 \
        -o "${SAMPLE_NON_ALIGNED_CALIBRATION_ALIGNED_GENOME_FILTERED}""
    fn_exec "${cmd}" "${LOGFILE}"
    cmd="samtools index -@ "${CPU}" "${SAMPLE_NON_ALIGNED_CALIBRATION_ALIGNED_GENOME_FILTERED}""
    fn_exec "${cmd}" "${LOGFILE}" 2>> "${LOGFILE}"

    fn_log "Filtering sample bam file of reads unmapped on reference genome and mapped to spikein genome" 2>&1 | tee -a "${LOGFILE}"
    cmd="samtools fixmate ${SAMTOOLS_OPTIONS} -m "${SAMPLE_NON_ALIGNED_GENOME_ALIGNED_CALIBRATION}" - \
        | samtools sort ${SAMTOOLS_OPTIONS} -T "${SAMPLE_NON_ALIGNED_GENOME_ALIGNED_CALIBRATION}"_sorting - \
        | samtools markdup ${SAMTOOLS_OPTIONS} ${REMOVE_DUPLICATES} -T "${SAMPLE_NON_ALIGNED_GENOME_ALIGNED_CALIBRATION}"_markdup - - \
        | samtools view ${SAMTOOLS_OPTIONS} ${FILTEROPTIONS} -1 -b - \
        | samtools sort ${SAMTOOLS_OPTIONS} -l 9 -T "${SAMPLE_NON_ALIGNED_GENOME_ALIGNED_CALIBRATION}"_sorting2 \
        -o "${SAMPLE_NON_ALIGNED_GENOME_ALIGNED_CALIBRATION_FILTERED}""
    fn_exec "${cmd}" "${LOGFILE}"
    cmd="samtools index -@ "${CPU}" "${SAMPLE_NON_ALIGNED_GENOME_ALIGNED_CALIBRATION_FILTERED}""
    fn_exec "${cmd}" "${LOGFILE}" 2>> "${LOGFILE}"

    fn_log "Filtering input bam file of reads unmapped on spikein genome and mapped to reference genome" 2>&1 | tee -a "${LOGFILE}"
    cmd="samtools fixmate ${SAMTOOLS_OPTIONS} -m "${INPUT_NON_ALIGNED_CALIBRATION_ALIGNED_GENOME}" - \
        | samtools sort ${SAMTOOLS_OPTIONS} -T "${INPUT_NON_ALIGNED_CALIBRATION_ALIGNED_GENOME}"_sorting - \
        | samtools markdup ${SAMTOOLS_OPTIONS} ${REMOVE_DUPLICATES} -T "${INPUT_NON_ALIGNED_CALIBRATION_ALIGNED_GENOME}"_markdup - - \
        | samtools view ${SAMTOOLS_OPTIONS} ${FILTEROPTIONS} -1 -b - \
        | samtools sort ${SAMTOOLS_OPTIONS} -l 9 -T "${INPUT_NON_ALIGNED_CALIBRATION_ALIGNED_GENOME}"_sorting2 \
        -o "${INPUT_NON_ALIGNED_CALIBRATION_ALIGNED_GENOME_FILTERED}""
    fn_exec "${cmd}" "${LOGFILE}"
    cmd="samtools index -@ "${CPU}" "${INPUT_NON_ALIGNED_CALIBRATION_ALIGNED_GENOME_FILTERED}""
    fn_exec "${cmd}" "${LOGFILE}" 2>> "${LOGFILE}"

    fn_log "Filtering input bam file of reads unmapped on reference genome and mapped to spikein genome" 2>&1 | tee -a "${LOGFILE}"
    cmd="samtools fixmate ${SAMTOOLS_OPTIONS} -m "${INPUT_NON_ALIGNED_GENOME_ALIGNED_CALIBRATION}" - \
        | samtools sort ${SAMTOOLS_OPTIONS} -T "${INPUT_NON_ALIGNED_GENOME_ALIGNED_CALIBRATION}"_sorting - \
        | samtools markdup ${SAMTOOLS_OPTIONS} ${REMOVE_DUPLICATES} -T "${INPUT_NON_ALIGNED_GENOME_ALIGNED_CALIBRATION}"_markdup - - \
        | samtools view ${SAMTOOLS_OPTIONS} ${FILTEROPTIONS} -1 -b - \
        | samtools sort ${SAMTOOLS_OPTIONS} -l 9 -T "${INPUT_NON_ALIGNED_GENOME_ALIGNED_CALIBRATION}"_sorting2 \
        -o "${INPUT_NON_ALIGNED_GENOME_ALIGNED_CALIBRATION_FILTERED}""
    fn_exec "${cmd}" "${LOGFILE}"
    cmd="samtools index -@ "${CPU}" "${INPUT_NON_ALIGNED_GENOME_ALIGNED_CALIBRATION_FILTERED}""
    fn_exec "${cmd}" "${LOGFILE}" 2>> "${LOGFILE}"

fi

if test "${MODE}" == MNase ; then 
    fn_log "Further filtering sample bam file of reads mapped to reference genome for fragment size ("${MNASE_MINSIZE}"-"${MNASE_MAXSIZE}" bp)" 2>&1 | tee -a "${LOGFILE}"
    cmd="samtools view -@ "${CPU}" -h "${SAMPLE_ALIGNED_GENOME_FILTERED}" \
        | mawk '/^@/ || (sqrt((\$9^2)) > "${MNASE_MINSIZE}" && sqrt((\$9^2)) < "${MNASE_MAXSIZE}")' \
        | samtools view -b - > "${SAMPLE_ALIGNED_GENOME_FILTERED_READSIZE}""
    fn_exec "${cmd}" "${LOGFILE}"
    cmd="samtools index -@ "${CPU}" "${SAMPLE_ALIGNED_GENOME_FILTERED_READSIZE}""
    fn_exec "${cmd}" "${LOGFILE}" 2>> "${LOGFILE}"
fi

## ------------------------------------------------------------------
## ------------------- GETTING STATS AND SCALING FACTOR -------------
## ------------------------------------------------------------------

if test "${DO_CALIBRATION}" == 0 ; then

    fn_log "Calculating scaling factor for spikein-based normalization of ${SAMPLE_BASE}" 2>&1 | tee -a "${LOGFILE}"
    N_SAMPLE_GENOME=`samtools flagstat "${SAMPLE_NON_ALIGNED_CALIBRATION_ALIGNED_GENOME_FILTERED}" | grep 'paired in sequencing' | sed 's, .*,,'`
    N_INPUT_GENOME=`samtools flagstat "${INPUT_NON_ALIGNED_CALIBRATION_ALIGNED_GENOME_FILTERED}" | grep 'paired in sequencing' | sed 's, .*,,'`
    N_SAMPLE_SPIKEIN=`samtools flagstat "${SAMPLE_NON_ALIGNED_GENOME_ALIGNED_CALIBRATION_FILTERED}" | grep 'paired in sequencing' | sed 's, .*,,'`
    N_INPUT_SPIKEIN=`samtools flagstat "${INPUT_NON_ALIGNED_GENOME_ALIGNED_CALIBRATION_FILTERED}" | grep 'paired in sequencing' | sed 's, .*,,'`
    ORI=`echo "(${N_SAMPLE_GENOME}/${N_SAMPLE_SPIKEIN})/(${N_INPUT_GENOME}/${N_INPUT_SPIKEIN})" | bc -l`
    SCALING_FACTOR=`echo "${ORI} / ( ${N_SAMPLE_GENOME} / 1000000)" | bc -l`
    echo -e ""${SAMPLE_BASE}_${GENOME}_filtered"\t${N_SAMPLE_GENOME}" > "${STATFILE}"
    echo -e ""${INPUT_BASE}_${GENOME}_filtered"\t${N_INPUT_GENOME}" >> "${STATFILE}"
    echo -e ""${SAMPLE_BASE}_${SPIKEIN}_filtered"\t${N_SAMPLE_SPIKEIN}" >> "${STATFILE}"
    echo -e ""${INPUT_BASE}_${SPIKEIN}_filtered"\t${N_INPUT_SPIKEIN}" >> "${STATFILE}"
    echo -e "ORI\t${ORI}" >> "${STATFILE}"
    echo -e "Scaling\t${SCALING_FACTOR}" >> "${STATFILE}"

fi

## ------------------------------------------------------------------
## ------------------- TRACKS ---------------------------------------
## ------------------------------------------------------------------

if test "${DO_CALIBRATION}" == 1 && test "${MODE}" != HiC ; then

    fn_log "Generating CPM track for ${SAMPLE_BASE}" 2>&1 | tee -a "${LOGFILE}"
    cmd="bamCoverage \
        --bam "${SAMPLE_ALIGNED_GENOME_FILTERED}" \
        --outFileName "${SAMPLE_RAW_TRACK}" \
        --binSize 1 \
        --numberOfProcessors "${CPU}" \
        "${BLACKLIST_OPTIONS}" \
        --normalizeUsing CPM \
        --skipNonCoveredRegions \
        --extendReads \
        --ignoreDuplicates"
    fn_exec "${cmd}" "${LOGFILE}" 2>> "${LOGFILE}"

    if test "${DO_INPUT}" == 0 ; then

        fn_log "Generating CPM track for ${INPUT_BASE}" 2>&1 | tee -a "${LOGFILE}"
        cmd="bamCoverage \
            --bam "${INPUT_ALIGNED_GENOME_FILTERED}" \
            --outFileName "${INPUT_RAW_TRACK}" \
            --binSize 1 \
            --numberOfProcessors "${CPU}" \
            "${BLACKLIST_OPTIONS}" \
            --normalizeUsing CPM \
            --skipNonCoveredRegions \
            --extendReads \
            --ignoreDuplicates"
        fn_exec "${cmd}" "${LOGFILE}" 2>> "${LOGFILE}"

        fn_log "Generating track for ${SAMPLE_BASE} divided by ${INPUT_BASE}" 2>&1 | tee -a "${LOGFILE}"
        cmd="bamCompare \
            -b1 "${SAMPLE_ALIGNED_GENOME_FILTERED}" \
            -b2 "${INPUT_ALIGNED_GENOME_FILTERED}" \
            --outFileName "${SAMPLE_INPUT_TRACK}" \
            --scaleFactorsMethod readCount \
            --operation ratio \
            --skipZeroOverZero \
            --skipNAs \
            --numberOfProcessors "${CPU}" \
            "${BLACKLIST_OPTIONS}" \
            --binSize 5 \
            --skipNonCoveredRegions \
            --ignoreDuplicates"
        fn_exec "${cmd}" "${LOGFILE}" 2>> "${LOGFILE}"

        fi

fi

if test "${DO_CALIBRATION}" == 0 ; then

    fn_log "Generating CPM track for ${SAMPLE_BASE}" 2>&1 | tee -a "${LOGFILE}"
    cmd="bamCoverage \
        --bam "${SAMPLE_NON_ALIGNED_CALIBRATION_ALIGNED_GENOME_FILTERED}" \
        --outFileName "${SAMPLE_RAW_TRACK}" \
        --binSize 1 \
        --numberOfProcessors "${CPU}" \
        "${BLACKLIST_OPTIONS}" \
        --normalizeUsing CPM \
        --skipNonCoveredRegions \
        --extendReads \
        --ignoreDuplicates"
    fn_exec "${cmd}" "${LOGFILE}" 2>> "${LOGFILE}"

    fn_log "Generating CPM track for ${INPUT_BASE}" 2>&1 | tee -a "${LOGFILE}"
    cmd="bamCoverage \
        --bam "${INPUT_NON_ALIGNED_CALIBRATION_ALIGNED_GENOME_FILTERED}" \
        --outFileName "${INPUT_RAW_TRACK}" \
        --binSize 1 \
        --numberOfProcessors "${CPU}" \
        "${BLACKLIST_OPTIONS}" \
        --normalizeUsing CPM \
        --skipNonCoveredRegions \
        --extendReads \
        --ignoreDuplicates"
    fn_exec "${cmd}" "${LOGFILE}" 2>> "${LOGFILE}"

    fn_log "Generating track for ${SAMPLE_BASE} divided by ${INPUT_BASE}" 2>&1 | tee -a "${LOGFILE}"
    cmd="bamCompare \
        -b1 "${SAMPLE_NON_ALIGNED_CALIBRATION_ALIGNED_GENOME_FILTERED}" \
        -b2 "${INPUT_NON_ALIGNED_CALIBRATION_ALIGNED_GENOME_FILTERED}" \
        --outFileName "${SAMPLE_INPUT_TRACK}" \
        --scaleFactorsMethod readCount \
        --operation ratio \
        --skipZeroOverZero \
        --skipNAs \
        --numberOfProcessors "${CPU}" \
        "${BLACKLIST_OPTIONS}" \
        --binSize 5 \
        --skipNonCoveredRegions \
        --ignoreDuplicates"
    fn_exec "${cmd}" "${LOGFILE}" 2>> "${LOGFILE}"

    fn_log "Generating CPM track for ${SAMPLE_BASE} scaled by spikein-derived factor (ORI)" 2>&1 | tee -a "${LOGFILE}"
    cmd="bamCoverage \
        --bam "${SAMPLE_NON_ALIGNED_CALIBRATION_ALIGNED_GENOME_FILTERED}" \
        --outFileName "${SAMPLE_SPIKEINSCALED_TRACK}" \
        --binSize 1 \
        --numberOfProcessors "${CPU}" \
        "${BLACKLIST_OPTIONS}" \
        --normalizeUsing 'None' \
        --scaleFactor "${SCALING_FACTOR}" \
        --skipNonCoveredRegions \
        --extendReads \
        --ignoreDuplicates"
    fn_exec "${cmd}" "${LOGFILE}" 2>> "${LOGFILE}"

fi

if test "${MODE}" == MNase ; then 

    fn_log "Generating CPM track for ${SAMPLE_BASE} filtered for reads >${MNASE_MINSIZE}bp & <${MNASE_MAXSIZE}bp" 2>&1 | tee -a "${LOGFILE}"
    cmd="bamCoverage \
        --bam "${SAMPLE_ALIGNED_GENOME_FILTERED_READSIZE}" \
        --outFileName "${SAMPLE_READSIZE_TRACK}" \
        --binSize 1 \
        --numberOfProcessors "${CPU}" \
        "${BLACKLIST_OPTIONS}" \
        --normalizeUsing CPM \
        --skipNonCoveredRegions \
        --extendReads \
        --ignoreDuplicates"
    fn_exec "${cmd}" "${LOGFILE}" 2>> "${LOGFILE}"

    fn_log "Generating nucleosome positioning track for ${SAMPLE_BASE}" 2>&1 | tee -a "${LOGFILE}"
    cmd="bamCoverage \
        --bam "${SAMPLE_ALIGNED_GENOME_FILTERED_READSIZE}" \
        --outFileName "${SAMPLE_NUCPOS_TRACK}" \
        --binSize 1 \
        --numberOfProcessors "${CPU}" \
        "${BLACKLIST_OPTIONS}" \
        --normalizeUsing CPM \
        --skipNonCoveredRegions \
        --extendReads \
        --ignoreDuplicates \
        --smoothLength 10 \
        --MNase"
    fn_exec "${cmd}" "${LOGFILE}" 2>> "${LOGFILE}"

fi

if test "${MODE}" == RNA ; then 

    fn_log "Generating forward track for ${SAMPLE_BASE}" 2>&1 | tee -a "${LOGFILE}"
    cmd="bamCoverage \
        --bam "${SAMPLE_ALIGNED_GENOME_FILTERED}" \
        --outFileName "${SAMPLE_TRACK_FWD}" \
        --binSize 1 \
        --numberOfProcessors 12 \
        "${BLACKLIST_OPTIONS}" \
        --normalizeUsing CPM \
        --skipNonCoveredRegions \
        --extendReads \
        --ignoreDuplicates \
        --filterRNAstrand forward"
    fn_exec "${cmd}" "${LOGFILE}" 2>> "${LOGFILE}"

    fn_log "Generating reverse track for ${SAMPLE_BASE}" 2>&1 | tee -a "${LOGFILE}"
    cmd="bamCoverage \
        --bam "${SAMPLE_ALIGNED_GENOME_FILTERED}" \
        --outFileName "${SAMPLE_TRACK_REV}" \
        --binSize 1 \
        --numberOfProcessors 12 \
        "${BLACKLIST_OPTIONS}" \
        --normalizeUsing CPM \
        --skipNonCoveredRegions \
        --extendReads \
        --ignoreDuplicates \
        --filterRNAstrand reverse"
    fn_exec "${cmd}" "${LOGFILE}" 2>> "${LOGFILE}"

fi

## ------------------------------------------------------------------
## ------------------- CALLING PEAKS --------------------------------
## ------------------------------------------------------------------

if test "${DO_PEAKS}" == 0 ; then

    fn_log "Calling peaks for ${SAMPLE_BASE}" 2>&1 | tee -a "${LOGFILE}"

    if test "${DO_INPUT}" == 1 ; then

        cmd="macs2 callpeak \
            -t "${SAMPLE_ALIGNED_GENOME_FILTERED}" \
            --format BAMPE \
            --gsize 13000000 \
            --outdir "${OUTDIR}"/peaks/"${SAMPLE_BASE}" \
            --name "${SAMPLE_BASE}_genome-${GENOME}_${HASH}""
        fn_exec "${cmd}" "${LOGFILE}" 2>> "${LOGFILE}"

    fi 

    if test "${DO_INPUT}" == 0 ; then

        cmd="macs2 callpeak \
            -t "${SAMPLE_NON_ALIGNED_CALIBRATION_ALIGNED_GENOME_FILTERED}" \
            -c "${INPUT_NON_ALIGNED_CALIBRATION_ALIGNED_GENOME_FILTERED}" \
            --format BAMPE \
            --gsize 13000000 \
            --outdir "${OUTDIR}"/peaks/"${SAMPLE_BASE}" \
            --name "${SAMPLE_BASE}_vs-${INPUT_BASE}_genome-${GENOME}_${HASH}""
        fn_exec "${cmd}" "${LOGFILE}" 2>> "${LOGFILE}"

    fi 

fi

fi # ------------------------------------- Exit HiC if...else statement

## ------------------------------------------------------------------
## ------------------- CHECK NB OF READS ----------------------------
## ------------------------------------------------------------------

echo -e "---" >> "${LOGFILE}"
fn_log "NUMBER OF SEQUENCED FRAGMENTS: ${SAMPLE_R1}" >> "${LOGFILE}"
echo `fastqfastcnt "${SAMPLE_R1}" frags.` >> "${LOGFILE}"
echo -e "---" >> "${LOGFILE}"

if test "${DO_CALIBRATION}" == 0 ; then
    files="${SAMPLE_ALIGNED_GENOME} ${SAMPLE_NON_ALIGNED_CALIBRATION_ALIGNED_GENOME_FILTERED}" 
elif test "${MODE}" == HiC ; then
    files="${OUTDIR}/tmp/${SAMPLE_BASE}.for.bam ${OUTDIR}/tmp/${SAMPLE_BASE}.rev.bam"
else 
    files="${SAMPLE_ALIGNED_GENOME} ${SAMPLE_ALIGNED_GENOME_FILTERED}"
fi

for file in ${files}
do
    fn_log "MAPPING STATS: ${file}" >> "${LOGFILE}"
    samtools flagstat "${file}" 2>&1 | tee -a "${LOGFILE}"
    echo -e "---" >> "${LOGFILE}"
done

if test "${DO_CALIBRATION}" == 0 ; then
    file="${STATFILE}"
    fn_log "CALIBRATION STATS: ${file}" >> "${LOGFILE}"
    cat "${file}" >> "${LOGFILE}"
    echo -e "---" >> "${LOGFILE}"
fi

## ------------------------------------------------------------------
## ------------------- CLEAN UP DIR ---------------------------------
## ------------------------------------------------------------------

if test "${MODE}" == HiC ; then
    rm --force "${OUTDIR}"/tmp/*bt2 "${OUTDIR}"/tmp/"${SAMPLE_BASE}".genome.fasta
    rm --force "${OUTDIR}"/"${SAMPLE_BASE}".frags.tsv "${OUTDIR}"/"${SAMPLE_BASE}".chr.tsv 
    rm --force "${OUTDIR}"/"${SAMPLE_BASE}".hicstuff*
    mv "${OUTDIR}"/"${SAMPLE_BASE}"_"${FIRSTREZ}".cool "${SAMPLE_COOL}"
    mv "${OUTDIR}"/tmp/"${SAMPLE_BASE}".for.bam "${SAMPLE_ALIGNED_GENOME_FWD}"
    mv "${OUTDIR}"/tmp/"${SAMPLE_BASE}".rev.bam "${SAMPLE_ALIGNED_GENOME_REV}"
    mv "${OUTDIR}"/tmp/"${SAMPLE_BASE}".valid.pairs "${OUTDIR}"/pairs/"${SAMPLE_BASE}"/"${SAMPLE_BASE}"^"${HASH}".valid.pairs
    mv "${OUTDIR}"/tmp/"${SAMPLE_BASE}".valid_idx.pairs "${OUTDIR}"/pairs/"${SAMPLE_BASE}"/"${SAMPLE_BASE}"^"${HASH}".valid_idx.pairs
    mv "${OUTDIR}"/tmp/"${SAMPLE_BASE}".valid_idx_filtered.pairs "${OUTDIR}"/pairs/"${SAMPLE_BASE}"/"${SAMPLE_BASE}"^"${HASH}".valid_idx_filtered.pairs
    mv "${OUTDIR}"/tmp/"${SAMPLE_BASE}".valid_idx_pcrfree.pairs "${OUTDIR}"/pairs/"${SAMPLE_BASE}"/"${SAMPLE_BASE}"^"${HASH}".valid_idx_pcrfree.pairs
    mv "${OUTDIR}"/plots/event_distance.pdf "${OUTDIR}"/pairs/"${SAMPLE_BASE}"/"${SAMPLE_BASE}"^"${HASH}".event_distance.pdf
    mv "${OUTDIR}"/plots/frags_hist.pdf "${OUTDIR}"/pairs/"${SAMPLE_BASE}"/"${SAMPLE_BASE}"^"${HASH}".frags_hist.pdf
fi

if test "${KEEPFILES}" == 1 ; then
    fn_log "Cleaning up temporary files for ${SAMPLE_BASE}" 2>&1 | tee -a "${LOGFILE}"
    rm --force "${SAMPLE_ALIGNED_GENOME}"
    rm --force "${SAMPLE_NON_ALIGNED_GENOME}"
    rm --force "${SAMPLE_ALIGNED_CALIBRATION}"
    rm --force "${SAMPLE_NON_ALIGNED_CALIBRATION}"
    rm --force "${SAMPLE_NON_ALIGNED_GENOME_ALIGNED_CALIBRATION}"
    rm --force "${SAMPLE_NON_ALIGNED_CALIBRATION_ALIGNED_GENOME}"
    rm --force "${INPUT_ALIGNED_GENOME}"
    rm --force "${INPUT_NON_ALIGNED_GENOME}"
    rm --force "${INPUT_ALIGNED_CALIBRATION}"
    rm --force "${INPUT_NON_ALIGNED_CALIBRATION}"
    rm --force "${INPUT_NON_ALIGNED_GENOME_ALIGNED_CALIBRATION}"
    rm --force "${INPUT_NON_ALIGNED_CALIBRATION_ALIGNED_GENOME}"
    rm --force "${OUTDIR}"/fastq/genome/"${SAMPLE_BASE}"/"${SAMPLE_BASE}"^unmapped_"${GENOME}"^"${HASH}"*
    rm --force "${OUTDIR}"/fastq/spikein/"${SAMPLE_BASE}"/"${SAMPLE_BASE}"^unmapped_"${SPIKEIN}"^"${HASH}"*
    rm --force "${OUTDIR}"/fastq/genome/"${INPUT_BASE}"/"${INPUT_BASE}"^unmapped_"${GENOME}"^"${HASH}"*
    rm --force "${OUTDIR}"/fastq/spikein/"${INPUT_BASE}"/"${INPUT_BASE}"^unmapped_"${SPIKEIN}"^"${HASH}"*
fi

## ------------------------------------------------------------------
## ------------------- WRAP THINGS UP -------------------------------
## ------------------------------------------------------------------

## -- Back up script file
fn_log "Backing up pipeline script" 2>&1 | tee -a "${LOGFILE}"
cp "${BASH_SOURCE}" "${SCRIPTFILE}"
chmod -x "${SCRIPTFILE}"

## -- List generated files
echo -e "---" 2>&1 | tee -a "${LOGFILE}"
fn_log "GENERATED FILES" >> "${LOGFILE}"
tree "${OUTDIR}" -P "*${HASH}*" --prune 2>&1 | tee -a "${LOGFILE}"

## -- List commands
echo -e "---" 2>&1 | tee -a "${LOGFILE}"
fn_log "COMMANDS EXECUTED" 2>&1 | tee -a "${LOGFILE}"
grep "\[EXEC\]" "${LOGFILE}" | sed 's,.*EXEC],,' | sed 's,| ,\\\n\t| ,g' 2>&1 | tee -a "${LOGFILE}"

## -- Done !
echo -e "---" 2>&1 | tee -a "${LOGFILE}"
fn_log "Pipeline achieved: `date`" 2>&1 | tee -a "${LOGFILE}"
echo -e "Do check the log file \`${LOGFILE}\` to make sure everything went ok!"
echo -e "You can re-run the commands by running \`./"${CMDFILE}"\`."

## -- Remove color decorators from log file
sed -i 's/\x1b\[[0-9;]*m//g' "${LOGFILE}"
grep "\[EXEC\]" "${LOGFILE}" | sed 's,.*EXEC],,' | sed 's,| ,\\\n\t| ,g' > "${CMDFILE}"
# sed -i 's/.*21m//g' "${CMDFILE}"

## -- Remove progress file and delete empty dirs
rm --force "${TMPFILE}"
if test `ls "${OUTDIR}"/*INPROGRESS 1> /dev/null 2>&1 ; echo $?` == 2 ; then
    find "${OUTDIR}" -type d -empty -delete
fi
