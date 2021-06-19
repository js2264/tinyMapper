#!/bin/bash

INVOC=$(printf %q "$BASH_SOURCE")$((($#)) && printf ' %q' "$@")
HASH=`LC_CTYPE=C tr -dc 'A-Z0-9' < /dev/urandom | head -c 6`

## ------------------------------------------------------------------
## -------- HELPER FUNCTIONS ----------------------------------------
## ------------------------------------------------------------------

function usage() {
    echo -e ""
    echo -e "# J. Serizay, C. Matthey-Doret, H. Bordelet"
    echo -e "# GPL-3.0"
    echo -e ""
    echo -e "Usage: tidyMapper.sh -m <MODE> -s <SAMPLE> -g <GENOME> -o <OUTPUT> [ -i <INPUT> | -c <CALIBRATION> | -t <THREADS> | -M <MEMORY> | -k <1, 0> ]"
    echo -e ""
    echo -e "-------------------------------"
    echo -e ""
    echo -e "   BASIC ARGUMENTS"
    echo -e ""
    echo -e "      -m | --mode                 Mapping mode (ChIP, MNase, ATAC, RNA) (Default: ChIP)"
    echo -e "      -s | --sample               Path prefix to sample <SAMPLE>_R*.fastq.gz (e.g. for ~/reads/JS001_R*.fastq.gz files: --sample ~/reads/JS001)"
    echo -e "      -g | --genome               Path prefix to reference genome (e.g. for ~/genome/W303/W303.fa fasta file: --genome ~/genome/W303/W303)"
    echo -e "      -o | --output               Path to store results (Default: ./results/)"
    echo -e "      -i | --input                (Optional) Path prefix to input <INPUT>_R*.fastq.gz"
    echo -e "      -c | --calibration          (Optional) Path prefix to genome used for calibration"
    echo -e "      -t | --threads              (Optional) Number of threads (Default: 8)"
    echo -e "      -M | --memory               (Optional) Memory in bits (Default: 12294967296, which is 12Gb)"
    echo -e "      -k | --keepIntermediate     (Optional) Keep intermediate mapping files (Default: 1 (i.e. 'false'))"
    echo -e "      -h | --help                 Print this message"
    echo -e ""
    echo -e "   -----------"
    echo -e ""
    echo -e "   ADVANCED ARGUMENTS"
    echo -e ""
    echo -e "      -f | --filter               Filtering options for \`samtools view\` (between single quotes)"
    echo -e "                                  Default: '-f 2 -q 10' (only keep paired reads and filter out reads with mapping quality score < 10)"
    echo -e ""
    echo -e "      -d | --duplicates           Keep duplicate reads"
    echo -e "                                  Default: 1 (i.e. 'false')"
    echo -e ""
    echo -e "-------------------------------"
    echo -e ""
    echo -e "Examples:"
    echo -e ""
    echo -e "   ChIP-seq mode:"
    echo -e ""
    echo -e "      Without input:               ./tidyMapper.sh -m ChIP -s ~/HB44 -g ~/genomes/R64-1-1/R64-1-1 -o ~/results"
    echo -e "      With input:                  ./tidyMapper.sh -m ChIP -s ~/HB44 -i ~/HB42 -g ~/genomes/R64-1-1/R64-1-1 -o ~/results"
    echo -e "      With input and calibration:  ./tidyMapper.sh -m ChIP -s ~/HB44 -i ~/HB42 -g ~/genomes/R64-1-1/R64-1-1 -c ~/genomes/Cglabrata/Cglabrata -o ~/results"
    echo -e ""
    echo -e "   RNA-seq mode:"
    echo -e ""
    echo -e "      ./tidyMapper.sh -m RNA -s ~/AB4 -g ~/genomes/W303/W303 -o ~/results"
    echo -e ""
    echo -e "   MNase-seq mode:"
    echo -e ""
    echo -e "      ./tidyMapper.sh -m MNase -s ~/CH266 -g ~/genomes/W303/W303 -o ~/results"
    echo -e ""
    echo -e "-------------------------------"
    echo -e ""
    echo -e "Required utilities:"
    echo -e ""
    echo -e "   bowtie2"
    echo -e "   samtools"
    echo -e "   deeptools"
    echo -e "   macs2"
    echo -e ""
}

function is_set() { 
    test -n "${1}" ; echo $?
}

function fn_log {
    date=`date "+%y-%m-%d %H:%M:%S"`
    BOLD="\e[1m"
    BOLDEND="\e[21m"
    GREEN="\e[32m"
    RED="\e[31m"
    BLUE="\e[96m"
    YELLOW="\e[33m"
    DEFAULT="\e[39m"
    echo -e "${BOLD}${BLUE}${date} | ${GREEN}[INFO]${DEFAULT} $@${BOLDEND}"
}

function fn_error {
    date=`date "+%y-%m-%d %H:%M:%S"`
    BOLD="\e[1m"
    BOLDEND="\e[21m"
    GREEN="\e[32m"
    RED="\e[31m"
    BLUE="\e[96m"
    YELLOW="\e[33m"
    DEFAULT="\e[39m"
    echo -e "${BOLD}${BLUE}${date} | ${RED}[ERR.]${DEFAULT} $@${BOLDEND}"
}

function fn_warning {
    date=`date "+%y-%m-%d %H:%M:%S"`
    BOLD="\e[1m"
    BOLDEND="\e[21m"
    GREEN="\e[32m"
    RED="\e[31m"
    ORANGE="\e[38;5;208m"
    BLUE="\e[96m"
    YELLOW="\e[33m"
    DEFAULT="\e[39m"
    echo -e "${BOLD}${BLUE}${date} | ${ORANGE}[WAR.]${DEFAULT} $@${BOLDEND}"
}

function fn_exec {
    date=`date "+%y-%m-%d %H:%M:%S"`
    BOLD="\e[1m"
    BOLDEND="\e[21m"
    GREEN="\e[32m"
    RED="\e[31m"
    BLUE="\e[96m"
    YELLOW="\e[33m"
    DEFAULT="\e[39m"
    cmd=`echo $1 | tr -s '' | sed 's,2>>.*,,' `
    echo -e "${BOLD}${BLUE}${date} | ${YELLOW}[EXEC]${DEFAULT} ${cmd}${BOLDEND}" >> $2
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

MODE=ChIP
SAMPLE=''
INPUT=''
GENOME=''
SPIKEIN=''
OUTDIR='results'
FILTEROPTIONS='-f 2 -q 10'
KEEPDUPLICATES=1
CPU=16
MEM=12294967296 # 12Gb
KEEPFILES=1

# MODE=ChIP
# SAMPLE=HB44
# INPUT=HB42
# GENOME=~/genomes/R64-1-1/R64-1-1
# SPIKEIN=~/genomes/Cglabrata_CBS138/Cglabrata_CBS138
# OUTDIR=results
# FILTEROPTIONS='-F 4 -q 5'
# KEEPDUPLICATES=0
# CPU=16
# MEM=12294967296 # 12Gb
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
        GENOME="${2}"
        shift 
        shift 
        ;;
        -c|--calibration)
        SPIKEIN="${2}"
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
        -f|--filter)
        FILTEROPTIONS=${2}
        shift 
        shift 
        ;;
        -d|--duplicates)
        KEEPDUPLICATES="${2}"
        shift 
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
        -M|--memory)
        MEM="${2}"
        shift 
        shift 
        ;;
        -k|--keepIntermediate)
        KEEPFILES="${2}"
        shift 
        shift 
        ;;
        -h|--help)
        usage && exit 0
        ;;
    esac
done

GENOME_DIR=`dirname "${GENOME}"`
GENOME=`basename "${GENOME}"`
GENOME_BASE="${GENOME_DIR}"/"${GENOME}"
GENOME_FA="${GENOME_BASE}.fa"
SPIKEIN_DIR=`dirname "${SPIKEIN}"`
SPIKEIN=`basename "${SPIKEIN}"`
SPIKEIN_BASE="${SPIKEIN_DIR}"/"${SPIKEIN}"
SPIKEIN_FA="${SPIKEIN_BASE}.fa"

SAMPLE_DIR=`dirname "${SAMPLE}"`
SAMPLE_BASE=`basename "${SAMPLE}"`
SAMPLE_R1="${SAMPLE}_R1.fastq.gz"
SAMPLE_R2="${SAMPLE}_R2.fastq.gz"
INPUT_DIR=`dirname "${INPUT}"`
INPUT_BASE=`basename "${INPUT}"`
INPUT_R1="${INPUT}_R1.fastq.gz"
INPUT_R2="${INPUT}_R2.fastq.gz"

LOGFILE="${OUTDIR}/`date "+%y%m%d"`-${HASH}-log.txt"
DO_INPUT=`is_set "${INPUT}"`
DO_CALIBRATION=`is_set "${SPIKEIN}"`
DO_PEAKS=`if test "${MODE}" == 'ChIP' || test "${MODE}" == 'ATAC'; then echo 0; else echo 1; fi`
REMOVE_DUPLICATES=`if test "${KEEPDUPLICATES}" == 1 ; then echo " -r " ; else echo " " ; fi`

## ------------------------------------------------------------------
## -------- CHECKING THAT REQUIRED FILES EXIST ----------------------
## ------------------------------------------------------------------

# Check the OUTDIR already exists
# if test -d "${OUTDIR}" ; then
#     fn_error "Output directory already exists (${OUTDIR}). Please erase or choose another directory to store output files."
#     fn_error "Aborting now."
#     exit 1
# fi

mkdir -p "${OUTDIR}"
touch "${LOGFILE}"

# Abort if trying to calibrate without input
if test "${DO_INPUT}" == 1 && test "${DO_CALIBRATION}" == 0 ; then
    fn_error "Calibration can only be done if an input is provided." 2>&1 | tee -a "${LOGFILE}"
    fn_error "Aborting now." 2>&1 | tee -a "${LOGFILE}"
    usage
    rm --force "${LOGFILE}"
    exit 1
fi

# Check if a sample is provided
if test `is_set "${SAMPLE_BASE}"` == 1 ; then
    fn_error "Please provide a sample." 2>&1 | tee -a "${LOGFILE}"
    fn_error "Aborting now." 2>&1 | tee -a "${LOGFILE}"
    usage
    rm --force "${LOGFILE}"
    exit 1
fi

# Check that sample files exist
if test ! -f "${SAMPLE_R1}" || test ! -f "${SAMPLE_R2}" ; then
    fn_error "Sample files are missing. Check sample directory: ${SAMPLE_DIR}." 2>&1 | tee -a "${LOGFILE}"
    fn_error "Files *must* be named as follows: ${SAMPLE_BASE}_R1.fastq.gz & ${SAMPLE_BASE}_R1.fastq.gz" 2>&1 | tee -a "${LOGFILE}"
    fn_error "Aborting now." 2>&1 | tee -a "${LOGFILE}"
    usage
    rm --force "${LOGFILE}"
    exit 1
fi

# Check if a genome is provided
if test `is_set "${GENOME}"` == 1 ; then
    fn_error "Please provide a genome." 2>&1 | tee -a "${LOGFILE}"
    fn_error "Aborting now." 2>&1 | tee -a "${LOGFILE}"
    usage
    rm --force "${LOGFILE}"
    exit 1
fi

# Check that the genome fasta file exists
if test ! -f "${GENOME_FA}" ; then
    fn_error ""${GENOME_FA}" does not exist." 2>&1 | tee -a "${LOGFILE}"
    fn_error "Aborting now." 2>&1 | tee -a "${LOGFILE}"
    usage
    rm --force "${LOGFILE}"
    exit 1
fi

# If providing input, check that the input files exist
if test "${DO_INPUT}" == 0 ; then
    if test ! -f "${INPUT_R1}" || test ! -f "${INPUT_R2}" ; then
        fn_error "Input files are missing. Check them in ${INPUT_DIR}." 2>&1 | tee -a "${LOGFILE}"
        fn_error "Files *must* be named as follows: ${INPUT_BASE}_R1.fastq.gz & ${INPUT_BASE}_R1.fastq.gz" 2>&1 | tee -a "${LOGFILE}"
        fn_error "Aborting now." 2>&1 | tee -a "${LOGFILE}"
        usage
        rm --force "${LOGFILE}"
        exit 1
    fi
fi

# If providing calibration, check that the calibration genome exists
if test "${DO_CALIBRATION}" == 0 ; then
    if test ! -f "${SPIKEIN_FA}" ; then
        fn_error ""${SPIKEIN_FA}" does not exist." 2>&1 | tee -a "${LOGFILE}"
        fn_error "Aborting now." 2>&1 | tee -a "${LOGFILE}"
        usage
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
mkdir -p "${OUTDIR}"/stats/

## ------------------------------------------------------------------
## -------- PREPARING FILE NAME VARIABLES ---------------------------
## ------------------------------------------------------------------

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

if test "${DO_CALIBRATION}" == 1 ; then
    SAMPLE_RAW_TRACK="${OUTDIR}"/tracks/"${SAMPLE_BASE}"/"${SAMPLE_BASE}"^mapped_"${GENOME}"^"${HASH}".CPM.bw
    INPUT_RAW_TRACK="${OUTDIR}"/tracks/"${INPUT_BASE}"/"${INPUT_BASE}"^mapped_"${GENOME}"^"${HASH}".CPM.bw
    SAMPLE_INPUT_TRACK="${OUTDIR}"/tracks/"${SAMPLE_BASE}"/"${SAMPLE_BASE}"^mapped_"${GENOME}"^"${HASH}".vs-"${INPUT_BASE}".bw
else
    SAMPLE_RAW_TRACK="${OUTDIR}"/tracks/"${SAMPLE_BASE}"/"${SAMPLE_BASE}"^unmapped_"${SPIKEIN}"^mapped_"${GENOME}"^"${HASH}".CPM.bw
    INPUT_RAW_TRACK="${OUTDIR}"/tracks/"${INPUT_BASE}"/"${INPUT_BASE}"^unmapped_"${SPIKEIN}"^mapped_"${GENOME}"^"${HASH}".CPM.bw
    SAMPLE_INPUT_TRACK="${OUTDIR}"/tracks/"${SAMPLE_BASE}"/"${SAMPLE_BASE}"^unmapped_"${SPIKEIN}"^mapped_"${GENOME}"^"${HASH}".vs-"${INPUT_BASE}".bw
    SAMPLE_SPIKEINSCALED_TRACK="${OUTDIR}"/tracks/"${SAMPLE_BASE}"/"${SAMPLE_BASE}"^unmapped_"${SPIKEIN}"^mapped_"${GENOME}"^"${HASH}".CPM.calibrated.bw
fi

if test "${MODE}" == MNase ; then 
    SAMPLE_ALIGNED_GENOME_FILTERED_READSIZE="${OUTDIR}"/bam/genome/"${SAMPLE_BASE}"/"${SAMPLE_BASE}"^mapped_"${GENOME}"^filtered^70-250^"${HASH}".bam
    SAMPLE_READSIZE_TRACK="${OUTDIR}"/tracks/"${SAMPLE_BASE}"/"${SAMPLE_BASE}"^mapped_"${GENOME}"^filtered^70-250^"${HASH}".70-250.CPM.bw
    SAMPLE_NUCPOS_TRACK="${OUTDIR}"/tracks/"${SAMPLE_BASE}"/"${SAMPLE_BASE}"^mapped_"${GENOME}"^filtered^70-250^"${HASH}".nucpos.CPM.bw
fi

if test "${MODE}" == RNA ; then 
    SAMPLE_RAW_TRACK="${OUTDIR}"/tracks/"${SAMPLE_BASE}"/"${SAMPLE_BASE}"^mapped_"${GENOME}"^"${HASH}".unstranded.CPM.bw
    SAMPLE_TRACK_FWD="${OUTDIR}"/tracks/"${SAMPLE_BASE}"/"${SAMPLE_BASE}"^mapped_"${GENOME}"^"${HASH}".fwd.CPM.bw
    SAMPLE_TRACK_REV="${OUTDIR}"/tracks/"${SAMPLE_BASE}"/"${SAMPLE_BASE}"^mapped_"${GENOME}"^"${HASH}".rev.CPM.bw
fi

SAMTOOLS_OPTIONS=" -@ ${CPU} --output-fmt bam "
STATFILE="${OUTDIR}/stats/sample-${SAMPLE_BASE}_input-${INPUT_BASE}_genome-${GENOME}_calibration-${SPIKEIN}_${HASH}".counts.tsv

## ------------------------------------------------------------------
## -------- PRINT STARTUP INFO --------------------------------------
## ------------------------------------------------------------------

fn_log "Pipeline started on `date`" 2>&1 | tee "${LOGFILE}"
fn_log "Command   : ${INVOC}" 2>&1 | tee -a "${LOGFILE}"
fn_log "Hash      : ${HASH}" 2>&1 | tee -a "${LOGFILE}"
fn_log "Log file  : ${LOGFILE}" 2>&1 | tee -a "${LOGFILE}"
echo -e "---" 2>&1 | tee -a "${LOGFILE}"
fn_log "MODE      : ${MODE}" 2>&1 | tee -a "${LOGFILE}"
fn_log "SAMPLE    : ${SAMPLE}" 2>&1 | tee -a "${LOGFILE}"
fn_log "GENOME    : ${GENOME}" 2>&1 | tee -a "${LOGFILE}"
if test "${DO_INPUT}" == 0 ; then
    fn_log "INPUT     : ${INPUT}" 2>&1 | tee -a "${LOGFILE}"
else 
    fn_warning "Input reads not provided. Processing without input." 2>&1 | tee -a "${LOGFILE}"
fi
if test "${DO_CALIBRATION}" == 0 ; then
    fn_log "SPIKEIN   : ${SPIKEIN}" 2>&1 | tee -a "${LOGFILE}"
else
    fn_warning "Spikein genome not provided. Processing without calibration." 2>&1 | tee -a "${LOGFILE}"
fi
fn_log "Keep dups.: ${KEEPDUPLICATES}" 2>&1 | tee -a "${LOGFILE}"
fn_log "Filt. opt.: ${FILTEROPTIONS}" 2>&1 | tee -a "${LOGFILE}"
fn_log "CPU       : ${CPU}" 2>&1 | tee -a "${LOGFILE}"
fn_log "MEM       : ${MEM}" 2>&1 | tee -a "${LOGFILE}"
fn_log "OUTDIR    : ${OUTDIR}" 2>&1 | tee -a "${LOGFILE}"
echo -e "---" 2>&1 | tee -a "${LOGFILE}"
fn_log "bowtie2   : `type -P bowtie2` (version: `bowtie2 --version | head -n1 | sed 's,.* ,,g'`)" 2>&1 | tee -a "${LOGFILE}"
fn_log "samtools  : `type -P samtools` (version: `samtools --version | head -n1 | sed 's,.* ,,'`)" 2>&1 | tee -a "${LOGFILE}"
fn_log "deeptools : `type -P deeptools` (version: `deeptools --version | head -n1 | sed 's,.* ,,g'`)" 2>&1 | tee -a "${LOGFILE}"
fn_log "macs2     : `type -P macs2` (version: `macs2 --version | head -n1 | sed 's,.* ,,g'`)" 2>&1 | tee -a "${LOGFILE}"
echo -e "---" 2>&1 | tee -a "${LOGFILE}"

## ------------------------------------------------------------------
## ------------------- MAPPING --------------------------------------
## ------------------------------------------------------------------

fn_log "Mapping sample reads to reference genome" 2>&1 | tee -a "${LOGFILE}"
cmd="bowtie2 \
    --threads "${CPU}" \
    -x "${GENOME_BASE}" \
    -1 "${SAMPLE_R1}" \
    -2 "${SAMPLE_R2}" \
    --un-conc-gz "${SAMPLE_NON_ALIGNED_GENOME}".gz \
    > "${SAMPLE_ALIGNED_GENOME}"" 
    fn_exec "${cmd}" "${LOGFILE}" 2>> "${LOGFILE}"

if test "${DO_CALIBRATION}" == 0 ; then
    fn_log "Mapping sample reads to spikein genome" 2>&1 | tee -a "${LOGFILE}"
    cmd="bowtie2 \
        --threads "${CPU}" \
        -x "${SPIKEIN_BASE}" \
        -1 "${SAMPLE_R1}" \
        -2 "${SAMPLE_R2}" \
        --un-conc-gz "${SAMPLE_NON_ALIGNED_CALIBRATION}".gz \
        > "${SAMPLE_ALIGNED_CALIBRATION}""
        fn_exec "${cmd}" "${LOGFILE}" 2>> "${LOGFILE}"
    
    fn_log "Mapping sample reads non-mapped on spikein genome to reference genome" 2>&1 | tee -a "${LOGFILE}"
    cmd="bowtie2 \
        --threads "${CPU}" \
        -x "${GENOME_BASE}" \
        -1 "${SAMPLE_NON_ALIGNED_CALIBRATION}".1.gz \
        -2 "${SAMPLE_NON_ALIGNED_CALIBRATION}".2.gz \
        > "${SAMPLE_NON_ALIGNED_CALIBRATION_ALIGNED_GENOME}""
    fn_exec "${cmd}" "${LOGFILE}" 2>> "${LOGFILE}"

    fn_log "Mapping sample reads non-mapped on reference genome to spikein genome" 2>&1 | tee -a "${LOGFILE}"
    cmd="bowtie2 \
        --threads "${CPU}" \
        -x "${SPIKEIN_BASE}" \
        -1 "${SAMPLE_NON_ALIGNED_GENOME}".1.gz \
        -2 "${SAMPLE_NON_ALIGNED_GENOME}".2.gz \
        > "${SAMPLE_NON_ALIGNED_GENOME_ALIGNED_CALIBRATION}""
    fn_exec "${cmd}" "${LOGFILE}" 2>> "${LOGFILE}"
fi

if test "${DO_INPUT}" == 0 ; then
    fn_log "Mapping input reads to reference genome" 2>&1 | tee -a "${LOGFILE}"
    cmd="bowtie2 \
        --threads "${CPU}" \
        -x "${GENOME_BASE}" \
        -1 "${INPUT_R1}" \
        -2 "${INPUT_R2}" \
        --un-conc-gz "${INPUT_NON_ALIGNED_GENOME}".gz \
        > "${INPUT_ALIGNED_GENOME}""
    fn_exec "${cmd}" "${LOGFILE}" 2>> "${LOGFILE}"

    if test "${DO_CALIBRATION}" == 0 ; then
        fn_log "Mapping input reads to spikein genome" 2>&1 | tee -a "${LOGFILE}"
        cmd="bowtie2 \
            --threads "${CPU}" \
            -x "${SPIKEIN_BASE}" \
            -1 "${INPUT_R1}" \
            -2 "${INPUT_R2}" \
            --un-conc-gz "${INPUT_NON_ALIGNED_CALIBRATION}".gz \
            > "${INPUT_ALIGNED_CALIBRATION}""
        fn_exec "${cmd}" "${LOGFILE}" 2>> "${LOGFILE}"

        fn_log "Mapping input reads non-mapped on spikein genome to reference genome" 2>&1 | tee -a "${LOGFILE}"
        cmd="bowtie2 \
            --threads "${CPU}" \
            -x "${GENOME_BASE}" \
            -1 "${INPUT_NON_ALIGNED_CALIBRATION}".1.gz \
            -2 "${INPUT_NON_ALIGNED_CALIBRATION}".2.gz \
            > "${INPUT_NON_ALIGNED_CALIBRATION_ALIGNED_GENOME}""
            fn_exec "${cmd}" "${LOGFILE}" 2>> "${LOGFILE}"

        fn_log "Mapping input reads non-mapped on reference genome to spikein genome" 2>&1 | tee -a "${LOGFILE}"
        cmd="bowtie2 \
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

if test "${DO_CALIBRATION}" == 1 ; then

    fn_log "Filtering sample bam file of reads mapped to reference genome" 2>&1 | tee -a "${LOGFILE}"
    cmd="samtools fixmate ${SAMTOOLS_OPTIONS} -m "${SAMPLE_ALIGNED_GENOME}" - \
        | samtools sort ${SAMTOOLS_OPTIONS} -T "${HASH}"_"${SAMPLE_ALIGNED_GENOME}"_sorting - \
        | samtools markdup ${SAMTOOLS_OPTIONS} ${REMOVE_DUPLICATES} -T "${HASH}"_"${SAMPLE_ALIGNED_GENOME}"_markdup - - \
        | samtools view ${SAMTOOLS_OPTIONS} ${FILTEROPTIONS} -1 -b - \
        | samtools sort ${SAMTOOLS_OPTIONS} -l 9 -T "${HASH}"_"${SAMPLE_ALIGNED_GENOME}"_sorting2 \
        -o "${SAMPLE_ALIGNED_GENOME_FILTERED}""
    fn_exec "${cmd}" "${LOGFILE}"
    cmd="samtools index -@ "${CPU}" "${SAMPLE_ALIGNED_GENOME_FILTERED}""
    fn_exec "${cmd}" "${LOGFILE}" 2>> "${LOGFILE}"

    if test "${DO_INPUT}" == 0 ; then

        fn_log "Filtering input bam file of reads mapped to reference genome" 2>&1 | tee -a "${LOGFILE}"
        cmd="samtools fixmate ${SAMTOOLS_OPTIONS} -m "${INPUT_ALIGNED_GENOME}" - \
            | samtools sort ${SAMTOOLS_OPTIONS} -T "${HASH}"_"${INPUT_ALIGNED_GENOME}"_sorting - \
            | samtools markdup ${SAMTOOLS_OPTIONS} ${REMOVE_DUPLICATES} -T "${HASH}"_"${INPUT_ALIGNED_GENOME}"_markdup - - \
            | samtools view ${SAMTOOLS_OPTIONS} ${FILTEROPTIONS} -1 -b - \
            | samtools sort ${SAMTOOLS_OPTIONS} -l 9 -T "${HASH}"_"${INPUT_ALIGNED_GENOME}"_sorting2 \
            -o "${INPUT_ALIGNED_GENOME_FILTERED}""
        fn_exec "${cmd}" "${LOGFILE}"
        cmd="samtools index -@ "${CPU}" "${INPUT_ALIGNED_GENOME_FILTERED}""
        fn_exec "${cmd}" "${LOGFILE}" 2>> "${LOGFILE}"

    fi

fi

if test "${DO_CALIBRATION}" == 0 ; then

    fn_log "Filtering sample bam file of reads unmapped on spikein genome and mapped to reference genome" 2>&1 | tee -a "${LOGFILE}"
    cmd="samtools fixmate ${SAMTOOLS_OPTIONS} -m "${SAMPLE_NON_ALIGNED_CALIBRATION_ALIGNED_GENOME}" - \
        | samtools sort ${SAMTOOLS_OPTIONS} -T "${HASH}"_"${SAMPLE_NON_ALIGNED_CALIBRATION_ALIGNED_GENOME}"_sorting - \
        | samtools markdup ${SAMTOOLS_OPTIONS} ${REMOVE_DUPLICATES} -T "${HASH}"_"${SAMPLE_NON_ALIGNED_CALIBRATION_ALIGNED_GENOME}"_markdup - - \
        | samtools view ${SAMTOOLS_OPTIONS} ${FILTEROPTIONS} -1 -b - \
        | samtools sort ${SAMTOOLS_OPTIONS} -l 9 -T "${HASH}"_"${SAMPLE_NON_ALIGNED_CALIBRATION_ALIGNED_GENOME}"_sorting2 \
        -o "${SAMPLE_NON_ALIGNED_CALIBRATION_ALIGNED_GENOME_FILTERED}""
    fn_exec "${cmd}" "${LOGFILE}"
    cmd="samtools index -@ "${CPU}" "${SAMPLE_NON_ALIGNED_CALIBRATION_ALIGNED_GENOME_FILTERED}""
    fn_exec "${cmd}" "${LOGFILE}" 2>> "${LOGFILE}"

    fn_log "Filtering sample bam file of reads unmapped on reference genome and mapped to spikein genome" 2>&1 | tee -a "${LOGFILE}"
    cmd="samtools fixmate ${SAMTOOLS_OPTIONS} -m "${SAMPLE_NON_ALIGNED_GENOME_ALIGNED_CALIBRATION}" - \
        | samtools sort ${SAMTOOLS_OPTIONS} -T "${HASH}"_"${SAMPLE_NON_ALIGNED_GENOME_ALIGNED_CALIBRATION}"_sorting - \
        | samtools markdup ${SAMTOOLS_OPTIONS} ${REMOVE_DUPLICATES} -T "${HASH}"_"${SAMPLE_NON_ALIGNED_GENOME_ALIGNED_CALIBRATION}"_markdup - - \
        | samtools view ${SAMTOOLS_OPTIONS} ${FILTEROPTIONS} -1 -b - \
        | samtools sort ${SAMTOOLS_OPTIONS} -l 9 -T "${HASH}"_"${SAMPLE_NON_ALIGNED_GENOME_ALIGNED_CALIBRATION}"_sorting2 \
        -o "${SAMPLE_NON_ALIGNED_GENOME_ALIGNED_CALIBRATION_FILTERED}""
    fn_exec "${cmd}" "${LOGFILE}"
    cmd="samtools index -@ "${CPU}" "${SAMPLE_NON_ALIGNED_GENOME_ALIGNED_CALIBRATION_FILTERED}""
    fn_exec "${cmd}" "${LOGFILE}" 2>> "${LOGFILE}"

    fn_log "Filtering input bam file of reads unmapped on spikein genome and mapped to reference genome" 2>&1 | tee -a "${LOGFILE}"
    cmd="samtools fixmate ${SAMTOOLS_OPTIONS} -m "${INPUT_NON_ALIGNED_CALIBRATION_ALIGNED_GENOME}" - \
        | samtools sort ${SAMTOOLS_OPTIONS} -T "${HASH}"_"${INPUT_NON_ALIGNED_CALIBRATION_ALIGNED_GENOME}"_sorting - \
        | samtools markdup ${SAMTOOLS_OPTIONS} ${REMOVE_DUPLICATES} -T "${HASH}"_"${INPUT_NON_ALIGNED_CALIBRATION_ALIGNED_GENOME}"_markdup - - \
        | samtools view ${SAMTOOLS_OPTIONS} ${FILTEROPTIONS} -1 -b - \
        | samtools sort ${SAMTOOLS_OPTIONS} -l 9 -T "${HASH}"_"${INPUT_NON_ALIGNED_CALIBRATION_ALIGNED_GENOME}"_sorting2 \
        -o "${INPUT_NON_ALIGNED_CALIBRATION_ALIGNED_GENOME_FILTERED}""
    fn_exec "${cmd}" "${LOGFILE}"
    cmd="samtools index -@ "${CPU}" "${INPUT_NON_ALIGNED_CALIBRATION_ALIGNED_GENOME_FILTERED}""
    fn_exec "${cmd}" "${LOGFILE}" 2>> "${LOGFILE}"

    fn_log "Filtering input bam file of reads unmapped on reference genome and mapped to spikein genome" 2>&1 | tee -a "${LOGFILE}"
    cmd="samtools fixmate ${SAMTOOLS_OPTIONS} -m "${INPUT_NON_ALIGNED_GENOME_ALIGNED_CALIBRATION}" - \
        | samtools sort ${SAMTOOLS_OPTIONS} -T "${HASH}"_"${INPUT_NON_ALIGNED_GENOME_ALIGNED_CALIBRATION}"_sorting - \
        | samtools markdup ${SAMTOOLS_OPTIONS} ${REMOVE_DUPLICATES} -T "${HASH}"_"${INPUT_NON_ALIGNED_GENOME_ALIGNED_CALIBRATION}"_markdup - - \
        | samtools view ${SAMTOOLS_OPTIONS} ${FILTEROPTIONS} -1 -b - \
        | samtools sort ${SAMTOOLS_OPTIONS} -l 9 -T "${HASH}"_"${INPUT_NON_ALIGNED_GENOME_ALIGNED_CALIBRATION}"_sorting2 \
        -o "${INPUT_NON_ALIGNED_GENOME_ALIGNED_CALIBRATION_FILTERED}""
    fn_exec "${cmd}" "${LOGFILE}"
    cmd="samtools index -@ "${CPU}" "${INPUT_NON_ALIGNED_GENOME_ALIGNED_CALIBRATION_FILTERED}""
    fn_exec "${cmd}" "${LOGFILE}" 2>> "${LOGFILE}"

fi

if test "${MODE}" == MNase ; then 
    fn_log "Further filtering sample bam file of reads mapped to reference genome for fragment size (70-250 bp)" 2>&1 | tee -a "${LOGFILE}"
    cmd="samtools view -@ "${CPU}" -h "${SAMPLE_ALIGNED_GENOME_FILTERED}" \
        | mawk '/^@/ || (sqrt(($9^2)) > 70 && sqrt(($9^2)) < 250)' \
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

if test "${DO_CALIBRATION}" == 1 ; then

    fn_log "Generating CPM track for ${SAMPLE_BASE}" 2>&1 | tee -a "${LOGFILE}"
    cmd="bamCoverage \
        --bam "${SAMPLE_ALIGNED_GENOME_FILTERED}" \
        --outFileName "${SAMPLE_RAW_TRACK}" \
        --binSize 1 \
        --numberOfProcessors "${CPU}" \
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
        --normalizeUsing 'None' \
        --scaleFactor "${SCALING_FACTOR}" \
        --skipNonCoveredRegions \
        --extendReads \
        --ignoreDuplicates"
    fn_exec "${cmd}" "${LOGFILE}" 2>> "${LOGFILE}"

fi

if test "${MODE}" == MNase ; then 

    fn_log "Generating CPM track for ${SAMPLE_BASE} filtered for reads >70bp & <250bp" 2>&1 | tee -a "${LOGFILE}"
    cmd="bamCoverage \
        --bam "${SAMPLE_ALIGNED_GENOME_FILTERED_READSIZE}" \
        --outFileName "${SAMPLE_READSIZE_TRACK}" \
        --binSize 1 \
        --numberOfProcessors "${CPU}" \
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
        --normalizeUsing CPM \
        --skipNonCoveredRegions \
        --extendReads \
        --ignoreDuplicates \
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

if test "${DO_PEAKS}" == 0 && test "${DO_INPUT}" == 1 ; then

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

## ------------------------------------------------------------------
## ------------------- CHECK NB OF READS ----------------------------
## ------------------------------------------------------------------

echo -e "---" >> "${LOGFILE}"
fn_log "NUMBER OF SEQUENCED FRAGMENTS: ${SAMPLE_R1}" >> "${LOGFILE}"
echo `fastqfastcnt "${SAMPLE_R1}" frags.` >> "${LOGFILE}"
echo -e "---" >> "${LOGFILE}"

if test "${DO_CALIBRATION}" == 0 ; then
    files="${SAMPLE_ALIGNED_GENOME} ${SAMPLE_NON_ALIGNED_CALIBRATION_ALIGNED_GENOME_FILTERED}" 
else
    files="${SAMPLE_ALIGNED_GENOME} ${SAMPLE_ALIGNED_GENOME_FILTERED}"
fi

for file in ${files}
do
    fn_log "MAPPING STATS: ${file}" >> "${LOGFILE}"
    cmd="samtools flagstat "${file}" >> "${LOGFILE}""
    fn_exec "${cmd}" "${LOGFILE}"
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
cp "${BASH_SOURCE}" "${OUTDIR}"/`date "+%y%m%d"`-"${HASH}"-script.txt
chmod -x "${OUTDIR}"/`date "+%y%m%d"`-"${HASH}"-script.txt

## -- List generated files
echo -e "---" 2>&1 | tee -a "${LOGFILE}"
fn_log "GENERATED FILES" >> "${LOGFILE}"
tree "${OUTDIR}" -P "*${HASH}*" --prune 2>&1 | tee -a "${LOGFILE}"

## -- List commands
echo -e "---" 2>&1 | tee -a "${LOGFILE}"
fn_log "COMMANDS EXECUTED" >> "${LOGFILE}"
grep '\[EXEC\]' "${LOGFILE}" | sed 's,.*\[EXEC\] ,,' | sed 's,| ,\\\n\t| ,g' > "${OUTDIR}"/`date "+%y%m%d"`-"${HASH}"-commands.txt
grep '\[EXEC\]' "${LOGFILE}" | sed 's,.*\[EXEC\] ,,' | sed 's,| ,\\\n\t| ,g' >> "${LOGFILE}"

## -- Done !
echo -e "---" 2>&1 | tee -a "${LOGFILE}"
fn_log "Pipeline achieved: `date`" 2>&1 | tee -a "${LOGFILE}"
echo -e "Do check the log file ("${LOGFILE}") to make sure everything went ok!"
echo -e "You can re-run the commands by running \`./"${OUTDIR}"/`date "+%y%m%d"`-"${HASH}"-commands.txt\`."

## -- Remove color decorators from log file
sed -i 's/\x1b\[[0-9;]*m//g' "${LOGFILE}"
