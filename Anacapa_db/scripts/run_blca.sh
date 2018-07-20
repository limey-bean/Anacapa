#!/bin/bash

### this script is run as follows
# ~/Anacapa_db/scripts/anacapa_bowtie2_blca.sh -o <out_dir_for_anacapa_QC_run> -d <database_directory> -m <metabarcode> -s <sam file name> -n <BLCA number of times to bootstrap> -x <Muscle alignment match score> -f <Muscle alignment mismatch score> -g <Muscle alignment gap penalty>
OUT=""
DB=""
MB=""
FN=""
BOOT=""
MATCH=""
MISMATCH=""
GAPP=""
B_VALUE=""
PER_MIN_LEN=""

while getopts "o:d:m:s:n:x:f:g:b:l:" opt; do
    case $opt in
        o) OUT="$OPTARG" # path to desired Anacapa output
        ;;
        d) DB="$OPTARG"  # path to Anacapa_db
        ;;
        m) MB="$OPTARG"  # metabarcode name
        ;;
        s) FN="$OPTARG"  # input filename
        ;;
        n) BOOT="$OPTARG" # number of times to bootstrap for taxonomy assignment
        ;;
        x) MATCH="$OPTARG" # score for muscle alignment matches
        ;;
        f) MISMATCH="$OPTARG" # penalty for muscle alignment mismathes
        ;;
        g) GAPP="$OPTARG" # penalty for muscle alignment gaps
        ;;
        b) B_VALUE="$OPTARG"  # percent sim req for match
        ;;
        l) PER_MIN_LEN="$OPTARG" # min percent length of subject relative to query
        ;;
    esac
done

####################################script & software
# This pipeline was developed and written by Emily Curd (eecurd@g.ucla.edu), Jesse Gomer (jessegomer@gmail.com), Baochen Shi (biosbc@gmail.com), and Gaurav Kandlikar (gkandlikar@ucla.edu), and with contributions from Zack Gold (zack.j.gold@gmail.com), Rachel Turba (rturba@ucla.edu) and Rachel Meyer (rsmeyer@ucla.edu).
# Last Updated 11-18-2017
#
# The purpose of this script is to process raw fastq.gz files from an Illumina sequencing and generate summarized taxonomic assignment tables for multiple metabarcoding targets.
#
# This script is currently designed to run on UCLA's Hoffman2 cluster.  Please adjust the code to work with your computing resources. (e.g. module / path names to programs, submitting jobs for processing if you have a cluster, etc)
#
# This script runs in two phases. The first is the qc phase that follows the anacapa_release.  The second phase follows the run_dada2_bowtie2.sh scripts, and includes dada2 denoising, merging (if reads are paired) and chimera detection / bowtie2 sequence assignment phase.
#
######################################

# location of the config and var files
source $DB/scripts/anacapa_vars.sh  # edit to change variables and parameters
source $DB/scripts/anacapa_config.sh # edit for proper configuration


##load modules / software
${MODULE_SOURCE} # use if you need to load modules from an HPC
${FASTX_TOOLKIT} #load fastx_toolkit
${ANACONDA_PYTHON} #load anaconda/python2-4.2
${BOWTIE2}
${ATS} #load ATS, Hoffman2 specific module for managing submitted jobs.
${PYTHONWNUMPY}
date
###

### blca

echo "Run blca on sam output"
python ${DB}/scripts/blca_from_bowtie.py -i ${FN} -r ${DB}/${MB}/${MB}_fasta_and_taxonomy/${MB}_taxonomy.txt -q ${DB}/${MB}/${MB}_fasta_and_taxonomy/${MB}_.fasta -b ${B_VALUE} -l ${PER_MIN_LEN} -p ${DB}/muscle -n ${BOOT:=$BOOTSTRAP} -m ${MATCH:=$MUSMATCH} -f ${MISMATCH:=$MUSMISMATCH} -g ${GAPP:=$MUSGAPP}


cat > ${FN}.blca.complete
echo "wrote ${FN}.blca.complete file"
exit
