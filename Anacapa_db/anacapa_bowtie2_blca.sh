#!/bin/bash

### this script is run as follows
# ~/Anacapa_db/anacapa_bowtie2_blca.sh -o <out_dir_for_anacapa_QC_run> -d <database_directory> -u <hoffman_account_user_name>
OUT=""
DB=""
UN=""


while getopts "o:d:u:" opt; do
    case $opt in
        o) OUT="$OPTARG" # path to desired Anacapa output
        ;;
        d) DB="$OPTARG"  # path to Anacapa_db
        ;;
        u) UN="$OPTARG"  # need username for submitting sequencing job
        ;;
    esac
done

####################################script & software
# This pipeline was developed and written by Emily Curd (eecurd@g.ucla.edu), Jesse Gomer (jessegomer@gmail.com), Baochen Shi (biosbc@gmail.com), and Gaurav Kandlikar (gkandlikar@ucla.edu), and with contributions from Zack Gold (zack.j.gold@gmail.com), Rachel Turba (rturba@ucla.edu) and Rachel Meyer (rsmeyer@ucla.edu).
# Last Updated 11-18-2017
#
# The purpose of these script is to process raw fastq.gz files from an Illumina sequencing run and generate summarized taxonomic assignment tables for multiple metabarcoding targets.
#
# This script is currently designed to run on UCLA's Hoffman2 cluster.  Please adjust the code to work with your computing resources. (e.g. module / path names to programs, submitting jobs for processing if you have a cluster, etc)
#
# This script runs in two phases. The first is a QC and dada2 sequence dereplication, denoising, merging (if reads are paired) and chimera detection.  The second phase runs bowtie2 and our blowtie2 modified blca run_scripts.
#
######################################

# Need to make a script to make sure dependencies are properly configured

# location of the config and var files
source $DB/scripts/anacapa_vars_nextV.sh  # edit to change variables and parameters
source $DB/scripts/anacapa_config.sh # edit for proper configuration
source $DB/scripts/anacapa_qsub_templates.sh # edit to change qsub job file format

##load modules / software
${MODULE_SOURCE} # use if you need to load modules from an HPC
${FASTX_TOOLKIT} #load fastx_toolkit
${ANACONDA_PYTHON} #load anaconda/python2-4.2
${PERL} #load perl
${ATS} #load ATS, Hoffman2-specific module for managing submitted jobs.
date
###

###############################
# Make sure unassembled reads are still paired and submit dada2 jobs
###############################

echo "Assign taxonomy!: 1) submit bowtie2 and blca for the dada2 output for each metabarcode"
for j in `ls ${OUT}/`
do
 if [[ "${j}" != "QC" && "${j}" != "Run_info" ]]; # ignore non-metabarcode folders...
 then
    #make folders for the metabarcode specific output of dada2 and bowtie2
 	echo "${j}"
    # generate runlogs that you can submit at any time!
    printf "$(BOWTIE2_BLCA_PAIRED_TEMPLATE)" > ${OUT}/Run_info/hoffman2/run_scripts/${j}_bowtie2_blca_job.sh
    echo ''
    qsub ${OUT}/Run_info/hoffman2/run_scripts/${j}_bowtie2_blca_job.sh
    if [ "${LOCALMODE}" = "TRUE"  ]  # if you are running locally (no hoffman2) you can run these jobs one after the other.
    then
        echo "Running Dada2 inline"
        bash ${OUT}/Run_info/hoffman2/run_scripts/${j}_bowtie2_blca_job.sh
    fi
 fi
done
date
echo "if a bowtie2/blca job fails you can find the job submission file in ${OUT}/Run_info/hoffman2/run_scripts"
date
echo "good luck!"
