#!/bin/bash

### this script is run as follows
# ~/Anacapa_db/scripts/anacapa_bowtie2_blca.sh -o <out_dir_for_anacapa_QC_run> -d <database_directory> -m <metabarcode> -u <hoffman_account_user_name> -b <percent of missmatch alloweb between the qury and subject, for blca> -l (run in local mode, no argument required) -c <path to text file with bcc cut off summary values> -n <BLCA number of times to bootstrap> -x <Muscle alignment match score> -f <Muscle alignment mismatch score> -g <Muscle alignment gap penalty> -k {hpc header file}

OUT=""
DB=""
MB=""
UN=""
B_VALUE=""
PER_MIN_LEN=""
LOCALMODE="FALSE"
BCC_CUT_OFF=""
BOOT=""
MATCH=""
MISMATCH=""
GAPP=""
HPC_HEADER_FILE=""

while getopts "o:d:m:u:b:p:l?:c:n:x:f:g:k:" opt; do
    case $opt in
        o) OUT="$OPTARG" # path to desired Anacapa output
        ;;
        d) DB="$OPTARG"  # path to Anacapa_db
        ;;
        m) MB="$OPTARG"  # metabarcode name
        ;;
        u) UN="$OPTARG"  # need username for submitting sequencing job
        ;;
        b) B_VALUE="$OPTARG"  # percent sim req for match
        ;;
        p) PER_MIN_LEN="$OPTARG" # min percent length of subject relative to query
        ;;
        l) LOCALMODE="TRUE"  # need sername for submitting sequencing job
        ;;
        c) BCC_CUT_OFF="$OPTARG" #path to text file with bcc cut off summary values
        ;;
        n) BOOT="$OPTARG" # number of times to bootstrap for taxonomy assignment
        ;;
        x) MATCH="$OPTARG" # score for muscle alignment matches
        ;;
        f) MISMATCH="$OPTARG" # penalty for muscle alignment mismathes
        ;;
        g) GAPP="$OPTARG" # penalty for muscle alignment gaps
        ;;
        k) HPC_HEADER_FILE="$OPTARG" # path to the HPC header file
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
###Local or HPC mode check for username
if [[ "${LOCALMODE}" = "TRUE"  ]];
then
  echo "Running in local mode"
elif [[ "${LOCALMODE}" = "FALSE" && ! -z ${UN} ]];
then
  echo "Running in HPC mode"
fi

# location of the config and var files
source $DB/scripts/anacapa_vars.sh  # edit to change variables and parameters
source $DB/scripts/anacapa_config.sh # edit for proper configuration
source "${HPC_HEADER_FILE}"
source "${BCC_CUT_OFF}"

##load modules / software
${MODULE_SOURCE} # use if you need to load modules from an HPC
${FASTX_TOOLKIT} #load fastx_toolkit
${ANACONDA_PYTHON} #load anaconda/python2-4.2
${BOWTIE2}
${ATS} #load ATS, Hoffman2 specific module for managing submitted jobs.
${PYTHONWNUMPY}
date
###

######################################
# Submit single read and un merged paired read ASV's to Bowtie2
######################################
mkdir -p ${OUT}/${MB}/${MB}bowtie2_out

echo "${MB}"
echo "Run Bowtie2 on merged, forward, and reverse dada2 ASV fasta file"
# find best taxonomic hits for single reads in using bowtie2's global alignment mode (end-to-end) first.  It only considers full length alighnments, and recovers up to -k returns (default = 100). The reads not assigned in global mode are then run in local alighnment mode. -p is the number of threads --no-hd and --no-sq suppress the header lines --no-unal means do not add unaligned reads to the sam file (output). --very-sensitive is a preset option in bowtie2 is slow but designed to be more accurate and sensitive

list="merged forward reverse"

for str in ${list}
do
  if [[ -e "${OUT}/${MB}/${MB}dada2_out/nochim_${str}${MB}.fasta" ]];
  then
    echo ""
    echo "nochim_${str}${MB}.fasta exists"
    echo "end-to-end"
    bowtie2 -x ${DB}/${MB}/${MB}_bowtie2_database/${MB}_bowtie2_index  -f -U ${OUT}/${MB}/${MB}dada2_out/nochim_${str}${MB}.fasta -S ${OUT}/${MB}/${MB}bowtie2_out/single_read_${str}_${MB}_end_to_end.sam --no-hd --no-sq --very-sensitive --end-to-end --no-unal -p 120 -k 100 --un ${OUT}/${MB}/${MB}bowtie2_out/single_read_${str}_${MB}_end_to_end_reject.fasta
    date
    #single read local mode
    # reads that do not get hits in global alignment mode are run through local alignment.  This mode does not require the hit to be the entire length of the query.  This includes partial matches to full length references, or full matches to short references.
    echo "local"
    bowtie2 -x ${DB}/${MB}/${MB}_bowtie2_database/${MB}_bowtie2_index  -f -U ${OUT}/${MB}/${MB}bowtie2_out/single_read_${str}_${MB}_end_to_end_reject.fasta -S ${OUT}/${MB}/${MB}bowtie2_out/single_read_${str}_${MB}_local.sam --no-hd --no-sq --very-sensitive --local --no-unal -p 120 -k 100 --un ${OUT}/${MB}/${MB}bowtie2_out/single_read_${str}_${MB}_end_to_end_and_local_reject.fasta
    date
  else
    echo ""
    echo "nochim_${str}${MB}.fasta does not exist"
  fi
done
############# paired reads global and local

#unmerged pair reads global

echo "Run Bowtie2 on unmerged dada2 ASV fasta files"

# same logic as above but applied to paired end reads.
if [[ -e "${OUT}/${MB}/${MB}dada2_out/nochim_unmerged${MB}F.fasta" ]];
then
  echo ""
  echo "nochim_unmerged${MB}.fasta files exists"
  echo "end-to-end"
  bowtie2 -x ${DB}/${MB}/${MB}_bowtie2_database/${MB}_bowtie2_index  -f -1 ${OUT}/${MB}/${MB}dada2_out/nochim_unmerged${MB}F.fasta -2 ${OUT}/${MB}/${MB}dada2_out/nochim_unmerged${MB}R.fasta -S ${OUT}/${MB}/${MB}bowtie2_out/paired_unmerged_read_${MB}_end_to_end.sam --no-hd --no-sq --very-sensitive --end-to-end --no-unal -p 120 -k 100 --fr --rf --no-mixed --un-conc ${OUT}/${MB}/${MB}bowtie2_out/paired_unmerged_read_${MB}_end_to_end_reject.fasta --no-discordant
  #unmerged pair reads local
  echo "local"
  bowtie2 -x ${DB}/${MB}/${MB}_bowtie2_database/${MB}_bowtie2_index -f -1 ${OUT}/${MB}/${MB}bowtie2_out/paired_unmerged_read_${MB}_end_to_end_reject.1.fasta -2 ${OUT}/${MB}/${MB}bowtie2_out/paired_unmerged_read_${MB}_end_to_end_reject.2.fasta -S ${OUT}/${MB}/${MB}bowtie2_out/paired_unmerged_read_${MB}_local.sam --no-hd --no-sq --very-sensitive --local --no-unal -p 120 -k 100 --un-conc ${OUT}/${MB}/${MB}bowtie2_out/paired_unmerged_read_${MB}_end_to_end_and_local_reject.fasta --no-discordant
else
  echo ""
  echo "nochim_unmerged${MB}.fasta files do not exist"
fi

######################################
# concatenate ASV read tables dada2 ASV's
######################################
echo ""
echo "Concatenate all ASV site frequency tables and sam output"
mkdir -p ${OUT}/${MB}/${MB}_taxonomy_tables
# combine all dada2 tables for a given metabarcode, match samples in ASV tables.
# make a breif summary table
python ${DB}/scripts/merge_asv1.py ${OUT}/${MB}/${MB}dada2_out/nochim_forward${MB}.txt ${OUT}/${MB}/${MB}dada2_out/nochim_merged${MB}.txt ${OUT}/${MB}/${MB}dada2_out/nochim_reverse${MB}.txt ${OUT}/${MB}/${MB}dada2_out/nochim_unmerged${MB}.txt -o ${OUT}/${MB}/${MB}_taxonomy_tables/${MB}_ASV_taxonomy_brief.txt
# make a detailed sumary table.
python ${DB}/scripts/merge_asv.py ${OUT}/${MB}/${MB}dada2_out/nochim_forward${MB}.txt ${OUT}/${MB}/${MB}dada2_out/nochim_merged${MB}.txt ${OUT}/${MB}/${MB}dada2_out/nochim_reverse${MB}.txt ${OUT}/${MB}/${MB}dada2_out/nochim_unmerged${MB}.txt -o ${OUT}/${MB}/${MB}_taxonomy_tables/${MB}_ASV_taxonomy_detailed.txt

######################################
# concatenate bowtie2 tables and run blca
######################################

# enrich summary file with bowtie2 data
python ${DB}/scripts/append_bowtie_to_summary.py ${OUT}/${MB}/${MB}_taxonomy_tables/${MB}_ASV_taxonomy_detailed.txt ${OUT}/${MB}/${MB}bowtie2_out/


##### if running big files or have a cluster do not run in local mode
echo ""
echo "Run BLCA"
if [ "${LOCALMODE}" = "TRUE"  ];  # if you are running loally (no hoffman2) you can run these jobs one after the other.
then
   echo "Running BLCA inline"
   ### concat all of the sam files for blca
   cat ${OUT}/${MB}/${MB}bowtie2_out/*.sam > ${OUT}/${MB}/${MB}bowtie2_out/${MB}_bowtie2_all.sam

   ### blca

   echo "Run blca on sam output locally"
   python ${DB}/scripts/blca_from_bowtie.py -i ${OUT}/${MB}/${MB}bowtie2_out/${MB}_bowtie2_all.sam -r ${DB}/${MB}/${MB}_fasta_and_taxonomy/${MB}_taxonomy.txt -q ${DB}/${MB}/${MB}_fasta_and_taxonomy/${MB}_.fasta -b ${B_VALUE} -l ${PER_MIN_LEN:=$BLCAperMINlen} -p ${DB}/muscle -n ${BOOT:=$BOOTSTRAP} -m ${MATCH:=$MUSMATCH} -f ${MISMATCH:=$MUSMISMATCH} -g ${GAPP:=$MUSGAPP}

else
  for str in ${OUT}/${MB}/${MB}bowtie2_out/*.sam
  do
    str1=${str%*.sam}
    j=${str1#${OUT}/${MB}/${MB}bowtie2_out/}
    echo "${str}"
    [ $# -eq 0 ] && { echo "Usage: $0 filename"; exit 1; }
    [ ! -f "${str}" ] && { echo "Error: $0 file not found."; exit 2; }
    if [ -s "${str}" ]
    then
     # generate runlogs that you can submit at any time!
     printf "${BLCA_HEADER} \n/${RUNNER} ${DB}/scripts/run_blca.sh -o ${OUT} -d ${DB} -b ${B_VALUE:=$BLCAB} -l ${PER_MIN_LEN:=$BLCAperMINlen} -m ${MB} -s ${str} -n ${BOOT:=$BOOTSTRAP} -x ${MATCH:=$MUSMATCH} -f ${MISMATCH:=$MUSMISMATCH} -g ${GAPP:=$MUSGAPP} \n\n" > ${OUT}/Run_info/run_scripts/${j}_blca_job.sh
     echo ''
     ${QUEUESUBMIT} ${OUT}/Run_info/run_scripts/${j}_blca_job.sh
     echo "if a blca job(s) fails you can find the job submission file in ${OUT}/Run_info/run_scripts"
     echo "${str}.blca.complete" >> ${OUT}/${MB}/${MB}bowtie2_out/${MB}_complete_outfiles.txt
    else
      echo "${str} is empty, nothing to submit"
    fi
  done
  # need to check if the array is done before moving on to the next step
  filename="${OUT}/${MB}/${MB}bowtie2_out/${MB}_complete_outfiles.txt"
  filelines=`cat $filename`
  echo Start
  for line in $filelines ; do
    while ! [ -f ${line} ];
    do
      echo "files not ready"
      sleep 600
    done
  done
  # concatonate sam files
  cat ${OUT}/${MB}/${MB}bowtie2_out/*.blca.out >> ${OUT}/${MB}/${MB}bowtie2_out/${MB}_bowtie2_all.sam.blca.out
fi
date

# add taxonomy to site frequency specturm tables (brief and detailed)
echo ""
echo "Add blca taxonomy to the ASV site frequency table"
python ${DB}/scripts/append_blca_to_summary.py ${OUT}/${MB}/${MB}_taxonomy_tables/${MB}_ASV_taxonomy_brief.txt ${OUT}/${MB}/${MB}bowtie2_out/${MB}_bowtie2_all.sam.blca.out
python ${DB}/scripts/append_blca_to_summary.py ${OUT}/${MB}/${MB}_taxonomy_tables/${MB}_ASV_taxonomy_detailed.txt ${OUT}/${MB}/${MB}bowtie2_out/${MB}_bowtie2_all.sam.blca.out

mkdir -p ${OUT}/${MB}/${MB}_taxonomy_tables/Summary_by_percent_confidence

mkdir -p ${OUT}/Run_info/taxon_summary_logs

for per in ${PERCENT}
do
  mkdir -p ${OUT}/${MB}/${MB}_taxonomy_tables/Summary_by_percent_confidence/${per}
  echo "${per}"
  python ${DB}/scripts/reformat_summary_for_r.py ${OUT}/${MB}/${MB}_taxonomy_tables/${MB}_ASV_taxonomy_brief.txt  ${OUT}/${MB}/${MB}_taxonomy_tables/Summary_by_percent_confidence/${per}/${MB}_ASV_raw_taxonomy_${per}.txt ${per}
  ${R}
  Rscript --vanilla ${DB}/scripts/sum_blca_for_R_by_taxon.R  ${OUT}/${MB}/${MB}_taxonomy_tables/Summary_by_percent_confidence/${per}/${MB}_ASV_raw_taxonomy_${per}.txt ${MB} ${OUT}/${MB}/${MB}_taxonomy_tables/Summary_by_percent_confidence/${per}/${MB}_ASV_sum_by_taxonomy_${per}.txt 2>> ${OUT}/Run_info/taxon_summary_logs/${MB}_ASV_sum_by_taxonomy.log.txt
done
date
