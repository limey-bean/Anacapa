#! /bin/bash

### this script is run as follows
# sh ~/Anacapa_db/scripts/anacapa_bowtie2_blca.sh -o <out_dir_for_anacapa_QC_run> -d <database_directory> -m <metabarcode> -l -u <hoffman_account_user_name>
OUT=""
DB=""
MB=""
UN=""


while getopts "o:d:m:u:l?" opt; do
    case $opt in
        o) OUT="$OPTARG" # path to desired Anacapa output
        ;;
        d) DB="$OPTARG"  # path to Anacapa_db
        ;;
        m) MB="$OPTARG"  # need username for submitting sequencing job
        ;;
        u) UN="$OPTARG"  # need username for submitting sequencing job
        ;;
        l) LOCALMODE="TRUE" #run dada2 locally (not on a cluster)
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
# This script runs in two phases, the first is the qc phase that follows the anacapa_release.  The second phase follows the run_dada2_bowtie2.sh scripts, and includes dada2 denoising, mergeing (if reads are paired) and chimera detection / bowtie2 sequence assignment phase.
#
######################################

# Need to make a script to make sure dependencies are properly configured

# location of the config and var files
source $DB/scripts/anacapa_vars_nextV.sh  # edit to change variables and parameters
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

######################################
# Submit single read and un merged paired read ASV's to Bowtie2
######################################
mkdir -p ${OUT}/${MB}/${MB}bowtie2_out
mkdir -p ${OUT}/${MB}/${MB}bowtie2_out/individual_out



echo "Run Bowtie2 on merged, forward, and reverse dada2 ASV fasta file"
# find best taxonomic hits for single reads in using bowtie2's global alighmnet mode (end-to-end) first.  I only considers full length alighnments, and recovers up to -k returns (default = 100). The reads not assigned in global mode are then run in local alighnment mode. -p is the number of threads --no-hd and --no-sq suppress the header lines --no-unal means do not add unaligned reads to the sam file (output). --very-sensitive is a preset option in bowtie2 is slow but designed to be more accurate and sensitive

list="merged forward reverse"

for str in ${list}
do
  bowtie2 -x ${DB}/${MB}/${MB}_bowtie2_database/${MB}_bowtie2_index  -f -U ${OUT}/${MB}/${MB}dada2_out/individual_out/nochim_${str}${MB}.fasta -S ${OUT}/${MB}/${MB}bowtie2_out/individual_out/single_read_${str}_${MB}_end_to_end.sam --no-hd --no-sq --very-sensitive --end-to-end --no-unal -p 120 -k 100 --un ${OUT}/${MB}/${MB}bowtie2_out/individual_out/single_read_${str}_${MB}_end_to_end_reject.fasta

#single read local mode
# reads that do not get hits in global alignment mode are run through local alignment.  This mode does not require the hit to be the entire length of the query.  This includes partial matches to full length references, or full matches to short references.
  bowtie2 -x ${DB}/${MB}/${MB}_bowtie2_database/${MB}_bowtie2_index  -f -U ${OUT}/${MB}/${MB}bowtie2_out/individual_out/single_read_${str}_${MB}_end_to_end_reject.fasta -S ${OUT}/${MB}/${MB}bowtie2_out/individual_out/single_read_${str}_${MB}_local.sam --no-hd --no-sq --very-sensitive --local --no-unal -p 120 -k 100 --un ${OUT}/${MB}/${MB}bowtie2_out/individual_out/single_read_${str}_${MB}_end_to_end_and_local_reject.fasta
done
############# paired reads global and local

#unmerged pair reads global

echo "Run Bowtie2 on unmerged dada2 ASV fasta files"

# same logic as above but applied to paired end reads.
bowtie2 -x ${DB}/${MB}/${MB}_bowtie2_database/${MB}_bowtie2_index  -f -1 ${OUT}/${MB}/${MB}dada2_out/individual_out/nochim_unmerged${MB}F.fasta -2 ${OUT}/${MB}/${MB}dada2_out/individual_out/nochim_unmerged${MB}R.fasta -S ${OUT}/${MB}/${MB}bowtie2_out/individual_out/paired_unmerged_read_${MB}_end_to_end.sam --no-hd --no-sq --very-sensitive --end-to-end --no-unal -p 120 -k 100 --fr --rf --no-mixed --un-conc ${OUT}/${MB}/${MB}bowtie2_out/individual_out/paired_unmerged_read_${MB}_end_to_end_reject.fasta --no-discordant

#unmerged pair reads local
bowtie2 -x ${DB}/${MB}/${MB}_bowtie2_database/${MB}_bowtie2_index -f -1 ${OUT}/${MB}/${MB}bowtie2_out/individual_out/paired_unmerged_read_${MB}_end_to_end_reject.1.fasta -2 ${OUT}/${MB}/${MB}bowtie2_out/individual_out/paired_unmerged_read_${MB}_end_to_end_reject.2.fasta -S ${OUT}/${MB}/${MB}bowtie2_out/individual_out/paired_unmerged_read_${MB}_local.sam --no-hd --no-sq --very-sensitive --local --no-unal -p 120 -k 100 --un-conc ${OUT}/${MB}/${MB}bowtie2_out/individual_out/paired_unmerged_read_${MB}_end_to_end_and_local_reject.fasta --no-discordant


######################################
# concatenate ASV read tables dada2 ASV's
######################################

echo "Concatenate all ASV site frequency tables and sam output"
mkdir -p ${OUT}/${MB}/${MB}_taxonomy_tables
# combine all dada2 tables for a given metabarcode, match samples in ASV tables.
# make a breif summary table
python ${DB}/scripts/merge_asv1.py ${OUT}/${MB}/${MB}dada2_out/individual_out/nochim_forward${MB}.txt ${OUT}/${MB}/${MB}dada2_out/individual_out/nochim_merged${MB}.txt ${OUT}/${MB}/${MB}dada2_out/individual_out/nochim_reverse${MB}.txt ${OUT}/${MB}/${MB}dada2_out/individual_out/nochim_unmerged${MB}.txt -o ${OUT}/${MB}/${MB}_taxonomy_tables/${MB}_ASV_taxonomy_brief.txt
# make a detailed sumary table.
python ${DB}/scripts/merge_asv.py ${OUT}/${MB}/${MB}dada2_out/individual_out/nochim_forward${MB}.txt ${OUT}/${MB}/${MB}dada2_out/individual_out/nochim_merged${MB}.txt ${OUT}/${MB}/${MB}dada2_out/individual_out/nochim_reverse${MB}.txt ${OUT}/${MB}/${MB}dada2_out/individual_out/nochim_unmerged${MB}.txt -o ${OUT}/${MB}/${MB}_taxonomy_tables/${MB}_ASV_taxonomy_detailed.txt

######################################
# concatenate bowtie2 tables and run blca
######################################

# enrich summary file with bowtie2 data
python ${DB}/scripts/append_bowtie_to_summary.py ${OUT}/${MB}/${MB}_taxonomy_tables/${MB}_ASV_taxonomy_detailed.txt ${OUT}/${MB}/${MB}bowtie2_out/individual_out/


##### if running big files or have a cluster do not run in local mode

echo "Run BLCA"
if [ "${LOCALMODE}" = "TRUE"  ]  # if you are running loally (no hoffman2) you can run these jobs one after the other.
then
   ### concat all of the sam files for blca
   cat ${OUT}/${MB}/${MB}bowtie2_out/individual_out/*.sam > ${OUT}/${MB}/${MB}bowtie2_out/${MB}_bowtie2_all.sam

   ### blca

   echo "Run blca on sam output locally"
   python ${DB}/scripts/blca_from_bowtie.py -i ${OUT}/${MB}/${MB}bowtie2_out/${MB}_bowtie2_all.sam -r ${DB}/${MB}/${MB}_fasta_and_taxonomy/${MB}_taxonomy.txt -q ${DB}/${MB}/${MB}_fasta_and_taxonomy/${MB}_.fasta -b 0.8 -p ${DB}/muscle

   ######################################
   # Add blca taxonomy to the asv table
   ######################################

   echo "Add blca taxonomy to the ASV site frequency table"
   python ${DB}/scripts/append_blca_to_summary.py ${OUT}/${MB}/${MB}_taxonomy_tables/${MB}_ASV_taxonomy_brief.txt ${OUT}/${MB}/${MB}bowtie2_blca_out/${MB}_bowtie2_all.sam.blca.out
   python ${DB}/scripts/append_blca_to_summary.py ${OUT}/${MB}/${MB}_taxonomy_tables/${MB}_ASV_taxonomy_detailed.txt ${OUT}/${MB}/${MB}bowtie2_blca_out/${MB}_bowtie2_all.sam.blca.out
else
  for str in ${OUT}/${MB}/${MB}bowtie2_out/individual_out/*.sam
  do
    str1=${str%*.sam}
    j=${str1#${OUT}/${MB}/${MB}bowtie2_out/individual_out/}
    echo "${str}"
    [ $# -eq 0 ] && { echo "Usage: $0 filename"; exit 1; }
    [ ! -f "${str}" ] && { echo "Error: $0 file not found."; exit 2; }
    if [ -s "${str}" ]
    then
 	   echo "${str}"
     # generate runlogs that you can submit at any time!
     printf "#!/bin/bash\n#$ -l highp,h_rt=200:00:00,h_data=20G\n#$ -N bowtie2_${j}_blca\n#$ -cwd\n#$ -m bea\n#$ -M ${UN} \n#$ -o ${OUT}/Run_info/hoffman2/run_logs/${j}_blca_$JOB_ID.out\n#$ -e ${OUT}/Run_info/hoffman2/run_logs/${j}_blca_$JOB_ID.err \n\necho _BEGIN_ [run_blca.sh]: `date`\n\nsh ${DB}/scripts/run_blca.sh -o ${OUT} -d ${DB} -m ${MB} -f ${str} \n\necho _END_ [run_blca.sh]" >> ${OUT}/Run_info/hoffman2/run_scripts/${j}_blca_job.sh
     echo ''
     qsub ${OUT}/Run_info/hoffman2/run_scripts/${j}_blca_job.sh
     echo "if a blca job(s) fails you can find the job submission file in ${OUT}/Run_info/hoffman2/run_scripts"
     echo "${str}.blca.complete" >> ${OUT}/${MB}/${MB}bowtie2_out/${MB}_complete_outfiles.txt
    else
      echo "${str} is empty, nothing to submit"
    fi
  done
fi
date

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
cat ${OUT}/${MB}/${MB}bowtie2_out/individual_out/*.blca.out >> ${OUT}/${MB}/${MB}bowtie2_out/${MB}_bowtie2_all.sam.blca.out
# add taxonomy to site frequency specturm tables (brief and detailed)
echo "Add blca taxonomy to the ASV site frequency table"
python ${DB}/scripts/append_blca_to_summary.py ${OUT}/${MB}/${MB}_taxonomy_tables/${MB}_ASV_taxonomy_brief.txt ${OUT}/${MB}/${MB}bowtie2_out/${MB}_bowtie2_all.sam.blca.out
python ${DB}/scripts/append_blca_to_summary.py ${OUT}/${MB}/${MB}_taxonomy_tables/${MB}_ASV_taxonomy_detailed.txt ${OUT}/${MB}/${MB}bowtie2_out/${MB}_bowtie2_all.sam.blca.out

mkdir -p ${OUT}/${MB}/${MB}_taxonomy_tables/Summary_by_percent_confidence

percent="60 70 80 90 95"

for per in ${percent}
do
  mkdir -p ${OUT}/${MB}/${MB}_taxonomy_tables/Summary_by_percent_confidence/${per}
  echo "${per}"
  python ${DB}/scripts/reformat_summary_for_r.py ${OUT}/${MB}/${MB}_taxonomy_tables/${MB}_ASV_taxonomy_brief.txt  ${OUT}/${MB}/${MB}_taxonomy_tables/Summary_by_percent_confidence/${per}/${MB}_ASV_raw_taxonomy_${per}.txt ${per}
  ${R}
  Rscript  --vanilla ${DB}/scripts/sum_blca_for_R_by_taxon.R  ${OUT}/${MB}/${MB}_taxonomy_tables/Summary_by_percent_confidence/${per}/${MB}_ASV_raw_taxonomy_${per}.txt ${MB} ${OUT}/${MB}/${MB}_taxonomy_tables/Summary_by_percent_confidence/${per}/${MB}_ASV_sum_by_taxonomy_${per}.txt
done
date
