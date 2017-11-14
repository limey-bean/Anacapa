#! /bin/bash

### this script is run as follows
# sh ~/Anacapa_db/scripts/anacapa_release_20171110.sh -i <input_dir> -o <out_dir> -d <database_directory> -u <hoffman_account_user_name> -f <fasta file of forward primers> -r <fasta file of reverse primers> -a <adapter type (nextera or truseq)>  -t <illumina run type HiSeq or MiSeq>
IN=""
OUT=""
DB=""
UN=""
FP=""
RP=""
ADPT=""
ILLTYPE=""

while getopts "i:o:d:u:f:r:a:t:" opt; do
    case $opt in
        i) IN="$OPTARG" # path to raw .fastq.gz files
        ;;
        o) OUT="$OPTARG" # path to desired Anacapa output
        ;;
        d) DB="$OPTARG"  # path to Anacapa_db
        ;;
        u) UN="$OPTARG"  # need username for submitting sequencing job
        ;;
        f) FP="$OPTARG"  # need forward reads for cutadapt
        ;;
        r) RP="$OPTARG"  # need reverse reads for cutadapt
        ;;
        a) ADPT="$OPTARG"  # need adapter for cutadapt
        ;;
        t) ILLTYPE="$OPTARG"  #need to know trim params cutadapt
        ;;
    esac
done

####################################script & software
# This pipeline was developed and written by Emily Curd (eecurd@g.ucla.edu), Jesse Gomer (jessegomer@gmail.com), Baochen Shi (biosbc@gmail.com), and Gaurav Kandlikar (gkandlikar@ucla.edu), and with contributions from Zack Gold (zack.j.gold@gmail.com), Rachel Turba (rturba@ucla.edu) and Rachel Meyer (rsmeyer@ucla.edu).
# Last Updated 11-10-2017
#
# The purpose of this script is to process raw fastq.gz files from an Illumina sequencing and generate summarized taxonomic assignment tables for multiple metabarcoding targets.
#
# This script is currently designed to run on UCLA's Hoffman2 cluster.  Please adjust the code to work with your computing resources. (e.g. module / path names to programs, submitting jobs for processing if you have a cluster, etc)
#
# This script runs in two phases, the first is the qc phase that follows the anacapa_release_20171110.sh.  The second phase follows the run_dada2_bowtie2.sh scripts, and includes dada2 denoising, mergeing (if reads are paired) and chimera detection / bowtie2 sequence assignment phase.
#
######################################

# Need to make a script to make sure dependencies are properly configured

# location of the config and var files
source $DB/scripts/anacapa_vars.sh  # edit to change variables and parameters
source $DB/scripts/anacapa_config.sh # edit for proper configuration


##load modules / software
${MODULE_SOURCE} # use if you need to load modules from an HPC
${FASTX_TOOLKIT} #load fastx_toolkit
${ANACONDA_PYTHON} #load anaconda/python2-4.2
${PERL} #load perl
${ATS} #load ATS, Hoffman2 specific module for managing submitted jobs.
date
###

################################
# Preprocessing .fastq files
################################
echo " "
echo " "
echo "Preprocessing: 1) Generate an md5sum file"  # user can check for file corruption
md5sum ${IN}/*fastq.gz > ${IN}/*fastq.gz.md5sum  
date
###
echo "Preprocessing: 2) Rename each file for readability" # remove the additional and less relevant information in an illumna fasta file name 
###################################
suffix1=R1_001.fastq
suffix2=R2_001.fastq
###################################
mkdir -p ${OUT}
mkdir -p ${OUT}/fastq
###
for str in `ls ${IN}/*_${suffix1}.gz`
do
 str1=${str%_S*}
 i=${str1#${IN}/}
 mod=${i//_/-} 
 cp ${IN}/${i}_*_${suffix1}.gz ${OUT}/fastq/${mod}_1.fastq.gz
 cp ${IN}/${i}_*_${suffix2}.gz ${OUT}/fastq/${mod}_2.fastq.gz
done
date
###

echo "Preprocessing: 3) Uncompress files"
gunzip ${OUT}/fastq/*
date
###

################################
# QC the preprocessed .fastq files
#############################

echo "QC: 1) Run cutadapt to remove 5'sequncing adapters and 3'primers + sequencing adapters, sort for length, and quality."

# Generate cut adapt primer files -> merge reverse complemented primers with adapters for cutting 3'end sequencing past the end of the metabarcode region, and add cutadapt specific characters to primers and primer/adapter combos so that the appropriate ends of reads are trimmed
mkdir -p ${DB}/adapters_and_PrimAdapt_rc
mkdir -p ${DB}/primers
echo " "
echo "Generating Primer and Primer + Adapter files for for cutadapt steps.  Your adapter type is ${ADPT}."
python ${DB}/scripts/anacapa_format_primers_cutadapt.py ${ADPT} ${FP} ${RP} ${DB}

# now use the formated cutadapt primer file to trim fastq reads
mkdir -p ${OUT}/cutadapt_fastq
mkdir -p ${OUT}/cutadapt_fastq/untrimmed
mkdir -p ${OUT}/primer_sort/
###
for str in `ls ${OUT}/fastq/*_1.fastq`
do
 # first chop of the 5' adapter and 3' adapter and primer combo (reverse complemented)
 str1=${str%_*}
 j=${str1#${OUT}/fastq/}
 echo ${j} "..."
 ${CUTADAPT} -e ${ERROR_QC1} -f ${FILE_TYPE_QC1} -g ${F_ADAPT} -a ${Rrc_PRIM_ADAPT} -G ${R_ADAPT} -A ${Frc_PRIM_ADAPT} -o ${OUT}/cutadapt_fastq/untrimmed/${j}_Paired_1.fastq -p ${OUT}/cutadapt_fastq/untrimmed/${j}_Paired_2.fastq ${str1}_1.fastq ${str1}_2.fastq >> ${OUT}/cutadapt_fastq/cutadapt-report.txt
 # stringent quality fileter to get rid of the junky sequence at the ends - modify in config file
 fastq_quality_trimmer -t ${MIN_QUAL} -l ${MIN_LEN}  -i ${OUT}/cutadapt_fastq/untrimmed/${j}_Paired_1.fastq -o ${OUT}/cutadapt_fastq/${j}_qcPaired_1.fastq -Q33
 fastq_quality_trimmer -t ${MIN_QUAL} -l ${MIN_LEN}  -i ${OUT}/cutadapt_fastq/untrimmed/${j}_Paired_2.fastq -o ${OUT}/cutadapt_fastq/${j}_qcPaired_2.fastq -Q33
 # sort by metabarcode but run additional trimming.  It makes a differnce in merging reads in dada2.  Trimming varies based on seqeuncing platform.
 echo "forward..."
 if [ "${ILLTYPE}" == "MiSeq"  ]; # if MiSeq chop more off the end than if HiSeq - modify length in the vars file
 then
  ${CUTADAPT} -e ${ERROR_PS} -f ${FILE_TYPE_PS} -g ${F_PRIM}  -u -${MS_F_TRIM} -o ${OUT}/primer_sort/{name}_${j}_Paired_1.fastq  ${OUT}/cutadapt_fastq/${j}_qcPaired_1.fastq >> ${OUT}/cutadapt_fastq/cutadapt-report.txt
  echo "check"
  echo "reverse..."
  ${CUTADAPT} -e ${ERROR_PS} -f ${FILE_TYPE_PS} -g ${R_PRIM}  -u -${MS_R_TRIM} -o ${OUT}/primer_sort/{name}_${j}_Paired_2.fastq   ${OUT}/cutadapt_fastq/${j}_qcPaired_2.fastq >> ${OUT}/cutadapt_fastq/cutadapt-report.txt
  echo "check"
 else
  ${CUTADAPT} -e ${ERROR_PS} -f ${FILE_TYPE_PS} -g ${F_PRIM}  -u -${HS_F_TRIM} -o ${OUT}/primer_sort/{name}_${j}_Paired_1.fastq  ${OUT}/cutadapt_fastq/${j}_qcPaired_1.fastq >> ${OUT}/cutadapt_fastq/cutadapt-report.txt
  echo "check"
  echo "reverse..."
  ${CUTADAPT} -e ${ERROR_PS} -f ${FILE_TYPE_PS} -g ${R_PRIM}  -u -${HS_R_TRIM} -o ${OUT}/primer_sort/{name}_${j}_Paired_2.fastq   ${OUT}/cutadapt_fastq/${j}_qcPaired_2.fastq >> ${OUT}/cutadapt_fastq/cutadapt-report.txt
  echo "check"
 fi
 date
 echo ${j} "...  check!"
done
date
###

###############################
# Make sure unassembled reads are still paired
###############################
makedir -p ${OUT}/paired/
makedir -p ${OUT}/unpaired/

echo "Checking that Paired reads are still paired: 1) Use  Armin PEYMANN perl script (https://www.biostars.org/p/56171/) to make sure that unassembled reads are still paired"
for str in `ls ${OUT}/primer_sort/*_Paired_1.fastq`
do
 str1=${str%_Paired_1.fastq}
 j=${str1#${OUT}/primer_sort/}
 echo ${j} "..."
 perl ${DB}/scripts/check_paired.pl ${OUT}/primer_sort/${j}_Paired_1.fastq ${OUT}/primer_sort/${j}_Paired_2.fastq
 echo ${j} "...check!" 
done
date

###############################
# Move files to paired and unpaired folders
###############################
mkdir -p ${OUT}/paired/
mkdir -p ${OUT}/unpaired_1/
mkdir -p ${OUT}/unpaired_2/

echo "Move paired and unpaired files to the correct folders"
for str in `ls ${OUT}/primer_sort/*_Paired_1.fastq`
do
 str1=${str%_*_Paired_1.fastq}
 j=${str1#${OUT}/primer_sort/}
 echo ${j} "..."
 mkdir -p ${OUT}/paired/${j}
 mkdir -p ${OUT}/unpaired_1/${j}
 mkdir -p ${OUT}/unpaired_2/${j}
 cp ${OUT}/primer_sort/${j}_*_Paired_1_singletons.fastq ${OUT}/unpaired_1/${j}
 cp ${OUT}/primer_sort/${j}_*_Paired_2_singletons.fastq ${OUT}/unpaired_2/${j}
 cp ${OUT}/primer_sort/${j}_*_sorted.fastq ${OUT}/paired/${j}
 echo ${j} "...check!" 
done
date

echo "Remove primer_sort, cutadapt_fastq, and fastq folders"
rm -r ${OUT}/primer_sort/
rm -r ${OUT}/cutadapt_fastq
rm -r ${OUT}/fastq

###############################
# Submit jobs for the paired and unpaired reads for each barcode
###############################

### Make directories for the dada2 results files
mkdir -p ${OUT}/dada2_out
mkdir -p ${OUT}/dada2_out/paired
mkdir -p ${OUT}/dada2_out/paired/merged
mkdir -p ${OUT}/dada2_out/paired/unmerged
mkdir -p ${OUT}/dada2_out/unpaired_F
mkdir -p ${OUT}/dada2_out/unpaired_R

### Make directories for the bowtie2 results files
mkdir -p ${OUT}/bowtie2_runs/
mkdir -p ${OUT}/bowtie2_runs/paired
mkdir -p ${OUT}/bowtie2_runs/paired/merged
mkdir -p ${OUT}/bowtie2_runs/paired/unmerged
mkdir -p ${OUT}/bowtie2_runs/unpaired_F
mkdir -p ${OUT}/bowtie2_runs/unpaired_R
mkdir -p ${OUT}/dada2_bowtie2/

### Make directories for the runlogs, and runscripts
mkdir -p ${OUT}/dada2_bowtie2/runscripts
mkdir -p ${OUT}/dada2_bowtie2/runlogs
mkdir -p ${OUT}/taxon_summaries

###
echo "Process metabarcode reads for taxonomy: 1) submit bowtie2 read for each metabarcode"
for j in `ls ${OUT}/paired/`
do
 if [ "${j}" != "unknown"  ]; # ignore all of the unknown reads...
 then
     #make folders for the metabarcode specific output of dada2 and bowtie2
 	mkdir -p ${OUT}/taxon_summaries/${j}
	mkdir -p ${OUT}/bowtie2_runs/paired/merged/${j}
	mkdir -p ${OUT}/bowtie2_runs/paired/unmerged/${j}
	mkdir -p ${OUT}/bowtie2_runs/unpaired_F/${j}
	mkdir -p ${OUT}/bowtie2_runs/unpaired_R/${j}
	mkdir -p ${OUT}/dada2_out/paired/merged/${j}
	mkdir -p ${OUT}/dada2_out/paired/unmerged/${j}
	mkdir -p ${OUT}/dada2_out/unpaired_F/${j}
	mkdir -p ${OUT}/dada2_out/unpaired_R/${j}
    echo "${j}"
    # generate runlogs that you can submit at any time!
    printf "#!/bin/bash\n#$ -l h_rt=02:00:00,h_data=8G\n#$ -N paired_${j}_dada2_bowtie2\n#$ -cwd\n#$ -m bea\n#$ -M ${UN}\n#$ -o ${OUT}/dada2_bowtie2/runlogs/${j}_paired.out\n#$ -e ${OUT}/dada2_bowtie2/runlogs/${j}_paired.err \n\necho _BEGIN_ [run_dada2_bowtie2_paired.sh]: `date`\n\nsh ${DB}/scripts/run_dada2_bowtie2_paired.sh  -o ${OUT} -d ${DB} -m ${j}\n\necho _END_ [run_dada2_bowtie2_paired.sh]" >> ${OUT}/dada2_bowtie2/runscripts/${j}_dada2_bowtie2_paired_job.sh
    printf "#!/bin/bash\n#$ -l h_rt=02:00:00,h_data=8G\n#$ -N unpaired_F_${j}_dada2_bowtie2\n#$ -cwd\n#$ -m bea\n#$ -M ${UN}\n#$ -o ${OUT}/dada2_bowtie2/runlogs/${j}_unpaired_F.out\n#$ -e ${OUT}/dada2_bowtie2/runlogs/${j}_unpaired_F.err \n\necho _BEGIN_ [run_dada2_bowtie2_unpaired_F.sh]: `date`\n\nsh ${DB}/scripts/run_dada2_bowtie2_unpaired_F.sh  -o ${OUT} -d ${DB} -m ${j}\n\necho _END_ [run_dada2_bowtie2_unpaired_F.sh]" >> ${OUT}/dada2_bowtie2/runscripts/${j}_dada2_bowtie2_unpaired_F_job.sh
    printf "#!/bin/bash\n#$ -l h_rt=02:00:00,h_data=8G\n#$ -N unpaired_R_${j}_dada2_bowtie2\n#$ -cwd\n#$ -m bea\n#$ -M ${UN}\n#$ -o ${OUT}/dada2_bowtie2/runlogs/${j}_unpaired_R.out\n#$ -e ${OUT}/dada2_bowtie2/runlogs/${j}_unpaired_R.err \n\necho _BEGIN_ [run_dada2_bowtie2_unpaired_R.sh]: `date`\n\nsh ${DB}/scripts/run_dada2_bowtie2_unpaired_R.sh  -o ${OUT} -d ${DB} -m ${j}\n\necho _END_ [run_dada2_bowtie2_unpaired_R.sh]" >> ${OUT}/dada2_bowtie2/runscripts/${j}_dada2_bowtie2_unpaired_R_job.sh
    echo ''
    # submit jobs to run dada2 and bowtie2
    qsub ${OUT}/dada2_bowtie2/runscripts/${j}_dada2_bowtie2_paired_job.sh
    qsub ${OUT}/dada2_bowtie2/runscripts/${j}_dada2_bowtie2_unpaired_F_job.sh
    qsub ${OUT}/dada2_bowtie2/runscripts/${j}_dada2_bowtie2_unpaired_R_job.sh
 fi
done
echo "check!"
echo "if a dada2 / bowtie2 job fails you can find the job submission file in ${OUT}/dada2_bowtie2/runscripts"
date
echo "good_luck!"



