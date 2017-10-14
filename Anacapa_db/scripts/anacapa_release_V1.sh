#! /bin/bash

### this script is run as follows
# sh ~/Anacapa_db/scripts/anacapa_release_V1.sh -i <input_dir> -o <out_dir> -d <database_directory> -u <hoffman_account_user_name> -f <fasta file of forward primers> -r <fasta file of reverse primers> -a <adapter type (nextera or truseq)>  
IN=""
OUT=""
DB=""
UN=""
FP=""
RP=""
ADPT=""

while getopts "i:o:d:u:f:r:a:" opt; do
    case $opt in
        i) IN="$OPTARG" # path to raw .fastq.gz files
        ;;
        o) OUT="$OPTARG" # path to desired Anacapa output
        ;;
        d) DB="$OPTARG"  # path to Anacapa_db
        ;;
        u) UN="$OPTARG"  # need username for submitting sequencing job
        ;;
        f) FP="$OPTARG"  # need a warning if this is not a proportion from 0 to 1
        ;;
        r) RP="$OPTARG"  # need a warning if this is not a proportion from 0 to 1
        ;;
        a) ADPT="$OPTARG"  # need a warning if this is not a proportion from 0 to 1
        ;;
    esac
done

####################################script & software
# This pipeline was developed and written by Emily Curd (eecurd@g.ucla.edu), Jesse Gomer (jessegomer@gmail.com), and Baochen Shi (biosbc@gmail.com), with contributions from Gaurav Kandlikar (gkandlikar@ucla.edu), Zack Gold (zack.j.gold@gmail.com), Rachel Turba (rturba@ucla.edu) and Rachel Meyer (rsmeyer@ucla.edu).
# Last Updated 9-15-2017
#
# The purpose of this script is to process raw fastq.gz files from an Illumina sequencing and generate summarized taxonomic assignment tables for multiple metabardocing targets.
#
# This script is currently designed to run on UCLA's Hoffman2 cluster.  Please adjust the code to work with your computing resources. (e.g. module / path names to programs, submitting jobs for processing if you have a cluster, etc)
# If you are submitting this as a job, the following qsub parameters usually run right away and is more than enough to get the job done:
#	qsub -l highp,h_rt=24:00:00,h_data=60G -N metabarcoding_pipeline -cwd -m bea -o ./1P_4.out -e ./1P_4.err -M eecurd <path_to_metabarcoding_script> <input_dir> <out_dir> <ref_dir> <hoffman_account_user_name> <uclust_percent>
# for multiple cores
#	qsub -l highp,h_rt=24:00:00,h_data=60G -pe shared 2 -N metabarcoding_pipeline -cwd -m bea -o ./1P_4.out -e ./1P_4.err -M eecurd <path_to_metabarcoding_script> <input_dir> <out_dir> <ref_dir> <hoffman_account_user_name> <uclust_percent>

#
# The steps of the pipeline are as follows: 
# 	preprocess the .fastq files: 1) Generate an md5sum file, 2) Rename each file for Qiime compatibility, 3) Uncompress files, 4) Rename each read in each file to reflect the sample ID
#	QC the .fastq files: 1) Run PEAR to filter low quality reads, and assemble paired reads where possible. Unassembled paired reads, and discarded reads will be retained and analyzed, 2) Run cutadapt to removal sequencing adapters, 3) Convert fastq to fasta files
#	Split reads by metabarcode: 1) Use split_on_primer2.py to sort reads by primer set.  This step requires three iterations of primer splitting due to the degenerate nature of 16S (PPM) and CO1.  
#	Processes metabarcode reads for taxonomy: 1) submit array job for pick open refs.  Use pick_open_reference_otus to implement uclust at the user determined threshold. 97% or .97 is typical, but other options are possible.

#
# In order to run the script you need a scripts and reference directory that contains: 1) Scripts that are integral for running the main scripts, 2) your reference library folders that contain .fasta and the accompanying taxonomy.txt tiles.
#
######################################

# Need to make a script to make sure dependencies are properly configured

# location of the config and var files
source $DB/scripts/anacapa_vars.sh
source $DB/scripts/anacapa_config.sh


##load module
${MODULE_SOURCE} # use if you need to load modules from an HPC

${PEAR} #load pear
${FASTX_TOOLKIT} #load fastx_toolkit
${QIIME} #load qiime
${ANACONDA_PYTHON} #load anaconda/python2-4.2
${BOWTIE2} #load bowtie2
${PERL} #load perl
${ATS} #load ATS
date
###


################################
# Preprocessing .fastq files
################################
echo " "
echo " "
echo "Preprocessing: 1) Generate an md5sum file"
md5sum ${IN}/*fastq.gz > ${IN}/*fastq.gz.md5sum  # need to add something to check the before and after md5sum
date
###
echo "Preprocessing: 2) Rename each file for Qiime compatibility"
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


echo "Preprocessing: 4) Rename each read in each file to reflect the sample ID"
#this is an awk program that does the following:
#if the line is the first of 4 (starting count of 1) do the substitution
#otherwise just output the line
fastqrenamer="(NR % 4) == 1 {sub(/^@[[:alnum:]]+/, filename); print }
              (NR % 4) != 1 { print }"
for str in `ls ${OUT}/fastq/*_1.fastq`
do
 str1=${str%_1*}
 FILE=${str1#${OUT}/fastq/}
 awk -v filename="@${FILE}_" "${fastqrenamer}" ${str1}_1.fastq > ${str1}_1.fastq.tmp
 mv ${str1}_1.fastq.tmp ${str1}_1.fastq
 awk -v filename="@${FILE}_" "${fastqrenamer}" ${str1}_2.fastq > ${str1}_2.fastq.tmp
 mv ${str1}_2.fastq.tmp ${str1}_2.fastq
done
date
###

################################
# Generate cut adapt primer files
#############################
mkdir -p ${DB}/adapters_and_PrimAdapt_rc
mkdir -p ${DB}/primers
echo " "
echo "Generating Primer and Primer + Adapter files for for cutadapt steps.  If not using nextera indexes, please check the primer seqeunces"
python ${DB}/scripts/anacapa_format_primers_cutadapt.py ${ADPT} ${FP} ${RP} ${DB}

################################
# QC the preprocessed .fastq files
#############################

echo "QC: 1) Run cutadapt to remove 5'sequncing adapters and 3'primers + sequencing adapters, sort for length, and quality. Paired reads that pass QC will go to the paired file, Paires that did not pass filter will go to the unpaired file"
mkdir -p ${OUT}/cutadapt_fastq
mkdir -p ${OUT}/cutadapt_fastq/paired
mkdir -p ${OUT}/cutadapt_fastq/unpaired
###
for str in `ls ${OUT}/fastq/*_1.fastq`
do
 str1=${str%_*}
 j=${str1#${OUT}/fastq/}
 echo ${j} "..."
 ${CUTADAPT} -e ${ERROR_QC1}  -f ${FILE_TYPE_QC1} -g ${F_ADAPT} -a ${Rrc_PRIM_ADAPT} -G ${R_ADAPT} -A ${Frc_PRIM_ADAPT} --minimum-length ${MIN_LEN} -q ${MIN_QUAL}  --too-short-output ${OUT}/cutadapt_fastq/unpaired/${j}_unPaired_1.fastq --too-short-paired-output ${OUT}/cutadapt_fastq/unpaired/${j}_unPaired_2.fastq -o ${OUT}/cutadapt_fastq/paired/${j}_Paired_1.fastq -p ${OUT}/cutadapt_fastq/paired/${j}_Paired_2.fastq ${str1}_1.fastq ${str1}_2.fastq >> ${OUT}/cutadapt_fastq/cutadapt-report.txt
 echo ${j} "...  check!"
done
date
###


echo "QC: 2) Run cutadapt on the unpaired reads to remove 5'sequncing adapters and 3'primers + sequencing adapters, and sort for length"
###
for str in `ls ${OUT}/cutadapt_fastq/unpaired/*_1.fastq`
do
 str1=${str%_*}
 j=${str1#${OUT}/cutadapt_fastq/unpaired/}
 echo ${j} "..."
 ${CUTADAPT} -e ${ERROR_QC1}  -f ${FILE_TYPE_QC1} -g ${F_ADAPT} -a ${Rrc_PRIM_ADAPT} --minimum-length ${MIN_LEN} -o ${OUT}/cutadapt_fastq/unpaired/${j}_clean_unpaired_1.fastq ${str1}_1.fastq >> ${OUT}/cutadapt_fastq/cutadapt-report.txt
 ${CUTADAPT} -e ${ERROR_QC1}  -f ${FILE_TYPE_QC1} -g ${R_ADAPT} -a ${Frc_PRIM_ADAPT} --minimum-length ${MIN_LEN} -o ${OUT}/cutadapt_fastq/unpaired/${j}_clean_unpaired_2.fastq ${str1}_2.fastq >> ${OUT}/cutadapt_fastq/cutadapt-report.txt
 echo ${j} "...  check!"
done
date
###

echo "QC: 3) Run PEAR on paired reads to assemble reads where possible."
mkdir -p ${OUT}/PEAR 
###
for str in `ls ${OUT}/cutadapt_fastq/paired/*_1.fastq`
do
 str1=${str%_1*}
 j=${str1#${OUT}/cutadapt_fastq/paired/}
 pear -f ${str1}_1.fastq -r ${str1}_2.fastq -o ${OUT}/PEAR/${j} -p ${P_VAL} -j ${THREADS} >> ${OUT}/PEAR/pear-report.txt
done
###
mkdir -p ${OUT}/primer_sort/
mkdir -p ${OUT}/primer_sort/assembled  
mkdir -p ${OUT}/primer_sort/unassembled_F  
mkdir -p ${OUT}/primer_sort/unassembled_R 
mkdir -p ${OUT}/primer_sort/discarded_F
mkdir -p ${OUT}/primer_sort/discarded_R
###
cat ${OUT}/PEAR/*.assembled.fastq > ${OUT}/primer_sort/assembled/all.assembled.fastq
cat ${OUT}/PEAR/*.unassembled.forward.fastq > ${OUT}/primer_sort/unassembled_F/all.unassembled.F.fastq 
cat ${OUT}/PEAR/*.unassembled.reverse.fastq > ${OUT}/primer_sort/unassembled_R/all.unassembled.R.fastq 
cat ${OUT}/cutadapt_fastq/unpaired/*_clean_unpaired_1.fastq > ${OUT}/primer_sort/discarded_F/all.discarded_F.fastq
cat ${OUT}/cutadapt_fastq/unpaired/*_clean_unpaired_2.fastq > ${OUT}/primer_sort/discarded_R/all.discarded_R.fastq
date
###

echo "QC: 4) Use cutadapt to remove the 3' primers from the assembled reads."
${CUTADAPT} -e ${ERROR_QC4} -f ${FILE_TYPE_QC4} -a ${R_PRIM_RC} -o ${OUT}/primer_sort/assembled/all.clean_assembled.fastq  ${OUT}/primer_sort/assembled/all.assembled.fastq >> ${OUT}/cutadapt_fastq/cutadapt-report.txt
date
###

echo "QC: 5) Reverse complement the unassembled reverse reads for down stream analysis."
fastx_reverse_complement -i ${OUT}/primer_sort/unassembled_R/all.unassembled.R.fastq  -o ${OUT}/primer_sort/unassembled_R/all.unassembled.R_rc.fastq  -Q33
echo "reverse complementation.....  check!"
date
###

 
###############################
# Split reads by metabarcode
###############################
date
echo "Sort reads by metabarcode: 1) Use cutadapt to sort reads by primer set and trim primer sequences from reads"
###
###
echo "assembled..."
${CUTADAPT} -e ${ERROR_PS} -f ${FILE_TYPE_PS} -g ${F_PRIM} -o ${OUT}/primer_sort/assembled/{name}_all.clean_assembled.fasta  ${OUT}/primer_sort/assembled/all.assembled.fastq >> ${OUT}/primer_sort/assembled/cutadapt-report.txt
echo "check"
echo "unassembled forward..."
${CUTADAPT} -e ${ERROR_PS} -f ${FILE_TYPE_PS} -g ${F_PRIM} -o ${OUT}/primer_sort/unassembled_F/{name}_pear_unassembled_F.fasta ${OUT}/primer_sort/unassembled_F/all.unassembled.F.fastq >> ${OUT}/primer_sort/unassembled_F/cutadapt-report.txt
echo "check"
echo "unassembled reverse..."
${CUTADAPT} -e ${ERROR_PS} -f ${FILE_TYPE_PS} -g ${R_PRIM} -o ${OUT}/primer_sort/unassembled_R/{name}_pear_unassembled_R.fasta ${OUT}/primer_sort/unassembled_R/all.unassembled.R_rc.fastq >> ${OUT}/primer_sort/unassembled_R/cutadapt-report.txt
echo "check"
echo "discarded forward..."
${CUTADAPT} -e ${ERROR_PS} -f ${FILE_TYPE_PS} -g ${F_PRIM} -o ${OUT}/primer_sort/discarded_F/{name}_discarded_F.fasta  ${OUT}/primer_sort/discarded_F/all.discarded_F.fastq >> ${OUT}/primer_sort/discarded_F/cutadapt-report.txt
echo "check"
echo "discarded reverse..."
${CUTADAPT} -e ${ERROR_PS} -f ${FILE_TYPE_PS} -g ${R_PRIM} -o ${OUT}/primer_sort/discarded_R/{name}_discarded_R.fasta  ${OUT}/primer_sort/discarded_R/all.discarded_R.fastq >> ${OUT}/primer_sort/discarded_R/cutadapt-report.txt
echo "check"
date

###############################
# Make sure unassembled reads are still paired
###############################
date
echo "Checking that Paired reads are still paired: 1) Use  Armin PEYMANN perl script (https://www.biostars.org/p/56171/) to make sure that unassembled reads are still paired"
for str in `ls ${OUT}/primer_sort/unassembled_F/*_pear_unassembled_F.fasta`
do
 str1=${str%_pear_unassembled_F.fasta}
 j=${str1#${OUT}/primer_sort/unassembled_F/}
 echo ${j} "..."
 perl ${DB}/scripts/check_paired.pl ${OUT}/primer_sort/unassembled_F/${j}_pear_unassembled_F.fasta ${OUT}/primer_sort/unassembled_R/${j}_pear_unassembled_R.fasta
 echo ${j} "...check!"
done
date
echo "Merge unpaired unassembled reads with discarded reads"
for str in `ls ${OUT}/primer_sort/unassembled_F/*_pear_unassembled_F.fasta`
do
 str1=${str%_pear_unassembled_F.fasta}
 j=${str1#${OUT}/primer_sort/unassembled_F/}
 echo ${j} "..."
 cat ${OUT}/primer_sort/discarded_F/${j}_discarded_F.fasta ${OUT}/primer_sort/unassembled_F/${j}_pear_unassembled_F_singletons.fastq > ${OUT}/primer_sort/discarded_F/${j}_all.discarded_F.fasta
 cat ${OUT}/primer_sort/discarded_R/${j}_discarded_R.fasta ${OUT}/primer_sort/unassembled_R/${j}_pear_unassembled_R_singletons.fastq > ${OUT}/primer_sort/discarded_R/${j}_all.discarded_R.fasta
 echo ${j} "...check!"
done
date


#################
#Processes metabarcode reads for taxonomy
#################
echo "Process metabarcode reads for taxonomy: 1) submit bowtie2 read for each metabarcode"
mkdir -p ${OUT}/bowtie2_runs/
mkdir -p ${OUT}/bowtie2_runs/runlog
mkdir -p ${OUT}/taxon_summaries
###
for str in `ls ${OUT}/primer_sort/assembled/*_all.clean_assembled.fasta`
do
 str1=${str%_all.clean_assembled.fasta}
 j=${str1#${OUT}/primer_sort/assembled/}
 	mkdir -p ${OUT}/taxon_summaries/${j}
 	mkdir -p ${OUT}/bowtie2_runs/${j}
 	mkdir -p ${OUT}/bowtie2_runs/${j}/assembled
 	mkdir -p ${OUT}/bowtie2_runs/${j}/unassembled
 	mkdir -p ${OUT}/bowtie2_runs/${j}/discarded_F
 	mkdir -p ${OUT}/bowtie2_runs/${j}/discarded_R
	cp ${OUT}/primer_sort/discarded_F/${j}_all.discarded_F.fasta ${OUT}/bowtie2_runs/${j}/discarded_F/${j}_all.discarded_F.fasta
	cp ${OUT}/primer_sort/discarded_R/${j}_all.discarded_R.fasta ${OUT}/bowtie2_runs/${j}/discarded_R/${j}_all.discarded_R.fasta
	cp ${OUT}/primer_sort/assembled/${j}_all.clean_assembled.fasta ${OUT}/bowtie2_runs/${j}/assembled/${j}_all.clean_assembled.fasta
    cp ${OUT}/primer_sort/unassembled_F/${j}_pear_unassembled_F_sorted.fastq ${OUT}/bowtie2_runs/${j}/unassembled/${j}_pear_unassembled_F.fasta
    cp ${OUT}/primer_sort/unassembled_R/${j}_pear_unassembled_R_sorted.fastq ${OUT}/bowtie2_runs/${j}/unassembled/${j}_pear_unassembled_R.fasta
 	qsub -l highp,h_rt=6:00:00,h_data=12G  -N pick${j}_open -cwd -m bea -o ${OUT}/bowtie2_runs/runlog/{j}.out -e ${OUT}/bowtie2_runs/runlog/${j}.err -M ${UN} ${DB}/scripts/run_bowtie2_make_3_Sfolders.sh  -o ${OUT} -d ${DB} -n ${str}
done
echo "check!"
date
echo "good_luck!"
