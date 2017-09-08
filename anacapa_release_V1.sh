#!/bin/bash

if [ $# -ne 5 ]; then
 echo "need 5 paras: <input_dir> <out_dir> <ref_dir> <hoffman_account_user_name> <uclust_percent>"; # for <uclust_percent> use proportion (e.g. 97% = .97)
 exit;

fi
####################################script & software
# This pipeline was developed and written by Emily Curd (eecurd@g.ucla.edu) and Baochen Shi (biosbc@gmail.com), with contributions from Gaurav Kandlikar (gkandlikar@ucla.edu), Zack Gold (zack.j.gold@gmail.com), and (Rachel Meyer (rsmeyer@ucla.edu).
# Last Updated 9-07-2017
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

script=./scripts
cutadapt=/u/local/apps/python/2.7.13/bin/cutadapt

##load module
source /u/local/Modules/default/init/bash

module load pear
module load qiime
module load fastx_toolkit
module load anaconda/python2-4.2
module load bowtie2
date
###


################################
# Preprocessing .fastq files
################################
echo " "
echo " "
echo "Preprocessing: 1) Generate an md5sum file"
md5sum $1/*fastq.gz > $1/*fastq.gz.md5sum  # need to add something to check the before and after md5sum
date
###
echo "Preprocessing: 2) Rename each file for Qiime compatibility"
###################################
suffix1=R1_001.fastq
suffix2=R2_001.fastq
###################################
mkdir -p $2
mkdir -p $2/fastq
###
for str in `ls $1/*_${suffix1}.gz`
do
 str1=${str%_S*}
 i=${str1#$1/}
 mod=${i//_/-} 
 cp $1/${i}_*_${suffix1}.gz $2/fastq/${mod}_1.fastq.gz
 cp $1/${i}_*_${suffix2}.gz $2/fastq/${mod}_2.fastq.gz
done
date
###

echo "Preprocessing: 3) Uncompress files"
gunzip $2/fastq/*
date
###

echo "Preprocessing: 4) Rename each read in each file to reflect the sample ID"
for str in `ls $2/fastq/*_1.fastq`
do
 str1=${str%_1*}
 FILE=${str1#$2/fastq/}
 sed -i -E "s/^@[[:alnum:]]*:([[:alnum:]]*:[[:alnum:]]*)/@${FILE}_:\1/g" ${str1}_1.fastq #### does not work properly but the following works if you change K00188 to match your fastq files:  sed -i "s/@K00188:/@${FILE}_:/g" ${str1}_1.fastq
 sed -i -E "s/^@[[:alnum:]]*:([[:alnum:]]*:[[:alnum:]]*)/@${FILE}_:\1/g" ${str1}_2.fastq #### does not work properly but the following works if you change K00188 to match your fastq files:  sed -i "s/@K00188:/@${FILE}_:/g" ${str1}_2.fastq
done
date
###


################################
# QC the preprocessed .fastq files
#############################


echo "QC: 1) Run PEAR to filter low quality reads, and assemble paired reads where possible. Unassembled paired reads, and discarded reads will be retained and analyzed"
mkdir -p $2/PEAR 
###
for str in `ls $2/fastq/*_1.fastq`
do
 str1=${str%_1*}
 j=${str1#$2/fastq/}
 ####change pear parameters below
 pear -f ${str1}_1.fastq -r ${str1}_2.fastq -o $2/PEAR/$j -q 30 -t 100 -j 100 >> $2/PEAR/pear-report.txt
done
###
mkdir -p $2/PEAR/assembled  $2/PEAR/unassembled  $2/PEAR/discarded
mv $2/PEAR/*.assembled.fastq     $2/PEAR/assembled
mv $2/PEAR/*.unassembled.*.fastq $2/PEAR/unassembled
mv $2/PEAR/*.discarded.fastq     $2/PEAR/discarded
date
###

echo "QC: 2) Run cutadapt to removal sequencing adapters, sort for length, and convert files to fasta format"
mkdir -p $2/fasta
mkdir -p $2/fasta/assembled/
mkdir -p $2/fasta/unassembled/
mkdir -p $2/fasta/discarded/
###
for str in `ls $2/PEAR/{assembled,discarded,unassembled}/*.fastq`
do
 str1=${str%.fastq}
 str2=${str%/*.fastq} 
 j=${str1#$2/PEAR/*/}
 k=${str2#$2/PEAR/}
 echo ${j} "..."
 ${cutadapt} -f fastq -a TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGGT -g GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG --minimum-length 125 -o $2/fasta/${k}/${j}.cutadapt.fasta ${str} >> $2/fasta/cutadapt-report.txt
 echo ${j} "...  check!"
done
date
###
###

echo "QC:3) Reverse complement the unassembled reverses for down stream analysis!"
for str in `ls $2/fasta/unassembled/*.unassembled.reverse.cutadapt.fasta`
do
 str1=${str%.cutadapt.fasta}
 j=${str1#$2/fasta/unassembled/}
 fastx_reverse_complement -i $2/fasta/unassembled/${j}.cutadapt.fasta -o $2/fasta/unassembled/${j}.rc.cutadapt.fasta -Q33
done
echo "reverse complementation.....  check!"
date
###
echo "concatenate reads by assembled, unassembled F or R, and discarded....!"
cat $2/fasta/assembled/*.assembled.cutadapt.fasta >> $2/fasta/all.assembled.fasta
cat $2/fasta/unassembled/*.unassembled.forward.cutadapt.fasta >> $2/fasta/all.unassembled.forward.fasta
cat $2/fasta/unassembled/*.unassembled.reverse.rc.cutadapt.fasta >> $2/fasta/all.unassembled.rc.reverse.fasta
cat $2/fasta/discarded/*.discarded.cutadapt.fasta >> $2/fasta/all.discarded.fasta
echo "concatenation complete"
date

###############################
#Split reads by metabarcode
###############################
date
echo "Split reads by metabarcode: 1) Use split_on_primer2.py to sort reads by primer set"
mkdir -p $2/primer_split
mkdir -p $2/primer_split/assembled
mkdir -p $2/primer_split/discarded
mkdir -p $2/primer_split/unassembled_forward
mkdir -p $2/primer_split/unassembled_reverse
###
cp $2/fasta/all.assembled.fasta  $2/primer_split/assembled/all.assembled_qc.fasta
cp $2/fasta/all.unassembled.forward.fasta  $2/primer_split/unassembled_forward/all.unassembled.forward_qc.fasta
cp $2/fasta/all.unassembled.rc.reverse.fasta  $2/primer_split/unassembled_reverse/all.unassembled.rc.reverse_qc.fasta
cp $2/fasta/all.discarded.fasta  $2/primer_split/discarded/all.discarded_qc.fasta
###
i=1
source $3/split_on_primer_files/split_primers_config.txt
for prim in $3/split_on_primer_files/primer_sort_file_*.txt
do
 echo "time to split round "$i #${primer_sort_file_1}
 for f in $2/primer_split/{assembled,unassembled_forward,unassembled_reverse,discarded}/
 do
  k=${f#$2/primer_split/}
   for str in $2/primer_split/${k}*qc.fasta
   do
    str1=${str%_qc.fasta}
    echo ${k}" reads...."
    python $3/scripts/Split_on_primers_2.py -f ${str1}_qc.fasta -p ${prim} ${parameters}
    echo ${k}" reads....   check!"
    cp ${f}unsorted.fasta ${f}unsorted_$i.fasta
    cp ${f}unsorted.fasta ${str1}_qc.fasta
    rm ${f}unsorted.fasta
   done
 done
 i=$((i+1))
 echo ${i}
done


#################
#Processes metabarcode reads for taxonomy
#################
echo "Process metabarcode reads for taxonomy: 1) submit pick_open_reference_otus job for each metabarcode"
mkdir -p $2/Qiime_open_ref/
mkdir -p $2/taxon_summaries
mkdir -p $2/Qiime_open_ref/err_log/
###
for str in `ls $2/primer_split/assembled/*_F.fasta`
do
 str1=${str%_F*}
 j=${str1#$2/primer_split/assembled/}
 	mkdir -p $2/taxon_summaries/${j}
 	mkdir -p $2/Qiime_open_ref/Params
 	printf "pick_otus:enable_rev_strand_match True \nassign_taxonomy:id_to_taxonomy_fp $3/${j}/${j}_taxonomy.txt \nassign_taxonomy:reference_seqs_fp $3/${j}/${j}_qiime.fasta \nassign_taxonomy:uclust_min_consensus_fraction .75 \nassign_taxonomy:uclust_max_accepts 50 \n" > "$2/Qiime_open_ref/Params/Params${j}_pick_open.txt" 
 	qsub -l highp,h_rt=48:00:00,h_data=60G -pe shared 2 -N pick${j}_open -cwd -m bea -o $2/Qiime_open_ref/err_log/${j}1P_4.out -e $2/Qiime_open_ref/err_log/${j}1P_4.err -M $4 $3/scripts/pick_open_otus_and_summ.sh ${j} $2 $3 $5
done
echo "check!"
date
echo "good_luck!"
