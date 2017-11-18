#! /bin/bash

### this script is run as follows
# sh ~/Anacapa_db/scripts/run_bowtie2_make_3_Sfolders.sh -o <out_dir> -d <database_directory> -m <metabarcode name>   
OUT=""
DB=""
MB=""

while getopts "o:d:m:" opt; do
    case $opt in
        o) OUT="$OPTARG" # path to desired Anacapa output
        ;;
        d) DB="$OPTARG"  # path to Anacapa_db
        ;;
        m) MB="$OPTARG"  # need username for submitting sequencing job
        ;;
    esac
done
####################################script & software

##load module
source $DB/scripts/anacapa_vars.sh
source $DB/scripts/anacapa_config.sh

##load module
${MODULE_SOURCE} # use if you need to load modules from an HPC

${BOWTIE2} #load bowtie2
${ANACONDA_PYTHON} # load anaconda python

#### critical or the dependency 'RcppParallel' will not install
${R}
${GCC}


##### add the single and paired bowtie 2 files to different folders. Turn the following code into a for loop for the single bowtie2 reads
length_var=LENGTH_${MB}
length=${!length_var}
echo ${length}

Rscript  --vanilla ${DB}/scripts/dada2_paired.R ${MB} ${OUT} ${length}
echo "moving on"

rm -r ${OUT}/paired/${MB}/filtered

########################################
# Run bowtie2 on merged paired reads
########################################
### end to end mode

bowtie2 -x ${DB}/${MB}/${MB}_bowtie2_database/${MB}_bowtie2_index  -f -U ${OUT}/${MB}/dada2_bowtie2/paired/merged/nochim_merged${MB}.fasta -S ${OUT}/${MB}/dada2_bowtie2/paired/merged/nochim_merged${MB}_end_to_end.sam --no-hd --no-sq --very-sensitive --end-to-end --no-unal -p 120 -k 100 --un ${OUT}/${MB}/dada2_bowtie2/paired/merged/nochim_merged${MB}_end_to_end_reject.fasta --omit-sec-seq

### local mode

bowtie2 -x ${DB}/${MB}/${MB}_bowtie2_database/${MB}_bowtie2_index  -f -U ${OUT}/${MB}/dada2_bowtie2/paired/merged/nochim_merged${MB}_end_to_end_reject.fasta -S ${OUT}/${MB}/dada2_bowtie2/paired/merged/nochim_merged${MB}_local.sam --no-hd --no-sq --very-sensitive --local --no-unal -p 120 -k 100 --un ${OUT}/${MB}/dada2_bowtie2/paired/merged/nochim_merged${MB}_end_to_end_and_local_reject.fasta --omit-sec-seq

########################################
# Summarize bowtie2 runs on merged paired reads and append to dada2 output
########################################

python ${DB}/scripts/append_bowtie_to_summary.py ${OUT}/${MB}/dada2_bowtie2/paired/merged/nochim_merged${MB}.txt ${DB}/${MB}/${MB}_fasta_and_taxonomy/${MB}_taxonomy.txt ${OUT}/${MB}/dada2_bowtie2/paired/merged/nochim_merged${MB}_local.sam ${OUT}/${MB}/dada2_bowtie2/paired/merged/nochim_merged${MB}_end_to_end.sam

########################################
# Run bowtie2 on unmerged paired reads
########################################

### end to end mode
bowtie2 -x ${DB}/${MB}/${MB}_bowtie2_database/${MB}_bowtie2_index  -f -1 ${OUT}/${MB}/dada2_bowtie2/paired/unmerged/nochim_unmerged${MB}F.fasta -2 ${OUT}/${MB}/dada2_bowtie2/paired/unmerged/nochim_unmerged${MB}R.fasta -S ${OUT}/${MB}/dada2_bowtie2/paired/unmerged/nochim_unmerged${MB}_end_to_end.sam --no-hd --no-sq --very-sensitive --end-to-end --no-unal -p 120 -k 100 --fr --rf --no-mixed --un-conc ${OUT}/${MB}/dada2_bowtie2/paired/unmerged/nochim_unmerged${MB}_end_to_end_reject.fasta --omit-sec-seq --no-discordant

### local mode

bowtie2 -x ${DB}/${MB}/${MB}_bowtie2_database/${MB}_bowtie2_index -f -1 ${OUT}/${MB}/dada2_bowtie2/paired/unmerged/nochim_unmerged${MB}_end_to_end_reject.1.fasta -2 ${OUT}/${MB}/dada2_bowtie2/paired/unmerged/nochim_unmerged${MB}_end_to_end_reject.2.fasta -S ${OUT}/${MB}/dada2_bowtie2/paired/unmerged/nochim_unmerged${MB}_local.sam --no-hd --no-sq --very-sensitive --local --no-unal -p 120 -k 100 --un-conc ${OUT}/${MB}/dada2_bowtie2/paired/unmerged/nochim_unmerged${MB}_end_to_end_and_local_reject.fasta --omit-sec-seq --no-discordant

########################################
# Summarize bowtie2 runs on unmerged paired reads and append to dada2 output
########################################

python ${DB}/scripts/append_bowtie_to_summary.py ${OUT}/${MB}/dada2_bowtie2/paired/unmerged/nochim_unmerged${MB}.txt ${DB}/${MB}/${MB}_fasta_and_taxonomy/${MB}_taxonomy.txt ${OUT}/${MB}/dada2_bowtie2/paired/unmerged/nochim_unmerged${MB}_local.sam ${OUT}/${MB}/dada2_bowtie2/paired/unmerged/nochim_unmerged${MB}_end_to_end.sam 
