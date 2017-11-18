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


Rscript  --vanilla ${DB}/scripts/dada2_reverse_only.R ${MB} ${OUT} 
echo "moving on"



rm -r ${OUT}/unpaired_2/${MB}/filtered

########################################
# Run bowtie2 on merged paired reads
########################################

### end to end mode
bowtie2 -x ${DB}/${MB}/${MB}_bowtie2_database/${MB}_bowtie2_index  -f -U ${OUT}/dada2_out/unpaired_R/${MB}/nochim_reverse${MB}.fasta -S ${OUT}/bowtie2_runs/unpaired_R/${MB}/nochim_reverse_only${MB}_end_to_end.sam --no-hd --no-sq -D 20 -R 3 -N 0 -L 15  -i S,1,0.50 --end-to-end --no-unal -p 120 -k 100 --un ${OUT}/dada2_out/unpaired_R/${MB}/nochim_reverse${MB}_end_to_end_reject.fasta --omit-sec-seq

### local mode

bowtie2 -x ${DB}/${MB}/${MB}_bowtie2_database/${MB}_bowtie2_index  -f -U ${OUT}/dada2_out/unpaired_R/${MB}/nochim_reverse${MB}_end_to_end_reject.fasta -S ${OUT}/bowtie2_runs/unpaired_R/${MB}/nochim_reverse_only${MB}_local.sam --no-hd --no-sq -D 20 -R 3 -N 0 -L 15  -i S,1,0.50 --local --no-unal -p 120 -k 100 --un ${OUT}/dada2_out/unpaired_F/${MB}/nochim_reverse${MB}_end_to_end_and_local_reject.fasta --omit-sec-seq
