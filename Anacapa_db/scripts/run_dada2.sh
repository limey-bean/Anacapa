#! /bin/bash

### this script is run as follows
# sh ~/Anacapa_db/scripts/run_bowtie2_make_3_Sfolders.sh -o <out_dir> -d <database_directory> -m <metabarcode name>   
OUT=""
DB=""
MB=""
TYP=""

while getopts "o:d:m:t:" opt; do
    case $opt in
        o) OUT="$OPTARG" # path to desired Anacapa output
        ;;
        d) DB="$OPTARG"  # path to Anacapa_db
        ;;
        m) MB="$OPTARG"  # need username for submitting sequencing job
        ;;
        t) TYP="$OPTARG"  # need username for submitting sequencing job
        ;;
    esac
done
####################################script & software

##load module
source $DB/scripts/anacapa_vars_nextV.sh
source $DB/scripts/anacapa_config.sh

##### add the single and paired bowtie 2 files to different folders. Turn the following code into a for loop for the single bowtie2 reads
echo ${MB} length:
length_var=LENGTH_${MB}
length=${!length_var}
echo ${length}

##load module
${MODULE_SOURCE} # use if you need to load modules from an HPC

${BOWTIE2} #load bowtie2
${ANACONDA_PYTHON} # load anaconda python
${PYTHON} # load python


#### critical or the dependency 'RcppParallel' will not install
${R}
${GCC}




mkdir -p ${OUT}/${MB}/${MB}dada2_out/individual_out

Rscript  --vanilla ${DB}/scripts/dada2_unified_script.R ${MB} ${OUT} ${length} ${TYP}
echo "moving on"
