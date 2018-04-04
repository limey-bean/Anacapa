#!/bin/bash

### this script is run as follows
# ~/Anacapa_db/scripts/run_bowtie2_make_3_Sfolders.sh -o <out_dir> -d <database_directory> -m <metabarcode name> -e <Minimum length for paired F and R reads to merge> -b <min ASV read count acceptable>
OUT=""
DB=""
MB=""
TYP=""
MIN_MERGE_LENGTH=""
MIN_ASV_ABUNDANCE=""

while getopts "o:d:m:t:e:b:" opt; do
    case $opt in
        o) OUT="$OPTARG" # path to desired Anacapa output
        ;;
        d) DB="$OPTARG"  # path to Anacapa_db
        ;;
        m) MB="$OPTARG"  # metabarcode name
        ;;
        t) TYP="$OPTARG"  # type of reads
        ;;
        e) MIN_MERGE_LENGTH="$OPTARG"  # File path to the minimum length reqired for paired F and R reads to overlap (length of the locus - primer size + 20 bp)
        ;;
        b) MIN_ASV_ABUNDANCE="$OPTARG"
        ;;
    esac
done
####################################script & software

##load module
source $DB/scripts/anacapa_vars.sh
source $DB/scripts/anacapa_config.sh
source ${MIN_MERGE_LENGTH}


##### add the single and paired bowtie 2 files to different folders. Turn the following code into a for loop for the single bowtie2 reads
length_var=LENGTH_${MB}
length=${!length_var}

mkdir -p ${OUT}/Run_info/dada2_out
mkdir -p ${OUT}/${MB}/${MB}dada2_out
echo ""
echo "Running dada2 on ${TYP} reads"

##load module
${MODULE_SOURCE} # use if you need to load modules from an HPC

${BOWTIE2} #load bowtie2
${ANACONDA_PYTHON} # load anaconda python
${PYTHON} # load python


#### critical or the dependency 'RcppParallel' will not install
${R} &> ${OUT}/Run_info/dada2_out/dada2_out_${TYP}
${GCC} &>> ${OUT}/Run_info/dada2_out/dada2_out_${TYP}

Rscript  --vanilla ${DB}/scripts/dada2_unified_script.R ${MB} ${OUT} ${length} ${TYP} ${MIN_ASV_ABUNDANCE} &>> ${OUT}/Run_info/dada2_out/dada2_out_${TYP}

echo "moving on"
