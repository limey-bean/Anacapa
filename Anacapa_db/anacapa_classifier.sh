#!/bin/bash

### this script is run as follows
# ~/Anacapa_db/anacapa_classifier.sh -o <out_dir_for_anacapa_QC_run> -d <database_directory> -u <hoffman_account_user_name> -l (if running in local mode, no argument required) -h (calls help screen, no argument required) -b <percent of missmatch alloweb between the qury and subject, for blca> -c <path to text file with bcc cut off summary values> -n <BLCA number of times to bootstrap> -x <Muscle alignment match score> -f <Muscle alignment mismatch score> -g <Muscle alignment gap penalty> -k path to custom HPC job submission header

OUT=""
DB=""
UN=""
B_VALUE=""
PER_MIN_LEN=""
HELP=""
LOCALMODE="FALSE"
BCC_CUT_OFF=""
BOOT=""
MATCH=""
MISMATCH=""
GAPP=""
HPC_HEADER=""

while getopts "o:d:u:b:p:h?:l?:c:n:x:f:g:k:" opt; do
    case $opt in
        o) OUT="$OPTARG" # path to desired Anacapa output
        ;;
        d) DB="$OPTARG"  # path to Anacapa_db
        ;;
        u) UN="$OPTARG"  # need username for submitting sequencing job
        ;;
        b) B_VALUE="$OPTARG"  # percent of missmatch allowed between the qury and subject for blca
        ;;
        p) PER_MIN_LEN="$OPTARG"  # minpercent of length of query relative to subject for blca
        ;;
        h) HELP="TRUE"  # calls help screen
        ;;
        l) LOCALMODE="TRUE"  # for running in local mode
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
        k) HPC_HEADER="$OPTARG" # path to custom HPC job submission header
        ;;
    esac
done

case "$DB" in
*/)
    DB=${DB%/}
    ;;
*)
    echo ""
    ;;
esac


case "$OUT" in
*/)
    OUT=${OUT%/}
    ;;
*)
    echo ""
    ;;
esac


if [ "${HELP}" = "TRUE" ]
then
  printf "<<< Anacapa: Taxonomic Assignment using Bowtie 2 and BLCA >>>\n\nThe purpose of this script is assign taxonomy to ASVs generated in the Sequence QC and ASV Parsing script. ASV files are mapped to CRUX reference libraries using Bowtie 2, taxonomy is assigned using BLCA, and summary tables are given. \n\n	For successful implementation \n		1. Make sure you have all of the dependencies and correct paths in the anacapa_config.sh file\n		2. Add the Metabarcode locus specific CRUX reference libraries to the Anacapa_db folder\n		3. All parameters can be modified using the arguments below.  Alternatively, all parameters can be altered in the anacapa_vars.sh folder\n\nArguments:\n- Required for either mode:\n	-o	path to output directory generated in the Sequence QC and ASV Parsing script\n	-d	path to Anacapa_db\n    \n- Optional:\n 	-u	If running on an HPC (e.g. UCLA's Hoffman2 cluster), this is your username: e.g. eecurd\n	-l	If running locally: -l  (no argument needed)\n	-b	Percent of missmatch allowed between the qury and subject for BLCA: 0.0 to 1.0 (default 0.8)\n	-p	Minimum percent of length of the subject reltive to the query for BLCA: 0.0 to 1.0 (default 0.8)\n -c	A list of BCC cut-off values to report taxonomy: \"0 to 100\" quotes required (default \"40 50 60 70 80 90 95\")	\n		The file must contain the following format: PERCENT=\"40 50 60 70 80 90 95 100\"\n		Where the value may differ but the PERCENT=\"values\" is required\n		see ~/Anacapa_db/scripts/BCC_default_cut_off.sh as an example\n	-n	BLCA number of times to bootstrap: integer value (default 100)\n	-m	Muscle alignment match score: default 1\n	-f	Muscle alignment mismatch score: default -2.5\n	-g	Muscle alignment gap penalty: default -2	\n	-k	Path to file with alternate HPC job submission parameters:  \n		default file = ~/Anacapa_db/scripts/Hoffman2_HPC_header.sh\n		modifiable template file = ~/Anacapa_db/scripts/anacapa_qsub_templates.sh\n		\n- Other:\n	-h	Shows program usage then quits\n\n\n\n "
  exit
else
  echo ""
fi


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

###Local or HPC mode check for username
if [[ "${LOCALMODE}" = "TRUE"  ]];
then
  echo "Running in local mode"
elif [[ "${LOCALMODE}" = "FALSE" && ! -z ${UN} ]];
then
  echo "Running in HPC mode"
elif [[ "${LOCALMODE}" = "FALSE" &&  -z ${UN} ]];
then
  echo "Running in HPC mode"
  echo "No username given..."
  echo ""
  exit
fi

#Check that user has all of the default flags set
if [[  -e ${OUT} && -e ${DB} ]];
then
  echo "Required Arguments Given"
  echo ""
else
  echo "Required Arguments Missing:"
  echo "check that you included arguments or correct paths for -o and -d"
  echo ""
  exit
fi


# location of the config and var files
source $DB/scripts/anacapa_vars.sh  # edit to change variables and parameters
source $DB/scripts/anacapa_config.sh # edit for proper configuration
DEF_BCC_CUT_OFF="${DB}/scripts/BCC_default_cut_off.sh"

if [[ -e ${HPC_HEADER} ]];
then
  source ${HPC_HEADER}
  echo ""
else
  source $DB/scripts/Hoffman2_HPC_header.sh
  HPC_HEADER_FILE="$DB/scripts/Hoffman2_HPC_header.sh"
  echo ""
fi

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
 if [[ "${j}" != "Run_info" ]]; # ignore non-metabarcode folders...
 then
    #make folders for the metabarcode specific output of dada2 and bowtie2
  echo ''
  echo "${j}"
    # generate runlogs that you can submit at any time!
    if [ "${LOCALMODE}" = "TRUE"  ]  # if you are running locally (no hoffman2) you can run these jobs one after the other.
    then
        echo "Running Bowtie 2 inline"
        printf "#!/bin/bash\n\n ${RUNNER} ${DB}/scripts/run_bowtie2_blca.sh -o ${OUT} -d ${DB} -m ${j} -l -b ${B_VALUE:=$BLCAB} -p ${PER_MIN_LEN:=$BLCAperMINlen} -c ${BCC_CUT_OFF:=$DEF_BCC_CUT_OFF} -n ${BOOT:=$BOOTSTRAP} -x ${MATCH:=$MUSMATCH} -f ${MISMATCH:=$MUSMISMATCH} -g ${GAPP:=$MUSGAPP} -k ${HPC_HEADER:=$HPC_HEADER_FILE}\n\necho _END_ [run_bowtie2_blca.sh]" > ${OUT}/Run_info/run_scripts/${j}_bowtie2_blca_job.sh
        chmod 755 ${OUT}/Run_info/run_scripts/*
        ${RUNNER} ${OUT}/Run_info/run_scripts/${j}_bowtie2_blca_job.sh
        date
    else
        printf "${B2_HEADER} \n\necho _BEGIN_ [run_bowtie2_blca_paired.sh]: `date`\n\n /bin/bash ${DB}/scripts/run_bowtie2_blca.sh -o ${OUT} -d ${DB} -m ${j} -u ${UN} -b ${B_VALUE:=$BLCAB} -p ${PER_MIN_LEN:=$BLCAperMINlen} -c ${BCC_CUT_OFF:=$DEF_BCC_CUT_OFF} -n ${BOOT:=$BOOTSTRAP} -x ${MATCH:=$MUSMATCH} -f ${MISMATCH:=$MUSMISMATCH} -g ${GAPP:=$MUSGAPP} -k ${HPC_HEADER:=$HPC_HEADER_FILE}\n\necho _END_ [run_bowtie2_blca.sh]" > ${OUT}/Run_info/run_scripts/${j}_bowtie2_blca_job.sh
        ${QUEUESUBMIT} ${OUT}/Run_info/run_scripts/${j}_bowtie2_blca_job.sh
        date
    fi
 fi
done

echo ''
echo "If a Bowtie 2 or BLCA step fails you can find the job script in ${OUT}/Run_info/run_scripts"

echo "good luck!"
