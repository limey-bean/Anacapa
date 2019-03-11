# function templates used to generate the job files that get submitted using qsub
# if your environment needs different arguments for your job scheduling system you can change these
# example usage:
# source anacapa_qsub_templates.sh
# $(B2_HEADER) >> some_file

DADA2_PAIRED_HEADER="#!/bin/bash\n#$ -q std.q\n#$ -cwd\n#$ -l mem_free=48G\n#$ -pe smp 1\n#$ -N paired_dada2\n#$ -M ${UN}\n#$ -o ${OUT}/Run_info/run_logs/paired.out\n#$ -e ${OUT}/Run_info/run_logs/paired.err"
DADA2_UNPAIRED_F_HEADER="#!/bin/bash\n#$ -q std.q\n#$ -cwd\n#$ -l mem_free=48G\n#$ -pe smp 1\n#$ -N unpaired_F_dada2\n#$ -M ${UN}\n#$ -o ${OUT}/Run_info/run_logs/unpaired_F.out\n#$ -e ${OUT}/Run_info/run_logs/unpaired_F.err"
DADA2_UNPAIRED_R_HEADER="#!/bin/bash\n#$ -q std.q\n#$ -cwd\n#$ -l mem_free=48G\n#$ -pe smp 1\n#$ -N unpaired_R_dada2\n#$ -M ${UN}\n#$ -o ${OUT}/Run_info/run_logs/unpaired_R.out\n#$ -e ${OUT}/Run_info/run_logs/unpaired_R.err"
B2_HEADER="#!/bin/bash\n#$ -q std.q\n#$ -cwd\n#$ -l mem_free=48G\n#$ -pe smp 1\n#$ -N bowtie2_blca \n#$ -M ${UN}\n#$ -o ${OUT}/Run_info/run_logs/bowtie2_blca${MB}.out \n#$ -e ${OUT}/Run_info/run_logs/bowtie2_blca${MB}.alignment_stats"
BLCA_HEADER="#!/bin/bash\n#$ -q std.q\n#$ -cwd\n#$ -l mem_free=48G\n#$ -pe smp 1\n#$ -N BLCA_${MB} \n#$ -M ${UN}\n#$ -o ${OUT}/Run_info/run_logs/blca${MB}.out\n#$ -e ${OUT}/Run_info/run_logs/blca${MB}.err"