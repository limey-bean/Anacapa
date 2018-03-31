# function templates used to generate the job files that get submitted using qsub
# if your environment needs different arguments for your job scheduling system you can change these
# example usage:
# source anacapa_qsub_templates.sh
# $(B2_HEADER) >> some_file

B2_HEADER="#!/bin/bash\n#$ -q std.q\n#$ -cwd\n#$ -l mem_free=48G\n#$ -pe smp 1\n#$ -N bowtie2_${j}_blca\n#$ -o ${OUT}/Run_info/run_logs/${j}_bowtie2_blca_$JOB_ID.out\n#$ -e ${OUT}/Run_info/run_logs/${j}_bowtie2_blca_$JOB_ID.err"
DADA2_PAIRED_HEADER="#!/bin/bash\n#$ -q std.q\n#$ -cwd\n#$ -l mem_free=48G\n#$ -pe smp 1\n#$ -N paired_${j}_dada2\n#$ -o ${OUT}/Run_info/run_logs/${j}_paired_$JOB_ID.out\n#$ -e ${OUT}/Run_info/run_logs/${j}_paired_$JOB_ID.err"
DADA2_UNPAIRED_F_HEADER="#!/bin/bash\n#$ -q std.q\n#$ -cwd\n#$ -l mem_free=48G\n#$ -pe smp 1\n#$ -N unpaired_F_${j}_dada2\n#$ -o ${OUT}/Run_info/run_logs/${j}_unpaired_F_$JOB_ID.out\n#$ -e ${OUT}/Run_info/run_logs/${j}_unpaired_F_$JOB_ID.err"
DADA2_UNPAIRED_R_HEADER="#!/bin/bash\n#$ -q std.q\n#$ -cwd\n#$ -l mem_free=48G\n#$ -pe smp 1\n#$ -N unpaired_R_${j}_dada2\n#$ -o ${OUT}/Run_info/run_logs/${j}_unpaired_R_$JOB_ID.out\n#$ -e ${OUT}/Run_info/run_logs/${j}_unpaired_R_$JOB_ID.err"
BLCA_HEADER="#!/bin/bash\n#$ -q std.q\n#$ -cwd\n#$ -l mem_free=48G\n#$ -pe smp 1\n#$ -N bowtie2_${j}_blca\n#$ -o ${OUT}/Run_info/run_logs/${j}_blca_$JOB_ID.out\n#$ -e ${OUT}/Run_info/run_logs/${j}_blca_$JOB_ID.err"
