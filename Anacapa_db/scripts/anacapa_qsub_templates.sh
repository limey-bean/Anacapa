# function templates used to generate the job files that get submitted using qsub
# if your environment needs different arguments for your job scheduling system you can change these
# example usage:
# source anacapa_qsub_templates.sh
# $(BOWTIE2_BLCA_PAIRED_TEMPLATE) >> some_file

BOWTIE2_BLCA_PAIRED_TEMPLATE() {
  echo "#!/bin/bash\n#$ -q std.q\n#$ -cwd\n#$ -l mem_free=48G\n#$ -pe smp 1\n#$ -N bowtie2_${j}_blca\n#$ -o ${OUT}/Run_info/hoffman2/run_logs/${j}_bowtie2_blca_$JOB_ID.out\n#$ -e ${OUT}/Run_info/hoffman2/run_logs/${j}_bowtie2_blca_$JOB_ID.err \n\necho _BEGIN_ [run_bowtie2_blca_paired.sh]: `date`\n\n${RUNNER} ${DB}/scripts/run_bowtie2_blca.sh  -o ${OUT} -d ${DB} -m ${j} -u ${UN}\n\necho _END_ [run_bowtie2_blca.sh]"
}

DADA2_PAIRED_TEMPLATE() {
  echo "#!/bin/bash\n#$ -q std.q\n#$ -cwd\n#$ -l mem_free=48G\n#$ -pe smp 1\n#$ -N paired_${j}_dada2\n#$ -o ${OUT}/Run_info/hoffman2/run_logs/${j}_paired_$JOB_ID.out\n#$ -e ${OUT}/Run_info/hoffman2/run_logs/${j}_paired_$JOB_ID.err\n\necho _BEGIN_ [run_dada2_bowtie2_paired.sh]: `date`\n\n${RUNNER} ${DB}/scripts/run_dada2.sh  -o ${OUT} -d ${DB} -m ${j} -t paired\n\necho _END_ [run_dada2_paired.sh]"

}
DADA2_PAIRED_F_TEMPLATE() {
  echo "#!/bin/bash\n#$ -q std.q\n#$ -cwd\n#$ -l mem_free=48G\n#$ -pe smp 1\n#$ -N unpaired_F_${j}_dada2\n#$ -o ${OUT}/Run_info/hoffman2/run_logs/${j}_unpaired_F_$JOB_ID.out\n#$ -e ${OUT}/Run_info/hoffman2/run_logs/${j}_unpaired_F_$JOB_ID.err \n\necho _BEGIN_ [run_dada2_bowtie2_unpaired_F.sh]: `date`\n\n${RUNNER} ${DB}/scripts/run_dada2.sh  -o ${OUT} -d ${DB} -m ${j} -t forward\n\necho _END_ [run_dada2_unpaired_F.sh]"
}

DADA2_PAIRED_R_TEMPLATE() {
  echo "#!/bin/bash\n#$ -q std.q\n#$ -cwd\n#$ -l mem_free=48G\n#$ -pe smp 1\n#$ -N unpaired_R_${j}_dada2\n#$ -o ${OUT}/Run_info/hoffman2/run_logs/${j}_unpaired_R_$JOB_ID.out\n#$ -e ${OUT}/Run_info/hoffman2/run_logs/${j}_unpaired_R_$JOB_ID.err \n\necho _BEGIN_ [run_dada2_bowtie2_unpaired_R.sh]: `date`\n\n${RUNNER} ${DB}/scripts/run_dada2.sh  -o ${OUT} -d ${DB} -m ${j} -t reverse\n\necho _END_ [run_dada2_unpaired_R.sh]"
}

BLCA_TEMPLATE() {
  echo "#!/bin/bash\n#$ -q std.q\n#$ -cwd\n#$ -l mem_free=48G\n#$ -pe smp 1\n#$ -N bowtie2_${j}_blca\n#$ -o ${OUT}/Run_info/hoffman2/run_logs/${j}_blca_$JOB_ID.out\n#$ -e ${OUT}/Run_info/hoffman2/run_logs/${j}_blca_$JOB_ID.err \n\necho _BEGIN_ [run_blca.sh]: `date`\n\n${RUNNER} ${DB}/scripts/run_blca.sh -o ${OUT} -d ${DB} -m ${MB} -f ${str} \n\necho _END_ [run_blca.sh]"
}