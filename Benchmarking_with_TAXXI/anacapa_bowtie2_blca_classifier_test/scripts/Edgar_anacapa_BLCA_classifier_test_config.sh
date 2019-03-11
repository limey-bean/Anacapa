# Edgar_anacapa_BLCA_classifier_test_config.sh contains the paths to programs required to run the Edgar_anacapa_BLCA_classifier_test_hoffman2.sh scripts

#############################
# Paths to programs / load programs
#############################

MODULE_SOURCE="source /u/local/Modules/default/init/bash" 	#if none, leave empty <- for HPC

### if not loading modules, need a way to not run lines ~65 - 70 (module load bit...) and a way to run programs from load names not paths...

#load qiime and python 2.3.7 at same time
QIIME="module load qiime/1.8.0"				#or whatever is used to load qiime 1.8 in a bash shell

#load bowtie2
BOWTIE2="module load bowtie2"							# version 2.3.4 or what ever code is used to load bowtie2 in a bash shell, or path to bowtie2

#load ATS
ATS="module load ATS"									#or what ever code is used to load ATS in a bash shell, or path to ATS.  ATS is a Hoffman2 module that allows the user to submit a job on the HPC from within a shell script
