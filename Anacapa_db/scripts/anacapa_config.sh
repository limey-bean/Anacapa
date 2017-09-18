# Config file for Anacapa_release_V1		09-15-2017
# Developed and written by Emily Curd (eecurd@g.ucla.edu), Jesse Gomer (jessegomer@gmail.com), and Baochen Shi (biosbc@gmail.com), with contributions from Gaurav Kandlikar (gkandlikar@ucla.edu), Zack Gold (zack.j.gold@gmail.com), Rachel Turba (rturba@ucla.edu) and Rachel Meyer (rsmeyer@ucla.edu).
# Developed at UCLA for the University of California Conservation Consortium's CALeDNA Program
# Last Updated 9-15-2017


#############################
# Paths to programs / load programs
#############################

MODULE_SOURCE="source /u/local/Modules/default/init/bash" 	#if none, leave empty <- for HPC

### if not loading modules, need a way to not not run lines ~65 - 70 (module load bit...) and a way to run programs from load names not paths...

#load cutadapt
CUTADAPT="/u/local/apps/python/2.7.13/bin/cutadapt" 		#path to cutadapt binary. see documentation for how to obtain this script

#load pear
PEAR="module load pear"									#or what ever code is used to load pear in a bash shell, or path to pear

#load fastx_toolkit
FASTX_TOOLKIT="module load fastx_toolkit"				#or what ever code is used to load fastx_toolkit in a bash shell, or path to fastx_toolkit 

#load anaconda/python2-4.2
ANACONDA_PYTHON2="module load anaconda/python2-4.2"				#or what ever code is used to load anaconda/python2-4.2 in a bash shell, or path to anaconda/python2-4.2 

#load bowtie2
BOWTIE2="module load bowtie2"							#or what ever code is used to load bowtie2 in a bash shell, or path to bowtie2 

#load ATS
ATS="module load ATS"									#or what ever code is used to load ATS in a bash shell, or path to ATS.  ATS is a Hoffman2 module that allows the user to submit a job on the HPC from within a shell script 

#load Qiime
QIIME="module load qiime" 								#or what ever code is used to load qiime in a bash shell (e.g. on a mac it might be "macqiime")


