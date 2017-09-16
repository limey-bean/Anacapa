# Variables file for anacapa_release_V1		09-15-2017
# Developed and written by Emily Curd (eecurd@g.ucla.edu), Jesse Gomer (jessegomer@gmail.com), and Baochen Shi (biosbc@gmail.com), with contributions from Gaurav Kandlikar (gkandlikar@ucla.edu), Zack Gold (zack.j.gold@gmail.com), Rachel Turba (rturba@ucla.edu) and Rachel Meyer (rsmeyer@ucla.edu).
# Developed at UCLA for the University of California Conservation Consortium's CALeDNA Program
# Last Updated 9-15-2017

###########################
# Parameters
###########################

#####
# cutadapt
#####

# cutadapt primer + adapter trim step "QC: 1) Run cutadapt to removal primers + sequencing adapters, sort for length, and quality"
# a or A adapters are 3' adapters and g or G adapters are 5' adapters
# in this script they are anchored to the ends of the read and allow up to 30% mismatch including indels

###QC 1
ERROR_QC1=".3"															# max errors allowed in trim step oligonucleotide (by percent)
FILE_TYPE_QC1="fastq"													# fasta or fastq. 
F_ADAPT="file:${DB}/adapters_and_PrimAdapt_rc/g_Forward_adapter.txt"   # path to file with the forward primer + adapter "g"
Rrc_PRIM_ADAPT="file:${DB}/adapters_and_PrimAdapt_rc/a_Reverse_PrimAdapt_rc.txt"   # path to the reverse complemented Reverse_adapter "a" that matches the forward primer + adapter "g"
R_ADAPT="file:${DB}/adapters_and_PrimAdapt_rc/G_Reverse_adapter.txt"   # path to G_Reverse_adapter.txt primer + adapter "G"
Frc_PRIM_ADAPT="file:${DB}/adapters_and_PrimAdapt_rc/A_Forward_PrimAdapt_rc.txt"   # path to the reverse complemented forward primer + adapter "A" that matches the reverse primer + adapter "G"
MIN_LEN="100"
MIN_QUAL="30"

###QC 4
ERROR_QC4=".3"															# max errors allowed in trim step oligonucleotide (by percent)
FILE_TYPE_QC4="fastq"													# fasta or fastq. 
R_PRIM_RC="file:${DB}/primers/a_reverse_rc_primers.txt"   # path to file with the forward primer + adapter "g"

###Primer sort
ERROR_PS=".3"															# max errors allowed in trim step oligonucleotide (by percent)
FILE_TYPE_PS="fastq"													# fasta or fastq. 
F_PRIM="file:${DB}/primers/g_forward_primers.txt"   # path to file with the forward primer + adapter "g"
R_PRIM="file:${DB}/primers/G_reverse_primers.txt"

#####
# pear
#####
P_VAL="0.0001" 	# Specify a p-value for the statistical test. If the computed p-value of a possible assembly exceeds the specified p-value then the paired-end read will not be assembled. Valid options are: 0.0001, 0.001, 0.01, 0.05 and 1.0. Setting 1.0 disables the test. (default: 0.01)
THREADS="100" 	#Number of threads to use

