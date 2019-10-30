#!/user/bin/env/ python

# Written by Emily Curd (eecurd@g.ucla.edu), with help from https://www.biostars.org/p/14614/
# for the University of California Conservation Consortium's CALeDNA Program
# Reverse complement primers
# python <path to script> <adapter_type> <input forward primer file> <input reverse primer file> <out put file path>

import sys
F_infile = open(sys.argv[2], "r") #fasta of forward primers
R_infile = open(sys.argv[3], "r") #fasta of reverse primers

out_path = sys.argv[4] +'/'
prim_adapt = out_path
prim = out_path

adapt = sys.argv[1]
next = "nextera"
true = "truseq"
neb = "NEBnext"

if adapt == next:
   adapter_Frc='ACACCTGTCTCTTATACACATCTGACGCTGCCGACGA'
   adapter_Rrc='CTGTCTCTTATACACATCTCCGAGCCCACGAGA'
elif adapt == true:
   adapter_Frc='AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT'
   adapter_Rrc='AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'
elif adapt == neb:
    adapter_Frc='GATCGGAAGAGCACACGTCTGAACTCCAGTCAC'
    adapter_Rrc='GATCGGAAGAGCACACGTCTGAACTCCAGTCAC'
else:
   print "\n warning \n not nextera or truseq add adater to anacapa_format_primers_cutadapt.py"


nuc_dict = {'A':'T','T':'A','U':'A','G':'C','C':'G','Y':'R','R':'Y','S':'S','W':'W','K':'M','M':'K','B':'V','D':'H','H':'D','V':'B','N':'N','I':'I','a':'T','t':'A','u':'A','g':'C','c':'G','y':'R','r':'Y','s':'S','w':'W','k':'M','m':'K','b':'V','d':'H','h':'D','v':'B','n':'N'}

def rComp(read):
    rc = ''
    for i in range(len(read) - 1,-1,-1):
        rc += nuc_dict[read[i]]

    return rc

####################################
# Process the F reads
####################################

### make regular primers with ^seq -> g
outfile = open(prim + "g_forward_primers.txt", "w+") # forwards with ^seq
header = ''
seq = ''
for line in F_infile:
    if line[0] == ">":
    	header = line.strip()
        outfile.write(header + "\n")
    else:
    	seq = line.strip()
    	outfile.write("^" + seq + "\n")
outfile.close()
F_infile.close()

### make reverse complement primers with seq$ -> a
F_infile = open(sys.argv[2], "r")
outfile = open(prim + "A_forward_rc_primers.txt", "w+") # forwards with ^seq
header = ''
seq = ''
for line in F_infile:
    if line[0] == ">":
    	header = line.strip()
        outfile.write(header + "_rc" + "\n")
    else:
    	seq = line.strip()
    	outfile.write(rComp(seq)+ "$" + "\n")
outfile.close()
F_infile.close()

### make adapter reverse complement primers with seq+adpter+$ ->
F_infile = open(sys.argv[2], "r")
outfile = open(prim_adapt + "A_Forward_PrimAdapt_rc.txt", "w+") # forwards with ^seq
header = ''
seq = ''
for line in F_infile:
    if line[0] == ">":
    	header = line.strip()
        outfile.write(header + "_rc" + "\n")
    else:
    	seq = line.strip()
    	outfile.write(rComp(seq) + adapter_Frc  + "\n")
outfile.close()
F_infile.close()


####################################
# Process the R reads
####################################

### make regular primers with ^seq -> G
outfile = open(prim+ "G_reverse_primers.txt", "w+") # forwards with ^seq
header = ''
seq = ''
for line in R_infile:
    if line[0] == ">":
    	header = line.strip()
        outfile.write(header + "\n")
    else:
    	seq = line.strip()
    	outfile.write("^" + seq + "\n")
outfile.close()
R_infile.close()

### make reverse complement primers with seq$ -> a
R_infile = open(sys.argv[3], "r")
outfile = open(prim + "a_reverse_rc_primers.txt", "w+") # forwards with ^seq
header = ''
seq = ''
for line in R_infile:
    if line[0] == ">":
    	header = line.strip()
        outfile.write(header + "_rc" + "\n")
    else:
    	seq = line.strip()
    	outfile.write(rComp(seq)+ "$" + "\n")
outfile.close()
R_infile.close()

### make reverse complement primers with seq$ -> a
R_infile = open(sys.argv[3], "r")
outfile = open(prim_adapt + "a_Reverse_PrimAdapt_rc.txt", "w+") # forwards with ^seq
header = ''
seq = ''
for line in R_infile:
    if line[0] == ">":
    	header = line.strip()
        outfile.write(header + "_rc" + "\n")
    else:
    	seq = line.strip()
    	outfile.write(rComp(seq) + adapter_Rrc + "\n")
outfile.close()
R_infile.close()
