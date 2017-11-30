# Adapted from https://github.com/qunfengdong/BLCA/ by Jesse Gomer
# for the University of California Conservation Consortium's CALeDNA Program

# Example usage python blca_from_bowtie.py  -i take_3_local.sam -r CO1_labeled_taxonomy.txt -q CO1_.fasta  -b 0.8
import sys
import os
import math

try:
    from Bio import AlignIO, SeqIO
except ImportError:
    sys.stderr.write("Error! BioPython is not detected!\n")
    sys.exit(1)

import random
import subprocess
import getopt
import re
from collections import namedtuple, defaultdict

# import argparse

'''
BLCA Core annotation tool
'''


class SamEntry(object):
    def __init__(self, raw_row):
        self.qname = raw_row[0]
        self.flag = raw_row[1]
        self.rname = raw_row[2]
        self.pos = raw_row[3]
        self.mapq = raw_row[4]
        self.cigar = raw_row[5]
        self.rnext = raw_row[6]
        self.pnext = raw_row[7]
        self.tlen = raw_row[8]
        self.seq = raw_row[9]
        self.qual = raw_row[10]
        self.alignment_scores = [int(score.split(':')[-1]) for score in raw_row[11:13]]
        # sometimes it is a column early  because there is no second best match
        if raw_row[17][:2] == 'MD':
            self.md_z_flags = raw_row[17].split(':')[-1]
        else:
            self.md_z_flags = raw_row[18].split(':')[-1]

        total_match_count, total_count = self.calulate_match_count(self.md_z_flags)
        self.identity_ratio = total_match_count / total_count
        self.total_match = total_match_count

    def calulate_match_count(self, md_z_flags):
        tokenizer = re.compile(r'(\d+)|(\^[A-Z])|([A-Z])')
        total_match_count = 0.0
        total_mismatch_count = 0.0
        for item in tokenizer.finditer(md_z_flags):
            match_count, deletion, mismatch = item.groups()
            if match_count:
                total_match_count += int(match_count)
            if deletion:
                total_mismatch_count += 1
            if mismatch:
                total_mismatch_count += 1

        return total_match_count, (total_match_count + total_mismatch_count)

    # I am not completely sure if this is the logic that you want to use for checking unmapped at the ends
    # This function returns the max S at either end of the CIGAR score, 0  if there is no S at both ends
    def cigar_max_s(self):
        tokenizer = re.compile(r'(?:\d+)|(?:[A-Z=])')
        cigar_elements = tokenizer.findall(self.cigar)

        if cigar_elements[1] == 'S':
            beginning_s = int(cigar_elements[0])
        else:
            beginning_s = 0

        if cigar_elements[-1] == 'S':
            ending_s = int(cigar_elements[-2])
        else:
            ending_s = 0

        return max(beginning_s, ending_s)


def usage():
    print "\n<< Bayesian-based LCA taxonomic classification method >>\n\n   Please make sure the following softwares are in your PATH:\n\t1.muscle (http://www.drive5.com/muscle/downloads.htm), muscle should be the program's name.\n\t2.ncbi-blast suite (ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)\n\t3.Biopython should be installed locally.\n"
    print 'Usage: python ' + sys.argv[0] + ' -i <sam file> [option]\n'
    print " \nArguments:\n - Required:"
    print "\t-i\t\tInput fasta file.\n - Taxonomy Profiling Options [filtering of hits]:"
    print "\t-n\t\tNumber of times to bootstrap. Default: 100"
    print "\t-j\t\tMaximum number of subjects to include for each query reads. Default: 50"
    print "\t-d\t\tProportion of hits to include from top hit. Default: 0.1 [0-1]"
    print "\t-e\t\tMinimum evalue to include for blastn. Default: 0.1"
    print "\t-a\t\tMinimum bitscore to include for blastn hits. Default: 100"
    print "\t-c\t\tMinimum coverage to include. Default: 0.85 [0-1]"
    print "\t-b\t\tMinimum identity score to include. Default: 0.90 [0-1.0]"
    print "\t-r\t\tReference Taxonomy file for the Database. Default: db/16SMicrobial.ACC.taxonomy"
    print "\t-q\t\tRefernece blast database. Default: db/16SMicrobial"
    print "\t-o\t\tOutput file name. Default: <fasta>.blca.out\n - Alignment Options:"
    print "\t-m\t\tAlignment match score. Default: 1"
    print "\t-f\t\tAlignment mismatch penalty. Default: -2.5"
    print "\t-g\t\tAlignment gap penalty. Default: -2\n - Other:"
    print "\t-t\t\tExtra number of nucleotides to include at the beginning and end of the hits. Default: 10"
    print "\t-h\t\tShow program usage and quit"


####set up default parameters #######
### bootstrap times ###
nper = 100  # number of bootstrap to permute
### Filter hits per query ###
iset = float(0.9)  # identify threshold
eset = float(0.1)  # evalue threshold
gap = 10  # number of nucleotides to include at the beginning and end of the hits
topper = float(0.1)  # top percentage to include hits
cvrset = float(0.80)  # coverage to include
bset = float(100)  # minimum bitscore to include
nsub = 50  # maximum number of subjects to include
### Alignment options ###
ngap = -2  # gap penalty
match = 1  # match score
mismatch = -2.5  # mismatch penalty
levels = ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]

opts, args = getopt.getopt(sys.argv[1:], "a:b:c:d:e:f:g:i:j:lm:n:o:r:q:t:s:h",
                           ['Minimum bitscore', 'Minimum Identity', 'Minimum Coverage', 'Top Proportion',
                            'Minimum evalue', 'Alignment Mismatch Penalty', 'Alignment Gap Penalty', 'Input File',
                            'Maximum Subjects', 'Long Output', 'Alignment Match Score', 'Number of Bootstrap to Run',
                            'Output File', 'Taxonomy File of Database', 'Reference FASTA file',
                            'Number of nt Length of Sequence', 'help'])
for o, a in opts:
    if o == "-i":
        sam_file_name = a
        outfile_name = sam_file_name + '.blca.out'
    elif o == "-n":
        nper = float(a)
    elif o == "-j":
        nsub = int(a)
    elif o == "-t":
        gap = float(a)
    elif o == "-c":
        cvrset = float(a)
    elif o == "-l":
        longout = True
    elif o == "-d":
        topper = float(a)
    elif o == "-e":
        eset = float(a)
    elif o == "-f":
        mismatch = float(a)
    elif o == "-g":
        ngap = float(a)
    elif o == "-q":
        reference_fasta = a
    elif o == "-m":
        match = float(a)
    elif o == "-b":
        iset = float(a)
    elif o == "-o":
        outfile_name = a
    elif o == "-r":
        tax = a
    elif o == "-a":
        bset = float(a)
    elif o in ('-h', '--help'):
        print usage()
        sys.exit(1)
    else:
        assert False, 'unhandle option'


def check_taxdb():
    ''' Check whether BLASTDB environmental variable has been setup'''
    if 'BLASTDB' not in os.environ.keys():
        print "ERROR: taxdb location has not been set up as the environmental variable BLASTDB. Please set up as \n\texport BLASTDB=/location/of/taxdb.bti/and/taxdb.btd/"
        sys.exit(1)


def check_program(prgname):
    '''Check whether a program has been installed and put in the PATH'''
    path = os.popen("which " + prgname).read().rstrip()
    if len(path) > 0 and os.path.exists(path):
        print prgname + " is located in your PATH!"
    else:
        print "ERROR: " + prgname + " is NOT in your PATH, please set up " + prgname + "!"
        sys.exit(1)


def get_dic_from_aln(aln):
    '''Read in alignment and convert it into a dictionary'''
    alignment = AlignIO.read(aln, "clustal")
    alndic = {}
    for r in alignment:
        alndic[r.id] = r.seq
    return alndic


def pairwise_score(alndic, query, match, mismatch, ngap):
    '''Calculate pairwise alignment score given a query'''
    nt = ["A", "C", "T", "G", "g", "a", "c", "t"]
    hitscore = {}
    for k, v in alndic.items():
        if k != query:
            hitscore[k] = 0
            for i in range(len(v)):
                if (alndic[query][i] in nt) and (v[i] in nt) and (alndic[query][i] == v[i]):
                    hitscore[k] += float(match)
                elif (alndic[query][i] not in nt) and (v[i] not in nt) and (alndic[query][i] == v[i]):
                    hitscore[k] += float(0)
                elif ((alndic[query][i] not in nt) or (v[i] not in nt)) and (alndic[query][i] != v[i]):
                    hitscore[k] += float(mismatch)
                elif (alndic[query][i] in nt) and (v[i] in nt) and (alndic[query][i] != v[i]):
                    hitscore[k] += float(ngap)
    total = float(sum(hitscore.values()))
    if total <= 0:
        total = 1
    for k, v in hitscore.items():
        hitscore[k] = v / total
    return hitscore


def random_aln_score(alndic, query, match, mismatch, ngap):
    '''Randomize the alignment, and calculate the score'''
    nt = ["A", "C", "T", "G", "g", "a", "c", "t"]
    idx = []
    for i in range(len(alndic.values()[0])):
        idx.append(random.choice(range(len(alndic.values()[0]))))
    hitscore = {}
    for k, v in alndic.items():
        if k != query:
            hitscore[k] = 0
            for i in idx:
                if (alndic[query][i] in nt) and (v[i] in nt) and (alndic[query][i] == v[i]):
                    hitscore[k] += float(match)
                elif (alndic[query][i] not in nt) and (v[i] not in nt) and (alndic[query][i] == v[i]):
                    hitscore[k] += float(0)
                elif ((alndic[query][i] not in nt) or (v[i] not in nt)) and (alndic[query][i] != v[i]):
                    hitscore[k] += float(mismatch)
                elif (alndic[query][i] in nt) and (v[i] in nt) and (alndic[query][i] != v[i]):
                    hitscore[k] += float(ngap)
    return hitscore


def get_gap_pos(query, alndic):
    '''Get the gap position in the alignment'''
    for i in range(len(alndic[query])):
        if alndic[query][i] != "-":
            start = i
            break
    for i in range(len(alndic[query]) - 1, 0, -1):
        if alndic[query][i] != "-":
            end = i
            break
    return start, end


def cut_gap(alndic, start, end):
    '''Given a start and end gap position, truncate the alignmnet'''
    trunc_alndic = {}
    for k_truc, v_truc in alndic.items():
        trunc_alndic[k_truc] = v_truc[start:end]
    return trunc_alndic


def read_tax_acc(taxfile):
    tx = open(taxfile)
    acctax = {}
    print "> 3 > Read in taxonomy information!"
    for l in tx:
        lne = l.rstrip().strip(";").split("\t")
        if len(lne) != 2:
            continue
        if (levels[0] + ':') not in l:
            acctax[lne[0].split('.')[0]] = dict(zip(levels, lne[1].split(';')))
        else:
            acctax[lne[0].split(".")[0]] = dict(x.split(":", 1) for x in lne[1].split(";"))
    tx.close()
    return acctax


################################################################
##
## 	Running Script Start
##
################################################################

## check taxdb
# check_taxdb()

## check whether blastdbcmd is located in the path
#check_program("blastdbcmd")

## check whether muscle is located in the path
check_program("muscle")

### read in pre-formatted lineage information ###
acc2tax = read_tax_acc(tax)
SequenceInfo = namedtuple('SequenceInfo', ['seq', 'hits'])
### read in input fasta file ###
input_sequences = {}
possible_rejects = set()
with open(sam_file_name) as sam_file:
    for line in sam_file:
        pieces = line.strip().split('\t')
        entry = SamEntry(pieces)
        # entry does not match filter
        if entry.identity_ratio < iset:
            possible_rejects.add(entry.qname)
        elif entry.qname in input_sequences:
            input_sequences[entry.qname].hits.append(entry.rname)
        else:
            input_sequences[entry.qname] = SequenceInfo(seq=entry.seq, hits=[entry.rname])

rejects = possible_rejects.difference(set(input_sequences))
print "> 1 > Read in bowtie2 output!"

reference_sequences = {}
with open(reference_fasta) as f:
    for r in SeqIO.parse(f, "fasta"):
        reference_sequences[r.id] = str(r.seq)

outfile = open(outfile_name, 'w')
for seqn, info in input_sequences.items():
    if seqn in acc2tax:
        print "[WARNING] Your sequence " + seqn + " has the same ID as the reference database! Please correct it!"
        print "...Skipping sequence " + seqn + " ......"
        outfile.write(seqn + "\tSkipped\n")
        continue

    ### Get all the hits list belong to the same query ###
    ### Add query fasta sequence to extracted hit fasta ###
    fifsa = open(seqn + ".hits.fsa", 'w')
    for hit in info.hits:
        if hit not in reference_sequences:
            print "Missing reference sequence for " + hit
            continue
        fifsa.write(">{}\n{}\n".format(hit, reference_sequences[hit]))
    fifsa.write(">" + seqn + "\n" + info.seq)
    fifsa.close()
    # os.system("rm " + seqn + ".dblist")
    ### Run muscle ###
    os.system("muscle -quiet -clw -in " + seqn + ".hits.fsa -out " + seqn + ".muscle")
    alndic = get_dic_from_aln(seqn + ".muscle")
    os.system("rm " + seqn + ".hits.fsa")
    os.system("rm " + seqn + ".muscle")
    #    	print "Processing:",k1
    ### get gap position and truncate the alignment###
    start, end = get_gap_pos(seqn, alndic)
    trunc_alndic = cut_gap(alndic, start, end)
    orgscore = pairwise_score(trunc_alndic, seqn, match, mismatch, ngap)
    ### start bootstrap ###
    perdict = {}  # record alignmet score for each iteration
    pervote = {}  # record vote after nper bootstrap
    ### If any equal score, average the vote ###
    for j in range(nper):
        random_scores = random_aln_score(trunc_alndic, seqn, match, mismatch, ngap)
        perdict[j] = random_scores
        max_score = max(random_scores.values())
        hits_with_max_score = [k3 for k3, v3 in random_scores.items() if v3 == max_score]
        vote_share = 1.0/len(hits_with_max_score)
        for hit in hits_with_max_score:
            if hit in pervote:
                pervote[hit] += vote_share
            else:
                pervote[hit] = vote_share

    ### normalize vote by total votes ###
    ttlvote = sum(pervote.values())
    for k4, v4 in pervote.items():
        pervote[k4] = v4 / ttlvote * 100
    ###

    votes_by_level = {}
    for level in levels:
        votes_by_level[level] = defaultdict(int)

    for hit in orgscore.keys():
        short_hit_name = hit.split(".")[0]
        if short_hit_name not in acc2tax:
            print "Missing taxonomy info for ", short_hit_name
            continue
        hit_taxonomy = acc2tax[short_hit_name]
        for level in levels:
            # deal with missing values in the taxonomy
            if level not in hit_taxonomy:
                hit_taxonomy[level] = "Not available"

            if hit in pervote:
                votes_by_level[level][hit_taxonomy[level]] += pervote[hit]
            else:
                votes_by_level[level][hit_taxonomy[level]] += 0

    outfile.write(seqn + "\t")
    for level in levels:
        levels_votes = votes_by_level[level]
        outfile.write(level + ":" + max(levels_votes, key=levels_votes.get) + ";" + str(max(levels_votes.values())) + ";")
    outfile.write("\n")

for seqn in rejects:
    outfile.write(seqn + "\tUnclassified\n")

outfile.close()
print ">> Taxonomy file generated!!"
