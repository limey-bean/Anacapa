# Adapted from https://github.com/qunfengdong/BLCA/ by Jesse Gomer
# for the University of California Conservation Consortium's CALeDNA Program

# Example usage:
# python blca_from_bowtie.py -i take_3_local.sam -r CO1_labeled_taxonomy.txt -q CO1_.fasta -b 0.8 -p path/to/muscle
from __future__ import print_function, division

import sys
import os

if sys.version_info < (3, 0):
    from StringIO import StringIO
else:
    from io import StringIO

try:
    from Bio import AlignIO, SeqIO
except ImportError:
    sys.stderr.write("Error! BioPython is not detected!\n")
    sys.exit(1)

import random
import subprocess
import re
from collections import namedtuple, defaultdict
import argparse

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
        self.identity_ratio_old = total_match_count / total_count
        self.identity_ratio = total_match_count / float(len(self.seq))

        self.soft_clipped_ratio = self.cigar_total_soft_clipping() / total_count

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
        beginning_s, ending_s = self.get_soft_clipping()
        return max(beginning_s, ending_s)

    def cigar_total_soft_clipping(self):
        beginning_s, ending_s = self.get_soft_clipping()
        return beginning_s + ending_s

    def get_soft_clipping(self):
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

        return (beginning_s, ending_s)


class NotAvailableHandler(object):
    def __init__(self):
        self.count = 0
        self.na_forms = {'na', 'NA', 'Not Available', 'not available', 'Not available'}
        self.encoded_na_form = 'NA;;'

    def encode_if_na(self, taxon):
        if taxon not in self.na_forms:
            return taxon
        self.count += 1
        return 'NA;;{}'.format(self.count)

    def decode_if_na(self, taxon):
        if not taxon.startswith(self.encoded_na_form):
            return taxon

        return 'NA'


parser = argparse.ArgumentParser(description='Bayesian-based LCA taxonomic classification method')
##### Required arguments #####
required = parser.add_argument_group('required arguments')
required.add_argument("-i", "--sam", help="Input SAM file", type=str, required=True)
required.add_argument('-q', '--reference', help="Reference fasta file", type=str, required=True)
required.add_argument("-r", "--tax", help="reference taxonomy file for the Database", type=str, required=True)

##### Taxonomy filtering arguments #####
taxoptions = parser.add_argument_group('taxonomy profiling options [filtering of hits]')
taxoptions.add_argument("-n", "--nper", help="number of times to bootstrap. Default: 100", type=int, default=100)
taxoptions.add_argument("-b", "--iset", help="minimum identity score to include", type=float, default=0.8)
taxoptions.add_argument('-l', '--length', help="minimum length of hit to include relative to query", type=float,
                        default=0.5)
taxoptions.add_argument('-s', '--softclipping', help='maximum soft clipped ratio to include', type=float, default=0.2)
##### Alignment control arguments #####
alignoptions = parser.add_argument_group('alignment control arguments')
alignoptions.add_argument("-m", "--match", default=1.0, help="alignment match score. Default: 1", type=float)
alignoptions.add_argument("-f", "--mismatch", default=-2.5, help="alignment mismatch penalty. Default: -2.5",
                          type=float)
alignoptions.add_argument("-g", "--ngap", default=-2.0, help="alignment gap penalty. Default: -2", type=float)
##### Other arguments #####
optional = parser.add_argument_group('other arguments')
optional.add_argument('-p', '--muscle', help='Path to call muscle default: muscle', default='muscle')

optional.add_argument("-o","--outfile",help="output file name. Default: <fasta>.blca.out",type=str)
optional.add_argument('--continue_mode', help="continue from a previous run by appending to the same output file", action='store_true')
optional.add_argument("--muscle_use_diags", help="pass the diag argument to muscle", action='store_true')
optional.add_argument("--muscle_max_iterations", help="set the max number of iterations for muscle", type=int, default=16)

##### parse arguments #####
args = parser.parse_args()

### bootstrap times ###
nper = args.nper  # number of bootstrap to permute
### Filter hits per query ###
iset = args.iset  # identify threshold
### Alignment options ###
ngap = args.ngap  # gap penalty
match = args.match  # match score
mismatch = args.mismatch  # mismatch penalty
min_length = args.length
sam_file_name = args.sam
outfile_name = args.outfile or (sam_file_name + '.blca.out')
reference_fasta = args.reference
tax = args.tax
muscle_path = args.muscle
max_soft_clipping_allowed = args.softclipping
continue_mode = args.continue_mode
muscle_use_diags = args.muscle_use_diags
muscle_max_iterations = args.muscle_max_iterations


levels = ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]


def check_program(prgname):
    '''Check whether a program has been installed and put in the PATH'''
    path = os.popen("which " + prgname).read().rstrip()
    if len(path) > 0 and os.path.exists(path):
        print(prgname + " is located in your PATH!")
    else:
        print("ERROR: " + prgname + " is NOT in your PATH, please set up " + prgname + "!")
        sys.exit(1)


def get_dic_from_aln(aln):
    '''Read in alignment and convert it into a dictionary'''
    alignment = AlignIO.read(aln, "clustal")
    alndic = {}
    for r in alignment:
        alndic[r.id] = list(r.seq)
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
    for i in range(len(list(alndic.values())[0])):
        idx.append(random.choice(range(len(list(alndic.values())[0]))))

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


def read_tax_acc(taxfile, not_available_handler):
    tx = open(taxfile)
    acctax = {}
    for l in tx:
        lne = l.rstrip().strip(";").split("\t")
        if len(lne) != 2:
            continue
        if (levels[0] + ':') not in l:
            taxons = [not_available_handler.encode_if_na(taxon) for taxon in lne[1].split(';')]
            acctax[lne[0].split('.')[0]] = dict(zip(levels, taxons))
        else:
            pairs = [x.split(":", 1) for x in lne[1].split(";")]
            encoded = [(level, not_available_handler.encode_if_na(taxon)) for level, taxon in pairs]
            acctax[lne[0].split(".")[0]] = dict(encoded)
    tx.close()
    return acctax


################################################################
##
## 	Running Script Start
##
################################################################

## check whether muscle is located in the path
# check_program("muscle")

### read in pre-formatted lineage information ###
na_handler = NotAvailableHandler()
acc2tax = read_tax_acc(tax, na_handler)
print("> 1 > Read in taxonomy information!")

reference_sequences = {}
with open(reference_fasta) as f:
    for r in SeqIO.parse(f, "fasta"):
        reference_sequences[r.id] = str(r.seq)

print("> 2 > Read in reference db")

SequenceInfo = namedtuple('SequenceInfo', ['seq', 'hits'])
### read in input fasta file ###
input_sequences = {}
possible_rejects = set()
with open(sam_file_name) as sam_file:
    for line in sam_file:
        pieces = line.strip().split('\t')
        entry = SamEntry(pieces)

        if entry.identity_ratio_old < iset:
            possible_rejects.add(entry.qname)
        if entry.soft_clipped_ratio > max_soft_clipping_allowed:
            possible_rejects.add(entry.qname)
        elif entry.rname not in reference_sequences:
            possible_rejects.add(entry.qname)
        elif len(reference_sequences[entry.rname]) / float(len(entry.seq)) < min_length:
            possible_rejects.add(entry.qname)
        elif entry.qname not in input_sequences:
            input_sequences[entry.qname] = SequenceInfo(seq=entry.seq, hits=[entry.rname])
        else:
            input_sequences[entry.qname].hits.append(entry.rname)

rejects = possible_rejects.difference(set(input_sequences))


print "> 3 > Read in bowtie2 output!"

already_assigned = set()
if continue_mode:
    for line in open(outfile_name):
        already_assigned.add(line.split('\t')[0])

    outfile = open(outfile_name, 'a')
else:
    outfile = open(outfile_name, 'w')

for seqn, info in input_sequences.items():
    if seqn in already_assigned:
        continue

print("> 3 > Read in bowtie2 output!")

count = 0
outfile = open(outfile_name, 'w')
for seqn, info in input_sequences.items():
    count += 1

    if seqn in acc2tax:
        print("[WARNING] Your sequence " + seqn + " has the same ID as the reference database! Please correct it!")
        print("...Skipping sequence " + seqn + " ......")
        outfile.write(seqn + "\tSkipped\n")
        continue

    ### Get all the hits list belong to the same query ###
    ### Add query fasta sequence to extracted hit fasta ###
    fifsa = []
    for hit in info.hits:
        if hit not in reference_sequences:
            print("Missing reference sequence for " + hit)
            continue
        fifsa.append(">{}\n{}\n".format(hit, reference_sequences[hit]))
    fifsa.append(">" + seqn + "\n" + info.seq)
    fifsa = "\n".join(fifsa)
    # os.system("rm " + seqn + ".dblist")
    ### Run muscle ###
    muscle_call = [muscle_path, '-quiet', '-clw', '-maxiters',  str(muscle_max_iterations)]
    if muscle_use_diags:
        muscle_call.append('-diags')
    proc = subprocess.Popen(muscle_call, stdout=subprocess.PIPE, stdin=subprocess.PIPE)
    outs, errs = proc.communicate(fifsa.encode('utf-8'))

    # print outs
    # print errs
    # print StringIO.StringIO(outs)
    alndic = get_dic_from_aln(StringIO(outs.decode('utf-8')))
    # os.system("rm " + seqn + ".hits.fsa")
    # os.system("rm " + seqn + ".muscle")
    #    	print "Processing:",k1
    ### get gap position and truncate the alignment###
    start, end = get_gap_pos(seqn, alndic)
    trunc_alndic = cut_gap(alndic, start, end)
    orgscore = pairwise_score(trunc_alndic, seqn, match, mismatch, ngap)
    ### start bootstrap ###
    perdict = {}  # record alignmet score for each iteration
    pervote = {}  # record vote after nper bootstrap

    for j in range(nper):
        random_scores = random_aln_score(trunc_alndic, seqn, match, mismatch, ngap)
        perdict[j] = random_scores
        max_score = max(random_scores.values())
        hits_with_max_score = [k3 for k3, v3 in random_scores.items() if v3 == max_score]
        vote_share = 1.0 / len(hits_with_max_score)
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
            print("Missing taxonomy info for ", short_hit_name)
            continue
        hit_taxonomy = acc2tax[short_hit_name]
        for level in levels:
            # deal with missing values in the taxonomy
            if level not in hit_taxonomy:
                hit_taxonomy[level] = na_handler.encode_if_na("NA")

            if hit in pervote:
                votes_by_level[level][hit_taxonomy[level]] += pervote[hit]
            else:
                votes_by_level[level][hit_taxonomy[level]] += 0

    outfile.write(seqn + "\t")
    for level in levels:
        levels_votes = votes_by_level[level]
        outfile.write(level + ":" + na_handler.decode_if_na(max(levels_votes, key=levels_votes.get)) + ";")
    outfile.write("\t")
    for level in levels:
        levels_votes = votes_by_level[level]
        outfile.write(level + ":" + str(max(levels_votes.values())) + ";")
    outfile.write("\t" + ";".join(info.hits))

    outfile.write("\n")

for seqn in rejects:
    outfile.write(seqn + "\tUnclassified\n")

outfile.close()
