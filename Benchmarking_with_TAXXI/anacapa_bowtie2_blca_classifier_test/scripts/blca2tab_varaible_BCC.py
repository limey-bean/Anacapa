#!/usr/bin/python2

import sys
import die
import fasta

FileName = sys.argv[1]		# BLCA output
FastaFileName = sys.argv[2]	# Original fasta file with utax2 annotations
Cutoff = int(sys.argv[3])

AccToTax = {}
DomainRank = '?'

def OnSeq(Label, Seq):
	global AccToTax
	global DomainRank

	Acc = fasta.GetAccFromLabel(Label)
	Tax = fasta.GetTaxFromLabel(Label)

	Acc = Acc.upper()
	if Acc.startswith("REF_"):
		Acc = Acc.replace("REF_", "")
	AccToTax[Acc] = Tax

	if DomainRank == '?':
		DomainRank = Tax[0]

fasta.ReadSeqsOnSeq(FastaFileName, OnSeq)

f = open(FileName)
while 1:
	Line = f.readline()
	if len(Line) == 0:
		break

# gi_219856860    superkingdom:Bacteria;100.0;phylum:Actinobacteria;100.0;class:Actinobacteria;100.0;order:Micrococcales;100.0;family:Microbacteriaceae;97.3333333333;genus:Microbacterium;45.8928571429;species:Sp_3690;39.5595238095;
	Fields = Line[:-1].split('\t')
	if len(Fields) != 2:
		die.Die("Got %u fields in line '%s'" % (len(Fields), Line[:-1]))

	Acc = Fields[0]
	Acc = Acc.upper()
	if Acc.startswith("REF_"):
		Acc = Acc.replace("REF_", "")

	QTax = AccToTax[Acc]

	Tax = Fields[1]
	Fields2 = Tax.split(";")
	n = len(Fields2)
	n -= 1
	assert n%2 == 0
	Pred = ""
	for i in range(0, n, 2):
		RankName = Fields2[i]
		BootStr = Fields2[i+1]
		Boot = float(BootStr)
		assert Boot >= 0.0 and Boot <= 100.0
		if Boot < Cutoff:
			break

		Fields3 = RankName.split(':')
		if RankName.find("Not Available") >= 0:
			break

		assert len(Fields3) == 2
		Rank = Fields3[0]
		Name = Fields3[1]
		if Rank == "superkingdom":
			r = DomainRank
		else:
			r = Rank[0]

		if Pred != "":
			Pred += ","

		Pred += r + ":" + Name

	if Pred == "":
		Pred = "*"

	print(Acc + ";tax=" + QTax + ";" + "\t" + Pred)