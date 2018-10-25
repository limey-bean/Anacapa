#!/usr/bin/python2

import sys
import fasta
import utax2

InputFileName = sys.argv[1]
FastaOutFileName = sys.argv[2]
TaxOutFileName = sys.argv[3]

fFa = open(FastaOutFileName, "w")
fTax = open(TaxOutFileName, "w")

SpCount = 0

def AppendRank(Tax, TaxStr, Rank):
	global SpCount

	if Rank == 'k' or Rank == 'd':
		RankName = "superkingdom"
	else:
		RankName = utax2.CharToLevelName(Rank)
	Name = utax2.GetNameFromTaxStr(TaxStr, Rank)
	if Name == "":
		if Rank == "s":
			SpCount += 1
			return Tax + "species:Sp_" + str(SpCount) + ";"
		assert False
	return Tax + RankName + ":" + Name[2:] + ";"

def OnSeq(Label, Seq):
	Acc = fasta.GetAccFromLabel(Label)
	TaxStr = utax2.GetTaxFromLabel(Label)
	k = 'k'
	if TaxStr.find("d:") >= 0:
		k = 'd'

# NR_117221.1     species:Mycobacterium arosiense;genus:Mycobacterium;family:Mycobacteriaceae;order:Corynebacteriales;class:Actinobacteria;phylum:Actinobacteria;superkingdom:Bacteria;

	Tax = ""
	if TaxStr.find(",s:") > 0:
		Tax = AppendRank(Tax, TaxStr, 's')
	Tax = AppendRank(Tax, TaxStr, 'g')
	Tax = AppendRank(Tax, TaxStr, 'f')
	Tax = AppendRank(Tax, TaxStr, 'o')
	Tax = AppendRank(Tax, TaxStr, 'c')
	Tax = AppendRank(Tax, TaxStr, 'p')
	Tax = AppendRank(Tax, TaxStr, k)

	Acc = Acc.upper()	
	NewLabel = "REF_" + Acc
	fasta.WriteSeq(fFa, Seq, NewLabel)

	print >> fTax, "%s\t%s" % (NewLabel, Tax)

fasta.ReadSeqsOnSeq(InputFileName, OnSeq)

