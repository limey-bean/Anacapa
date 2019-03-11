#!/usr/bin/python2

import sys
import fasta
import utax2
import die

# Input format (utax2 format):
# >EF115542|S001020530;tax=d:Bacteria,p:Cyanobacteria/Chloroplast,c:Chloroplast,f:Chloroplast,g:Streptophyta;

# Two output files:
# FASTA: >EF115542|S001020530
# Taxonomy: EF115542|S001020530	p__Bacteroidetes; c__Flavobacteriia; o__Flavobacteriales; f__Flavobacteriaceae; g__Flavobacterium; s__

FastaFileName = sys.argv[1]
FastaOutFileName = sys.argv[2]
TaxOutFileName = sys.argv[3]
KeepTaxAnnot = False
if len(sys.argv) > 4:
	assert sys.argv[4] == "-k"
	KeepTaxAnnot = True

SeqDict = fasta.ReadSeqsDict(FastaFileName)

Labels = SeqDict.keys()
LevelChars = utax2.GetTopLevelChars(Labels)
RequiredLevelCount = len(LevelChars)

fFa = open(FastaOutFileName, "w")
fTax = open(TaxOutFileName, "w")

def OnSeq(Label, Seq):
	Acc = Label.split(";")[0]
	Tax = utax2.GetTaxFromLabel(Label)
	if KeepTaxAnnot:
		Acc = Label

	Chars, Names = utax2.TaxToVecs(Tax)
	N = len(Chars)
	assert len(Names) == N

	N = len(Chars)
	OutLine = Acc + "\t"
	for i in range(0, N):
		Char = Chars[i]
		Name = Names[i]
		OutLine += Char + "__" + Name
		if i+1 < N:
			OutLine += "; "
			
	print >> fTax, OutLine

	if KeepTaxAnnot:
		fasta.WriteSeq(fFa, Seq, Label)
	else:
		fasta.WriteSeq(fFa, Seq, Acc)

for Label in Labels:
	Seq = SeqDict[Label]
	OnSeq(Label, Seq)
