#!/usr/bin/python2

import sys
import fasta
import utax2

# Input is FASTA in utax2 ref format:
#	>Acc;tax=d:Bacteria,...;

# Output is mothur format
#	FASTA
#		 >Acc
#
#	Tax	is tabbed text, e.g.:
#		AB002527.1      Bacteria;Firmicutes;Bacilli;Lactobacillales;Streptococcaceae;Streptococcus;

FastaFileName = sys.argv[1]
FastaOutFileName = sys.argv[2]
TaxOutFileName = sys.argv[3]

SeqDict = fasta.ReadSeqsDict(FastaFileName)

Labels = SeqDict.keys()
LevelChars = utax2.GetTopLevelChars(Labels)
RequiredLevelCount = len(LevelChars)

fFa = open(FastaOutFileName, "w")
fTax = open(TaxOutFileName, "w")

SeqCount = 0
MissingRankCount = 0
EmptyTaxCount = 0
BadLevelCount = 0

def OnSeq(Label, Seq):
	global MissingRankCount
	global EmptyTaxCount
	global BadLevelCount
	global SeqCount

	SeqCount += 1

	Acc = utax2.GetAccFromLabel(Label)
	Acc = Acc.replace("_", "")
	TaxStr = utax2.GetTaxFromLabel(Label)

	Fields = TaxStr.split(',')
	if len(Fields) != RequiredLevelCount:
		BadLevelCount += 1
		return
	MotTaxStr = ""
	LastRankIndex = -1
	LastLevelChar = ""
	for Field in Fields:
		assert Field[1] == ':'
		LevelChar = Field[0]
		RankIndex = utax2.LevelCharToRankIndex(LevelChar)
		if LastRankIndex != -1 and LastRankIndex > 1 and not RankIndex == LastRankIndex + 1:
			MissingRankCount += 1
#			print >> sys.stderr, "Missing rank: %c,%d %c,%d %s" % (LastLevelChar, LastRankIndex, LevelChar, RankIndex, TaxStr)
			break
		LastRankIndex = RankIndex
		LastLevelChar = LevelChar
		if MotTaxStr != "":
			MotTaxStr += ';'
		MotTaxStr += LevelChar + ":" + Field[2:]

	if MotTaxStr == "":
		EmptyTaxCount += 1
		return

	MotTaxStr += ";"

	fasta.WriteSeq(fFa, Seq, Acc)

	print >> fTax, Acc + "\t" + MotTaxStr

for Label in Labels:
	Seq = SeqDict[Label]
	OnSeq(Label, Seq)

print >> sys.stderr, "Seqs %u, missing rank %u, empty %u, BadLevelCount %u" % (SeqCount, MissingRankCount, EmptyTaxCount, BadLevelCount)
