#!/usr/bin/python2

import sys
import sortdict
import utax2
import fasta
import die

#python2 ${maindir}scripts/taxbench.py $preddir/pred_BLCA.out.100 blca o ${maindir} ${name} ${namecountsdir}

Set = sys.argv[1]
Algo = sys.argv[2]
Rank = sys.argv[3]
name = sys.argv[4]
outdir = sys.argv[5]

REPORT = (len(sys.argv) > 6 and sys.argv[6] == "-report")
MIN = 10

Fields = Set.split(".")
n = len(Fields)
LastField = Fields[n-1]
sPctId = LastField
PctId = int(sPctId)
assert PctId in [ 100, 99, 97, 95, 90, 80 ]


FileName_NameCounts = outdir + "/" + name + "_namecount.txt"

print  FileName_NameCounts

FileName_Pred = Set


def ReadNameCounts(FileName):
	NameToCount = {}
	f = open(FileName)
	while 1:
		Line = f.readline()
		if len(Line) == 0:
			return NameToCount
		Fields = Line[:-1].split('\t')
		assert len(Fields) == 2
		Count = int(Fields[0])
		Name = Fields[1]
		assert Name[1] == ':'
		Rank = Name[0]
		NameToCount[Name] = Count		


NameToCount = ReadNameCounts(FileName_NameCounts)

def GetPct(x, y):
	if x > y:
		die.Die("GetPct(%u, %u)" % (x, y))
	if y == 0:
		return 0.0
	return (100.0*x)/y

def DoPredFile(FileName, NameToCount):
	f = open(FileName)
	TP = 0
	TN = 0
	FN = 0
	OC = 0
	MC = 0
	NC = 0
	Known = 0
	Novel = 0

	while 1:
		Line = f.readline()
		if len(Line) == 0:
			return Known, Novel, TP, TN, FN, OC, MC, NC
		
		Fields = Line[:-1].split('\t')
		assert len(Fields) == 2

		QueryLabel = Fields[0]
		Pred = Fields[1]
		if Pred.endswith(';'):
			Pred = Pred[:-1]

		QueryName = utax2.GetNameFromLabel(QueryLabel, Rank)
		if QueryName == "":
			die.Die("Missing name for rank %c: in query label >%s" % (Rank, QueryLabel))

		if Pred.find("tax=") >= 0:
			Pred = Pred.split("tax=")[1]

		if Pred == "*":
			PredName = ""
		else:
			PredName = utax2.GetNameFromTaxStr(Pred, Rank)
		if PredName != "":
			NC += 1

		Count = sortdict.GetCount(NameToCount, QueryName)
		IsKnown = (Count > 0)
		if IsKnown:
			Known += 1
		else:
			Novel += 1

		if PredName == QueryName and not IsKnown:
			die.Die("QueryName=%s, PredName=%s >%s" % (QueryName, PredName, QueryLabel))

		if PredName == QueryName:
			XX = "TP"
			TP += 1
		elif PredName == "":
			if Count == 0:
				XX = "TN"
				TN += 1
			else:
				XX = "FN"
				FN += 1
		else:
			if Count == 0:
				XX = "OC"
				OC += 1
			else:
				XX = "MC"
				MC += 1

		if REPORT:
			Acc = fasta.GetAccFromLabel(QueryLabel)
			if IsKnown:
				k = "known"
			else:
				k = "novel"

			PredNameStr = "-"
			if PredName != "":
				PredNameStr = PredName
			s = Acc
			s += "\t" + XX
			s += "\t" + k
			s += "\t" + QueryName
			s += "\t" + PredNameStr
			s += "\t" + str(Count)
			print s

Known, Novel, TP, TN, UC, OC, MC, NC = DoPredFile(FileName_Pred, NameToCount)

Total = Known + Novel

Resolved = Known + OC

CR = GetPct(NC, Total)

TPR = GetPct(TP, Known)
UCR = GetPct(UC, Known)
MCR = GetPct(MC, Known)
OCR = GetPct(OC, Novel)

Acc = GetPct(TP, Known + OC)

def ParseSet(Set):
	if Set == "rdp_its":
		DB = "rdp_its"
		Seg = "fl"
	elif Set.find("_v4") > 0:
		DB = Set.split("_v4")[0]
		Seg ="v4"
	elif Set.find("_v35") > 0:
		Seg ="v35"
		DB = Set.split("_v35")[0]
	else:
		DB = Set
		Seg ="fl"
	return DB, Seg

DB, Seg = ParseSet(Set)

def Out(Metric, Value):
	print(DB + "\t" + Seg + "\t" + str(PctId) + "\t" + Metric + "\t" + Algo + "\t" + Rank + "\t" + str(Value))

sCR = "%.1f" % CR

sTPR = "%.1f" % TPR
sUCR = "%.1f" % UCR
sMCR = "%.1f" % MCR
sOCR = "%.1f" % OCR
sAcc = "%.1f" % Acc

if Novel < MIN:
	sOCR = "."

if Known < MIN:
	sTPR = "."
	sUCR = "."
	sMCR = "."

if Known + OC < MIN:
	sAcc = "."

Out("TPR", sTPR)
Out("UCR", sUCR)
Out("MCR", sMCR)
Out("OCR", sOCR)
Out("Acc", sAcc)