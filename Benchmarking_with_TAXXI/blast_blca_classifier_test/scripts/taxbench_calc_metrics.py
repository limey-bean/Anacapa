#!/usr/bin/python2

import sys
import sortdict
import utax2
import fasta
import die

Set = sys.argv[1]
Algo = sys.argv[2]
Dist = sys.argv[3]
Rank = sys.argv[4]

REPORT = (len(sys.argv) > 5 and sys.argv[5] == "-report")
MIN = 10

def DistToPctId(Dist):
	iPctId = 100 - int(Dist)
	return "%s" % iPctId

PctId = DistToPctId(Dist)

FileName_NameCountsA = "/e/res/taxbench/namecounts/" + Set + "." + Dist + ".a"
FileName_NameCountsB = "/e/res/taxbench/namecounts/" + Set + "." + Dist + ".b"

FileName_PredAB = "/e/res/taxbench/pred/" + Algo + "/" + Set + "." + Dist + ".ab"
if Dist == "00":
	FileName_PredBA = FileName_PredAB
else:
	FileName_PredBA = "/e/res/taxbench/pred/" + Algo + "/" + Set + "." + Dist + ".ba"

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

NameToCountA = ReadNameCounts(FileName_NameCountsA)
NameToCountB = ReadNameCounts(FileName_NameCountsB)

def GetPct(x, y):
	if x > y:
		die.Die("GetPct(%u, %u)" % (x, y))
	if y == 0:
		return 0.0
	return (100.0*x)/y

def DoPredFile(FileName, NameToCountA, NameToCountB):
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
		assert QueryName != ""

		if Pred.find("tax=") >= 0:
			Pred = Pred.split("tax=")[1]

		if Pred == "*":
			PredName = ""
		else:
			PredName = utax2.GetNameFromTaxStr(Pred, Rank)
		if PredName != "":
			NC += 1

		Count = sortdict.GetCount(NameToCountB, QueryName)
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

Known_AB, Novel_AB, TP_AB, TN_AB, FN_AB, OC_AB, MC_AB, NC_AB = DoPredFile(FileName_PredAB, NameToCountA, NameToCountB)
Known_BA, Novel_BA, TP_BA, TN_BA, FN_BA, OC_BA, MC_BA, NC_BA = DoPredFile(FileName_PredBA, NameToCountB, NameToCountA)

Known = Known_AB + Known_BA
Novel = Novel_AB + Novel_BA
Total = Known + Novel

Correct = TP_AB + TP_BA
Overclassified = OC_AB + OC_BA
Resolved = Known + Overclassified

NC = NC_AB + NC_BA

CR = GetPct(NC, Total)

TPR_AB = GetPct(TP_AB, Known_AB)
TPR_BA = GetPct(TP_BA, Known_BA)

FNR_AB = GetPct(FN_AB, Known_AB)
FNR_BA = GetPct(FN_BA, Known_BA)

MCR_AB = GetPct(MC_AB, Known_AB)
MCR_BA = GetPct(MC_BA, Known_BA)

TNR_AB = GetPct(TN_AB, Novel_AB)
TNR_BA = GetPct(TN_BA, Novel_BA)

OCR_AB = GetPct(OC_AB, Novel_AB)
OCR_BA = GetPct(OC_BA, Novel_BA)

if 0:
	print "TPR_AB  %5.1f%%  %7u / %7u" % (TPR_AB, TP_AB, Known_AB)
	print "TPR_BA  %5.1f%%  %7u / %7u" % (TPR_BA, TP_BA, Known_BA)

	print "FNR_AB  %5.1f%%  %7u / %7u" % (FNR_AB, FN_AB, Known_AB)
	print "FNR_BA  %5.1f%%  %7u / %7u" % (FNR_BA, FN_BA, Known_BA)

	print "MCR_AB  %5.1f%%  %7u / %7u" % (MCR_AB, MC_AB, Known_AB)
	print "MCR_BA  %5.1f%%  %7u / %7u" % (MCR_BA, MC_BA, Known_BA)

	print "TNR_AB  %5.1f%%  %7u / %7u" % (TNR_AB, TN_AB, Novel_AB)
	print "TNR_BA  %5.1f%%  %7u / %7u" % (TNR_BA, TN_BA, Novel_BA)

	print "OCR_AB  %5.1f%%  %7u / %7u" % (OCR_AB, OC_AB, Novel_AB)
	print "OCR_BA  %5.1f%%  %7u / %7u" % (OCR_BA, OC_BA, Novel_BA)

TPR = (TPR_AB + TPR_BA)/2.0
FNR = (FNR_AB + FNR_BA)/2.0
MCR = (MCR_AB + MCR_BA)/2.0
TNR = (TNR_AB + TNR_BA)/2.0
OCR = (OCR_AB + OCR_BA)/2.0

Acc = GetPct(Correct, Resolved)

if 0:
	s = Set
	s += "\tAlgo"
	s += "\tPctId"
	s += "\tKnown"
	s += "\tNovel"
	s += "\tTPR"
	s += "\tFNR"
	s += "\tMCR"
	s += "\tTNR"
	s += "\tOCR"
	s += "\tCR"
	s += "\tAcc"
	print s

	s = Set
	s += "\t" + Algo
	s += "\t%.0f" % PctId
	s += "\t%u" % Known
	s += "\t%u" % Novel

	if Known >= MIN:
		s += "\t%.1f" % TPR
		s += "\t%.1f" % FNR
		s += "\t%.1f" % MCR
	else:
		s += "\t."
		s += "\t."
		s += "\t."

	if Novel >= MIN:	
		s += "\t%.1f" % TNR
		s += "\t%.1f" % OCR
	else:
		s += "\t."
		s += "\t."

	s += "\t%.1f" % CR
	s += "\t%.1f" % Acc
	print s

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
	print(DB + "\t" + Seg + "\t" + PctId + "\t" + Metric + "\t" + Algo + "\t" + Rank + "\t" + str(Value))

sCR = "%.1f" % CR

sTPR = "%.1f" % TPR
sFNR = "%.1f" % FNR
sTNR = "%.1f" % TNR
sMCR = "%.1f" % MCR
sOCR = "%.1f" % OCR
sAcc = "%.1f" % Acc

if Novel < MIN:
	sTNR = "."
	sOCR = "."
if Known < MIN:
	sTPR = "."
	sFNR = "."
	sMCR = "."

Out("CR", sCR)
Out("TPR", sTPR)
Out("FNR", sFNR)
Out("MCR", sMCR)
Out("TNR", sTNR)
Out("OCR", sOCR)
Out("Acc", sAcc)
