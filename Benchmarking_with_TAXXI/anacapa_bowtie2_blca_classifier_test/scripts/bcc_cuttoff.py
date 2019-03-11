#!/usr/bin/python2

import sys
import utax2
import die

FileName = sys.argv[1]
Cutoff = float(sys.argv[2])

f = open(FileName)
WarningDone = False
while 1:
	Line = f.readline()
	if len(Line) == 0:
		break
	Fields = Line[:-1].split('\t')
	assert len(Fields) >= 2

	Label = Fields[0]
	PredStr = Fields[1]

	if PredStr.endswith(",") and not WarningDone:
		die.Warning("Pred ends in comma '%s'" % FileName)
		WarningDone = True

	Preds = utax2.GetPredsFromPredStr(PredStr)
	s = ""
	for Pred in Preds:
		Name, Score = utax2.ParsePred(Pred)
		if Score < Cutoff:
			break
		if s != "":
			s += ","
		s += Name
	print "%s\t%s" % (Label, s) 