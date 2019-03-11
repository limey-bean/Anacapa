#!/usr/bin/python2

import sys
import utax2
import sortdict
import progress

# Distance matrix in mothur format (can generate
# with usearch -calc_distmx).
# Sequence labels must have tax= annotations in
# usearch format.

FileName = sys.argv[1]

LCRToCount = {}
PctIdToCount = {}
PctIdToLCRCount = {}

def Add(LCR, PctId):
	sortdict.IncCount(LCRToCount, LCR)
	sortdict.IncCount(PctIdToCount, PctId)
	sortdict.IncCount2(PctIdToLCRCount, PctId, LCR)

MinPctId = 100
File = open(FileName)
progress.FileInit(File, FileName)
while 1:
	Line = File.readline()
	progress.FileStep()
	if len(Line) == 0:
		progress.FileDone()
		break

	Fields = Line.strip().split('\t')
	assert len(Fields) == 3
	Label1 = Fields[0]
	Label2 = Fields[1]
	Dist = float(Fields[2])
	assert Dist >= 0.0 and Dist <= 1.0
	if Label1 == Label2:
		continue
	LCR = utax2.GetLCRFromLabels(Label1, Label2)
	PctId = int(100.0*(1.0 - Dist))
	if PctId < MinPctId:
		MinPctId = PctId
	Add(LCR, PctId)

LCRsK = LCRToCount.keys()
LCRs = []
for LCR in utax2.LevelChars:
	if LCR in LCRsK:
		LCRs.append(LCR)
LCRCount = len(LCRs)

s = "P(LCR)"
for i in range(0, LCRCount):
	LCR = LCRs[LCRCount-i-1]
	s += "\t%s" % utax2.CharToLevelName(LCR)
print(s)

for PctId in range(100, MinPctId-1, -1):
	s = "%u" % PctId
	N = PctIdToCount[PctId]

	SumFreq = 0.0
	Sumn = 0
	for i in range(0, LCRCount):
		LCR = LCRs[LCRCount-i-1]
		n = sortdict.GetCount2(PctIdToLCRCount, PctId, LCR)
		Freq = 0.0
		if n > 0:
			Freq = float(n)/N
		SumFreq += Freq
		Sumn += n
		s += "\t%.4f" % Freq
	assert Sumn == N
	assert SumFreq > 0.99 and SumFreq < 1.01
	print(s)

print("")
s = "P(CR)"
for i in range(0, LCRCount):
	# LCR = LCRs[i]
	LCR = LCRs[LCRCount-i-1]
	s += "\t%s" % utax2.CharToLevelName(LCR)
print(s)

for PctId in range(100, MinPctId-1, -1):
	s = "%u" % PctId
	N = PctIdToCount[PctId]

	SumFreq = 0.0
	Sumn = 0
	for i in range(0, LCRCount):
		# LCR = LCRs[i]
		LCR = LCRs[LCRCount-i-1]
		n = sortdict.GetCount2(PctIdToLCRCount, PctId, LCR)
		Freq = 0.0
		if n > 0:
			Freq = float(n)/N
		SumFreq += Freq
		Sumn += n
		s += "\t%.4f" % SumFreq
	assert Sumn == N
	assert SumFreq > 0.99 and SumFreq < 1.01
	print(s)
