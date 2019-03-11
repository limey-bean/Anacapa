#!/usr/bin/python2

import sys

LCRFileName = sys.argv[1]
IdDistFileName = sys.argv[2]	# PctId	Count
FieldNr = int(sys.argv[3])		# Column nr. with Count in id dist tabbed file
Total = int(sys.argv[4])		# Total nr seqs (because some may be < search threshold)
REPORT = False

def ReadLCRProbs(FileName):
	PctIdToRankToPLCR = {}
	PctIdToRankToPCR = {}
	f = open(FileName)
	Line = f.readline()
	HdrFields = Line.strip().split('\t')
	assert HdrFields[0] == "P(LCR)"
	Ranks = HdrFields[1:]
	N = len(Ranks)
	for PctId in range(100, 0, -1):
		Line = f.readline().strip()
		if len(Line) == 0:
			f.close()
			return Ranks, PctIdToRankToPLCR, PctIdToRankToPCR

		PctIdToRankToPLCR[PctId] = {}
		PctIdToRankToPCR[PctId] = {}
		Fields = Line.split('\t')
		len(Fields) == N + 1
		assert int(Fields[0]) == PctId
		SumProb = 0.0
		for i in range(0, N):
			Prob = float(Fields[i+1])
			SumProb += Prob
			assert Prob >= 0.0 and Prob <= 1.0
			Rank = Ranks[i]
			PctIdToRankToPLCR[PctId][Rank] = Prob
			PctIdToRankToPCR[PctId][Rank] = SumProb

def ReadIdDist(FileName):
	SumCount = 0
	PctIdToCount = {}
	f = open(FileName)
	for PctId in range(100, 0, -1):
		Line = f.readline().strip()
		if len(Line) == 0 or Line.startswith("0\t"):
			f.close()
			return PctIdToCount

		Fields = Line.split('\t')
		assert int(Fields[0]) == PctId
		Count = int(Fields[FieldNr-1])
		SumCount += Count
		PctIdToCount[PctId] = Count

Ranks, PctIdToRankToPLCR, PctIdToRankToPCR = ReadLCRProbs(LCRFileName)
PctIdToCount = ReadIdDist(IdDistFileName)

SumCount = 0
for PctId in PctIdToCount.keys():
	Count = PctIdToCount[PctId]
	SumCount += Count
assert SumCount <= Total

if REPORT:
	s = "PctId"
	for Rank in Ranks:
		s += "\t" + Rank
	print s

RankToNovelTotal = {}
for Rank in Ranks:
	RankToNovelTotal[Rank] = 0

MinPctId = 100
for PctId in range(100, 0, -1):
	if PctId in PctIdToRankToPCR.keys() and PctId in PctIdToCount.keys():
		MinPctId = PctId

for PctId in range(100, MinPctId-1, -1):
	s = str(PctId)
	Count = PctIdToCount[PctId]
	if PctId == MinPctId:
		Count += Total - SumCount
	for Rank in Ranks:
		PLCR = PctIdToRankToPLCR[PctId][Rank]
		PCR = PctIdToRankToPCR[PctId][Rank]
		Known = Count*PCR
		Novel = int(float(Count)*(1 - PCR))
		RankToNovelTotal[Rank] += Novel
		s += "\t%u" % (Novel)
	if REPORT:
		print s

s = ""
for Rank in Ranks:
	if s != "":
		s += "\t"
	s += Rank
print s

s = ""
for Rank in Ranks:
	Novel = RankToNovelTotal[Rank]
	Pct = (100.0*Novel)/Total
	if s != "":
		s += "\t"
	s += "%u (%.1f%%)" % (Novel, Pct)
Fields = IdDistFileName.split('/')
Name = Fields[len(Fields)-1]
s += "\t" + Name
print s
