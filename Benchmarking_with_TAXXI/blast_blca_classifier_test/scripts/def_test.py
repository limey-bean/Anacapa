#!/usr/bin/python2

FileName_NameCounts='/Users/limeybean/Downloads/test.txt'

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
