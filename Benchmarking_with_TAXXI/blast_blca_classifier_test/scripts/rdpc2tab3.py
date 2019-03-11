#!/usr/bin/python2

# Like rdpc2tab2.py except outputs c:Name instead of Name.

import sys
import die

FileName = sys.argv[1]
Cutoff = -1.0
if len(sys.argv) > 2:
	Cutoff = float(sys.argv[2])

File = open(FileName)

while 1:
	Line = File.readline()
	if len(Line) == 0:
		break

	Fields = Line.strip().split('\t')
	n = len(Fields)
	# two tabs after query! yuck
	assert((n-2)%3 == 0)
	QueryLabel = Fields[0]
	TaxStr = ""
	for i in range(5, n, 3):
		Name = Fields[i]
		LevelName = Fields[i+1]
		if LevelName.startswith("sub"):
			continue
		LevelChar = LevelName[0].lower()
		P = float(Fields[i+2])*100.0
		assert P >= 0.0 and P <= 100.0
		if P < Cutoff:
			break
		if TaxStr != "":
			TaxStr += ','
		sP = "%.0f" % P
		if Cutoff < 0:
			TaxStr += LevelChar + ":" + Name + '(' + sP + ')'
		else:
			TaxStr += LevelChar + ":" + Name
#	TaxStr += ';'
	print QueryLabel + '\t' + TaxStr
