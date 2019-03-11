#!/usr/bin/python2

import sys
import die

FileName = sys.argv[1]
SHOW_SCORE = False

# AJ698860	k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Sporolactobacillaceae; g__Sporolactobacillus; s__	1.00	3

File = open(FileName)
while 1:
	Line = File.readline()
	if len(Line) == 0:
		break
	if Line.startswith("#"):
		continue

	Fields = Line[:-1].split('\t')
	assert len(Fields) >= 3
#	if len(Fields) != 4:
#		die.Die("Expected 4 fields got %d in '%s'" % (len(Fields), Line[:-1]))

	Label = Fields[0]
	QiimeTax = Fields[1]
	if QiimeTax == "unassigned" or Fields[2] == "None":
		Score = 0.0
	else:
		try:
			Score = float(Fields[2])
		except:
			print sys.stderr, "Warning -- bad score '%s'" % Fields[2]
			Score = 0.0

	if QiimeTax == "Unassigned" or QiimeTax == "No blast hit":
		print Label + "\t*"
		continue

	Fields2 = QiimeTax.split("; ")
	if Fields2[0] == "No blast hit":
		print Label + "\t*"
		continue
	
	Tax = ""
	for Field in Fields2:
		if Field[1:3] != "__":
			print >> sys.stderr
			print >> sys.stderr, Line
			die.Die("expected __, Field='%s'" % Field)
		c = Field[0]
		Name = Field[3:]
		if Name == "":
			break
		if Tax != "":
			Tax += ","
		if c == "k":
			c = "d"
		Tax += c + ":" + Name
		if SHOW_SCORE:
			Tax += "(%.2f)" % Score

	s = Label
	s += "\t" + Tax
	print s
