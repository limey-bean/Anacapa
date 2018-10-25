#!/usr/bin/python2

import sys
import die

FileName = sys.argv[1]

# AJ698860	k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Sporolactobacillaceae; g__Sporolactobacillus; s__	1.00	3

def FixQuotes(s):
	if not s.startswith('"'):
		return s
	assert s.endswith('"')
	s = s[1:-1]
	s = s.replace('""', '"')
	return s

File = open(FileName)
Line = File.readline()
assert Line.startswith("Feature ID\tTaxon")

while 1:
	Line = File.readline()
	if len(Line) == 0:
		break
	if Line.startswith("#"):
		continue

	Fields = Line[:-1].split('\t')
	assert len(Fields) >= 3

	Label = FixQuotes(Fields[0])
	QiimeTax = FixQuotes(Fields[1])

	Tax = ""
	Fields2 = QiimeTax.split("; ")
	for Field in Fields2:
		if Field.lower() == "unassigned":
			continue
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

	s = Label
	s += "\t" + Tax
	print(s)
