#!/usr/bin/python2

import sys
import utax2
import die

Levels = [ 'd', 'p', 'c', 'o', 'f', 'g' ]

File = open(sys.argv[1])
while 1:
	Line = File.readline()
	if len(Line) == 0:
		break

# Tabbed, 2 fields, x: inserted by mothur_make_taxtrainfiles.py
# QueryLabel d:Bacteria;p:Proteobacteria;c:Epsilonproteobacteria;o:Campylobacterales;f:Helicobacteraceae;g:Helicobacter;unclassified;
	Fields = Line.strip().split('\t')
	assert len(Fields) == 2

	QueryLabel = Fields[0]
	MotPred = Fields[1]
	if not MotPred.endswith(";"):
		die.Die("doesn't end with ; '" + MotPred + "'")

	PredStr = ""
	Fields2 = MotPred.split(';')
	for Field in Fields2:
		if Field == "" or Field == "unclassified" or Field == "unknown":
			continue
		if PredStr != "":
			PredStr += ","
		PredStr += Field
		
	print QueryLabel + '\t' + PredStr
