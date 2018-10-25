import sys
import die
import sortdict
import re

LevelChars = "rkdpcofgsq"
FakeNameCounter = 0

def GetParentLevelChar(LevelChar):
	if LevelChar == LevelChars[0]:
		return "-"
	n = LevelChars.find(LevelChar)
	assert n >= 0
	return LevelChars[n-1]

def SplitLevelIsTrue(SplitLevelChar, TestLevelChar):
	if len(SplitLevelChar) != 1 or len(TestLevelChar) != 1:
		die.Die("SplitLevelIsTrue(Split=%s, Test=%s)" % (SplitLevelChar, TestLevelChar))

	Splitn = LevelChars.find(SplitLevelChar)
	Testn = LevelChars.find(TestLevelChar)
	assert Splitn >= 0
	assert Testn >= 0
	return Testn <= Splitn

def SplitLevelToTrueLevels(SplitLevel):
	if SplitLevel == "ref":
		return LevelChars	
	s = ""
	for TestLevel in utax2.LevelChars:
		if utax2.SplitLevelIsTrue(SplitLevel, TestLevel):
			s += TestLevel
	return s

def SplitLevelToFalseLevels(SplitLevel):
	if SplitLevel == "ref":
		return ""	
	s = ""
	for TestLevel in utax2.LevelChars:
		if not utax2.SplitLevelIsTrue(SplitLevel, TestLevel):
			s += TestLevel
	return s

def GetIsTrue(LevelChar, tof, TestLevelChar):
	if tof == 't':
		SplitLevel = LevelChar
	elif tof == 'f':
		SplitLevel = GetParentLevelChar(LevelChar)
	else:
		assert False
	return SplitLevelIsTrue(SplitLevel, TestLevelChar)

def CharToLevelName(c):
	if c == 'r':
		return "root"
	if c == 'd':
		return "domain"
	if c == 'k':
		return "kingdom"
	if c == 'p':
		return "phylum"
	if c == 'c':
		return "class"
	if c == 'o':
		return "order"
	if c == 'f':
		return "family"
	if c == 'g':
		return "genus"
	if c == 's':
		return "species"
	if c == 't':
		return "strain"
	die.Die("CharToLevel(%s)" % c)

def LevelToChar(Level):
	if Level == "domain":
		return 'd'
	if Level == "kingdom":
		return 'k'
	if Level == "phylum":
		return 'p'
	if Level == "class":
		return 'c'
	if Level == "order":
		return 'o'
	if Level == "family":
		return 'f'
	if Level == "genus":
		return 'g'
	if Level == "species":
		return 's'
	if Level == "subclass":
		return 'l'
	if Level == "suborder":
		return 'u'
	if Level == 'subfamily':
		return 'b'
	die.Die("LevelToChar(%s)", Level)

def GetLevelName(Path, Level, Default = ""):
	Fields = Path.split(',')
	for Field in Fields:
		if len(Field) < 2 or Field[1] != ':':
			die.Die("Invalid path '%s'" % Path)
		if Field[0] == Level:
			s = Field[2:]
			n = s.find('(')
			if n > 0:
				return s[:n]
			return s
	if Default == "":
		die.Die("GetLevelName, '%c:' not found in path '%s'" % (Level, Path))
	return Default

def GetLevelMatch(Path1, Path2, Level):
	Name1 = GetLevelName(Path1, Level, "?1")
	Name2 = GetLevelName(Path2, Level, "?2")
	return Name1 == Name2

def StripScores(Path):
	return re.sub(r"\([0-9.]*\)", "", Path)

# d:*,p:467,c:508,o:569,f:756,g:863

def ParseNNs(NNs):
	LevelToCount = {}
	Fields = NNs.split(',')
	for Field in Fields:
		assert Field[1] == ':'
		Level = Field[0]
		sCount = Field[2:]
		if sCount == "*":
			sCount = "0"
		LevelToCount[Level] = int(sCount)
	return LevelToCount

def GetPredScore(Path, LevelChar):
	n = Path.find(LevelChar + ":")
	if n < 0:
		return -1.0
	s = Path[n+2:]
	m = s.find("(")
#	if not s[m+1].isdigit():
#		s = s[m+1:]
#		m = s.find("(")
	if m <= 0:
		die.Die("Missing score %c: in %s" % (LevelChar, Path))
	s = s[m+1:]
	m = s.find(")")
	if m <= 0:
		die.Die("Missing ) for level %c: in %s" % (LevelChar, Path))
	s = s[:m]
	return float(s)

def ParsePred(Pred):
	Fields = Pred.split("(")
	n = len(Fields)
	if n == 1:
		die.Die("Missing score in %s" % Pred)
	Name = Fields[0]

	s = Fields[n-1]
	if s.endswith(';'):
		s = s.replace(";", "")
	assert s.endswith(')')
	ss = s[:-1]
	try:
		Score = float(ss)
	except:
		die.Die("utax2.ParsePred, bad pred '%s'" % Pred)
	return Name, Score

def GetPred(PredStr, LevelChar):
	Fields = PredStr.split(',')
	for Field in Fields:
		if len(Field) < 2 or Field[1] != ':':
			die.Die("Invalid predstr '%s'" % PredStr)
		if Field[0] == LevelChar:
			if not Field.endswith(')') and not Field.endswith(");"):
				die.Die("Pred must end with ) '%s'" % Field)
			m = Field.rfind("(")
			if m <= 0:
				die.Die("Missing ( in '%s'" % (PredStr))
			Name = Field[0:m]
			s = Field[m+1:]
			m = s.find(")")
			if m <= 0:
				die.Die("Missing ) for level %c: in %s" % (LevelChar, PredStr))
			s = s[:m]
			Score = float(s)
			return Name, Score
	return "*", -1.0

def GetScoreLo(Level, QueryUniqueWordCount, TopWordCount, LevelToCount):
	TopWordId = (100.0*TopWordCount)/QueryUniqueWordCount
	Count = 0
	try:
		NNCount = LevelToCount[Level]
	except:
		NNCount = 0
	NNWordId = (100.0*NNCount)/QueryUniqueWordCount
	Diff = TopWordId - NNWordId
	Score = (TopWordId + Diff)/2.0
#	Score = (TopWordId + 2*Diff)/3.0
	return Score

def GetTwidScore(Level, QueryUniqueWordCount, TopWordCount, LevelToCount):
	TopWordId = (100.0*TopWordCount)/QueryUniqueWordCount
	return TopWordId

def GetScore(Level, QueryUniqueWordCount, TopWordCount, LevelToCount):
	Scoref = GetScoreLo('f', QueryUniqueWordCount, TopWordCount, LevelToCount)
	Scoreg = GetScoreLo('g', QueryUniqueWordCount, TopWordCount, LevelToCount)
#	return (2*Scoref + Scoreg)/3.0
	TopWordId = (100.0*TopWordCount)/QueryUniqueWordCount
	try:
		NNCount = LevelToCount[Level]
	except:
		NNCount = 0
	NNWordId = (100.0*NNCount)/QueryUniqueWordCount
	Diff = TopWordId - NNWordId
	Score = (TopWordId + Diff)/2.0
#	Score = (TopWordId + 2*Diff)/3.0
	return Score

def GetLKRFromLabel(Label):
	Fields = Label.split(";")
	for Field in Fields:
		if Field.startswith("lkr="):
			assert len(Field) == 5
			return Field[4]
	die.Die("lkr= not found in >" + Label)

def GetPctIdFromLabel(Label):
	Fields = Label.split(";")
	for Field in Fields:
		if Field.startswith("pctid="):
			PctId = int(Field[6:])
			return PctId
	die.Die("pctid= not found in >" + Label)

def GetNamesFromLabel(Label):
	Tax = GetTaxFromLabel(Label)
	Names = Tax.split(',')
	for Name in Names:
		assert len(Name) > 2 and Name[1] == ':'
	return Names

def GetRanksFromLabel(Label):
	Tax = GetTaxFromLabel(Label)
	Names = Tax.split(',')
	Ranks = ""
	for Name in Names:
		assert len(Name) > 2 and Name[1] == ':'
		Rank = Name[0]
		Ranks += Rank
	return Ranks

def GetPredsFromPredStr(PredStr):
	Preds = PredStr.split(',')
	for Pred in Preds:
		assert len(Pred) > 2 and Pred[1] == ':'
	return Preds

def GetDictFromLabel(Label):
	Tax = GetTaxFromLabel(Label)
	Names = Tax.split(',')
	Dict = {}
	for Name in Names:
		if len(Name) <= 2 or Name[1] != ':':
			die.Die("GetDictFromLabel(%s)" % Label)
		Dict[Name[0]] = Name
	return Dict

def GetNameFromTaxStr(Tax, Level):
	if Tax == "":
		return ""
	Names = Tax.split(',')
	for Name in Names:
		if len(Name) <= 2 or Name[1] != ':':
			die.Die("GetNameFromTaxStr(" + Tax + ", " + Level + "), bad name '" + Name + "'")
		if Name[0] == Level:
			return Name
	return ""

def GetNameFromLabel(Label, Level):
	Tax = GetTaxFromLabel(Label)
	Names = Tax.split(',')
	for Name in Names:
		assert len(Name) > 2 and Name[1] == ':'
		if Name[0] == Level:
			return Name
	return ""

def GetNameIndex(Names, Level):
	for i in range(0, len(Names)):
		if Names[i][0] == Level:
			return i
	return -1

def GetTaxFromLabel(Label):
	Fields = Label.split(";")
	for Field in Fields:
		if Field.startswith("tax="):
			Tax = Field[4:]
			return Tax
	die.Die("tax= not found in label >" + Label)

def GetAccFromLabel(Label):
	if Label.find(";_") > 0:
		Fields = Label.split(";_")
	else:
		Fields = Label.split("\t")
	f = Fields[0]
	Fields2 = f.split(";")
	return Fields2[0]

def TaxToVec(Tax):
	return Tax.split(",")

def TaxToVecs(Tax):
	Chars = []
	Names = []
	Fields = Tax.split(",")
	for Field in Fields:
		assert Field[1] == ':'
		Char = Field[0]
		Name = Field[2:]

		Chars.append(Char)
		Names.append(Name)

	return Chars, Names

def TaxToDict(Tax):
	Dict = {}
	Fields = Tax.split(",")
	for Field in Fields:
		assert Field[1] == ':'
		Char = Field[0]
		Name = Field[2:]
		Dict[Char] = Name
	return Dict

def TaxToRDP(Tax, LevelChars): # e.g. LevelChars=dpcofgs
	global FakeNameCounter

	s = "Root"
	Dict = TaxToDict(Tax)
	for LevelChar in LevelChars:
		try:
			Name = Dict[LevelChar]
		except:
			Name = "?"

		if Name == "?":
			FakeNameCounter += 1
			Name = "FAKE_" + str(FakeNameCounter)
			print >> sys.stderr, Name

		s += ";" + Name
	return s

def TttTaxToTax(TttTax):
	Tax = ""
	Fields = TttTax.split("; ")
	n = len(Fields)
	for i in range(0, n):
		Field = Fields[i]
		if i > 0:
			Tax += ','
		assert Field[1:3] == ':'
		c = Field[0]
		assert c.isalpha()
		Name = FIeld[3:]
		if Name == "":
			return Tax
		Tax += c + ':' + Name
	return Tax

def LevelCharToRankIndex(c):
	n = LevelChars.find(c)
	if n < 0:
		die.Die("Unknown level char '%c'" % c)
	return n

def IsHigher(Level1, Level2):
	Ix1 = LevelCharToRankIndex(Level1)
	Ix2 = LevelCharToRankIndex(Level2)
	return Ix1 < Ix2

def IsKnown(Level, LKR):
	Ix = LevelCharToRankIndex(Level)
	IxLKR = LevelCharToRankIndex(LKR)
	return Ix <= IxLKR

def LevelCharToName(c):
	if c == 's':
		return "Species"
	elif c == 'g':
		return "Genus"
	elif c == 'f':
		return "Family"
	elif c == 'c':
		return "Class"
	elif c == 'o':
		return "Order"
	elif c == 'p':
		return "Phylum"
	elif c == 'd':
		return "Domain"
	elif c == 'k':
		return "Kingdom"
	die.Die("Unknown level " + c)

def TaxToTttTax(Tax, LevelChars):
	LevelCharToName = {}
	Fields = Tax.split(",")
	n = len(Fields)
	for i in range(0, n):
		Field = Fields[i]
		assert Field[1] == ':'
		c = Field[0]
		Name = Field[2:]
		LevelCharToName[c] = Name

	TttTax = ""
	for i in range(0, len(LevelChars)):
		c = LevelChars[i]
		try:
			Name = LevelCharToName[c]
		except:
			Name = ""
		if i > 0:
			TttTax += "; "
		TttTax += (c + "__" + Name)
	return TttTax

def GetLeafNameFromLabel(Label):
	Names = GetNamesFromLabel(Label)
	n = len(Names)
	assert n > 0
	return Names[n-1]

def GetNameWeights(NameToCount):
	Names = NameToCount.keys()
	UniqueCount = len(Names)

	N = 0
	for i in range(0, UniqueCount):
		Name = Names[i]
		n = NameToCount[Name]
		N += n

	NameToWeight = {}
	SumWeight = 0.0
	for i in range(0, UniqueCount):
		Name = Names[i]
		n = NameToCount[Name]
		Weight = float(n)*UniqueCount/float(N)
		NameToWeight[Name] = Weight
		SumWeight += Weight

	r = SumWeight/UniqueCount
	assert r > 0.99 and r < 1.01
	return NameToWeight

def FixName(Name):
	Name = Name.replace(' ', '_')
	Name = Name.replace('(', '[')
	Name = Name.replace(')', ']')
	Name = Name.replace(',', '_')
	Name = Name.replace(';', '_')
	return Name

def GetTopLevelChars(Labels):
	RanksToCount = {}
	for Label in Labels:
		Ranks = GetRanksFromLabel(Label)
		sortdict.IncCount(RanksToCount, Ranks)
	Order = sortdict.GetOrder(RanksToCount)
	RankVec = RanksToCount.keys()
	TopLevelChars = RankVec[Order[0]]
	return TopLevelChars

def GetRankNameFromDict(Dict, Rank, Default = None):
	try:
		Name = Dict[Rank]
	except:
		Name = Default
	return Name

def GetLCR(Tax1, Tax2):
	d1 = TaxToDict(Tax1)
	d2 = TaxToDict(Tax2)
	n = len(LevelChars)
	for i in range(0, n):
		Rank = LevelChars[n-i-1]
		Name1 = GetRankNameFromDict(d1, Rank)
		Name2 = GetRankNameFromDict(d2, Rank)
		if Name1 != None and Name1 == Name2:
			return Rank
	return None

def GetLCRFromLabels(Label1, Label2):
	Tax1 = GetTaxFromLabel(Label1)
	Tax2 = GetTaxFromLabel(Label2)
	LCR = GetLCR(Tax1, Tax2)
	return LCR
