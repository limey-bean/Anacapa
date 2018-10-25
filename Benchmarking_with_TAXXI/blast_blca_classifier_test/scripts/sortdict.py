d = {}
Order = []
Keys = []

def VectToDict(v):
	d = {}
	for x in v:
		try:
			d[x] += 1
		except:
			d[x] = 1
	return d

def Cmp__(i, j):
	global d, Keys

	ki = Keys[i]
	kj = Keys[j]

	ni = d[ki]
	nj = d[kj]

	if ni < nj:
		return 1
	elif ni > nj:
		return -1
	return 0

def GetOrder(Dict):
	global d, Order, Keys

	d = Dict
	Keys = d.keys()
	N = len(Keys)
	Order = range(0, N)
	Order.sort(Cmp__)
	return Order

def IncCount(Dict, Key, n=1):
	try:
		Dict[Key] += n
	except:
		Dict[Key] = n
	return Dict[Key]

def IncCount2(Dict, Key1, Key2, n=1):
	try:
		x = Dict[Key1]
	except:
		Dict[Key1] = {}
	IncCount(Dict[Key1], Key2, n)

def GetCount(Dict, Key):
	try:
		return Dict[Key]
	except:
		return 0

def GetCount2(Dict, Key1, Key2):
	try:
		return Dict[Key1][Key2]
	except:
		return 0

def GetRanks(Dict):
	Order = GetOrder(Dict)
	Keys = Dict.keys()
	KeyToRank = {}
	N = len(Order)
	for Rank in range(0, N):
		i = Order[Rank]
		Key = Keys[i]
		KeyToRank[Key] = Rank
	return KeyToRank

def CmpVec__(i, j):
	global g_Vec
	if g_Vec[i] < g_Vec[j]:
		return 1
	elif g_Vec[i] > g_Vec[j]:
		return -1
	return 0

def GetOrderDesc(Vec):
	global g_Vec
	g_Vec = Vec
	Order = range(0, len(Vec))
	Order.sort(CmpVec__)
	return Order