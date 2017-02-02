#some import
import numpy
import itertools
import random
import math
import os.path
import urllib2
import sys

#global variables
n = None
N = None
Dk = None
Du = None
Dc = None
e = None
r = None
x = None
xa = None
ASW = None
EXC = None
INS = None
pASW = None
pEXC = None
pINS = None

#help function
def help():
	print("""
VARIABLES
---------
n                     permutation size
N                     total number of permutations
Dk                    diameter with Kendall's tau distance
Du                    diameter with Ulam distance
Dc                    diameter with Cayley distance
e                     identity permutation
r                     identity reversed permutation
x                     a random permutation
xa                    antipodal permutation of x with respect to exchanges
ASW                   set of adjacent swap generators (in tuple normal form)
EXC                   set of exchange generators (in tuple normal form)
INS                   set of insertion generators (in tuple normal form)
pASW                  ASW in permutation form
pEXC                  EXC in permutation form
pINS                  INS in permutation form

FUNCTIONS
---------
reset(n)              re-initialize environment with a new permutation size
fact(n)               factorial of n
inv(x)                inverse of x
rev(x)                reverse of x
compl(x)              complement of x
dot(x,y)              composition x*y
prand()               a random permutation
isPerm(x)             check if the list x is a permutation
asw(i)                adjacent swap (i,i+1) [callable also with a tuple arg]
exc(i,j)              exchange (i,j)        [callable also with a tuple arg]
ins(i,j)              insertion (i,j)       [callable also with a tuple arg]
asw_nf(i)             tuple normal form for adjacent swap (i,j) with j=i+1
exc_nf(t)             tuple normal form for exchange (t[0],t[1]) with t[1]>t[0]
ins_nf(t)             tuple normal form for insertion (t[0],t[1]) with t[1]>t[0] only for case t[1]=t[0]+1 
aswToTuple(x)         convert an adjacent swap from permutation form to tuple form
excToTuple(x)         convert an exchange from permutation form to tuple form
insToTuple(x)         convert an insertion from permutation form to tuple form
swap(x,i,j)           swap items at position i and j in x (inplace)
insert(x,i,j)         shift item at position i to position j in x (inplace)
dk(x,y)               Kendall's tau distance between x and y
du(x,y)               Ulam distance between x and y
dc(x,y)               Cayley distance between x and y
ninver(x)             number of inversion in x
inver(x)              set of inversions in x
ainver(x)             list of adjacent inversions in x
lis(x)                standard lis of x
llis(x)               length of a lis of x
alis(x)               all the lis of x
urlis(x)              unfirom random lis of x
ind(l,x)              indexes of x where values in list l appear
cycles(x)             cycle decomposition of x
ncycles(x)            number of cycles of the decomposition of x
cycleToPerm(c)        build a permutation corresponding to cycle c
prandUnredLis()       random permutation whose lis is "unreducible"
hasUnredLis(x)        check if the lis of x can be reduced or not
lisldsRedCases(x)     print and return lis/lds reduction cases after applying all insertion generators (if additional False is passed it doesnt print)
printLisLdsRedCases(d)print the result of "lisldsRedCases"
bothLisLdsIncrease()  get a random permutation x + ins generator g such that x*g increases both lis and lds lengths
prand1c()             random permutation whith only one cycle
stirling1u(n,k)       [n,k] unsigned stirling number of the 1st kind
npermWithCycles(k)    number of permutations with k cycles
lds(x)                standard lds of x
llds(x)               length of a lds of x
alds(x)               all the lds of x
urlds(x)              uniform random lds of x
mahonian(n,k)         [n,k] mahonian number
npermWithInvers(k)    number of permutations with k inversions
seqA047874(n,k)       [n,k] number of the sequence A047874 (it works only till n=60 and requires the file at https://oeis.org/A047874/b047874.txt or internet)
npermWithLisLength(k) number of permutations with a lis of length k (it works only tille n=60 the file at https://oeis.org/A047874/b047874.txt or internet)
applySeq(x,s)         return x*s[0]*s[1]*...
composeSeq(s)         return s[0]*s[1]*...
mapAswSeq(s)          from a sequence of ASW tuples return a sequence of permutations
mapExcSeq(s)          from a sequence of EXC tuples return a sequence of permutations
mapInsSeq(s)          from a sequence of INS tuples return a sequence of permutations
randbs(x)             return a sequence of ASW tuples that sorts x
randDecAsw(x)         return a ASW decomposition of x
randss(x)             return a sequence of EXC tuples that sorts x
randmergess(x)        return a sequence of EXC tuples that UNsorts x
randDecExc(x)         return a EXC decomposition of x
randis(x,randlis)     return a sequence of INS tuples that sorts x (UNIFORM STEP NOT IMPLEMENTED) (the randlis function as parameter is optional)
randDecIns(x,randlis) return a INS decomposition of x (see randis)
checkAllInsDiamRev()  return true if for all permutations x the Ulam distance between x and rev(x) equals the Ulam diameter
ssort(x)              return the sequence of EXC using classical selection sort
expInertia(nexp=1000) write how many inertia anomalies are over nexp random experiments
	""")

#test function
#def test():
#	pass

#default permutation size
DEFAULT_n = 10

#reset/init functions (with global variables declaration)
def reset(size=DEFAULT_n):
	#global variables
	global n,N,Dk,Du,Dc,e,r,ASW,EXC,INS,pASW,pEXC,pINS
	#permutation size
	n = size
	#total number of permutations
	N = fact(n)
	#diameters
	Dk = n*(n-1)/2
	Du = n-1
	Dc = n-1
	#useful permutations
	e = range(n)
	r = e[::-1]
	x = prand()
	xa = applySeq(x,mapExcSeq(randmergess(x)))
	#generators sets
	ASW = set()
	for i in range(n-1):
		ASW.add((i,i+1))
	pASW = sorted(map(lambda p : asw(p), ASW))
	EXC = set()
	for i in range(n):
		for j in range(i+1,n):
			EXC.add((i,j))
	pEXC = sorted(map(lambda p : exc(p), EXC))
	INS = set()
	for i in range(n):
		for j in filter(lambda j : j!=i and j!=i-1,range(n)):
			INS.add((i,j))
	pINS = sorted(map(lambda p : ins(p), INS))
	#copy variables to the main module scope
	import __main__
	__main__.n = n
	__main__.N = N
	__main__.Dk = Dk
	__main__.Du = Du
	__main__.Dc = Dc
	__main__.e = e
	__main__.r = r
	__main__.x = x
	__main__.xa = xa
	__main__.ASW = ASW
	__main__.pASW = pASW
	__main__.EXC = EXC
	__main__.pEXC = pEXC
	__main__.INS = INS
	__main__.pINS = pINS
			
init = reset

#some basic functions
def fact(n):
	return math.factorial(n)

def inv(x):
	z = [None]*n
	for i in range(n):
		z[x[i]] = i
	return z

def rev(x):
	return x[::-1]
    
def dot(x,y):
	return [x[v] for v in y]
    
def prand():
	return numpy.random.permutation(n).tolist()
	
def isPerm(x):
	return sorted(x)==e
	
def compl(x):
	return [n-1-v for v in x]

#generators to permutation functions
def asw(g1,g2=-1):
	if type(g1) is tuple:
		return exc(g1[0],g1[0]+1)
	else:
		return exc(g1,g1+1)

def exc(g1,g2=-1):
	if type(g1) is tuple:
		i,j = g1
	else:
		i,j = g1,g2
	z = e[:]
	z[i],z[j] = z[j],z[i]
	return z
    
def ins(g1,g2=-1):
	if type(g1) is tuple:
		i,j = g1
	else:
		i,j = g1,g2
	if i<j:
		return range(i) + range(i+1,j+1) + [i] + range(j+1,n)
	else:
		return range(j) + [i] + range(j,i) + range(i+1,n)

def asw_nf(t):
	if type(t) is not tuple:
		return (t,t+1)
	return exc_nf(t)

def exc_nf(t):
	return tuple(sorted(t))

def ins_nf(t):
	return tuple(sorted(t)) if t[0]==t[1]+1 else t
	
def aswToTuple(x):
	t = excToTuple(x)
	if t[1]!=t[0]+1:
		print("It is not an adjacent swap!!!")
	return (t[0],t[0]+1)

def excToTuple(x):
	diff = [i==x[i] for i in range(n)]
	if diff.count(False)!=2:
		print("It is not an exchange!!!")
	return tuple([i for i,v in enumerate(diff) if not v])
	
def insToTuple(x):
	diff = [i==x[i] for i in range(n)]
	if diff.count(False)<2:
		print("It is not an insertion!!!")
	first,last = diff.index(False),n-1-diff[::-1].index(False)
	if any(diff[first:last]):
		print("It is not an insertion!!!")
	if x[first]==first+1: #i<j
		if x[first:last-1]!=range(first+1,last):
			print("It is not an insertion!!!")
		return (first,last)
	else: #i>j
		if x[first+1:last]!=range(first,last-1) or x[first]!=last:
			print("It is not an insertion!!!")
		return (last,first)
	
	
#swap and insert inplace
def swap(x,i,j=-1):
	if j==-1:
		j = i+1
	x[i],x[j] = x[j],x[i]
    
def insert(x,i,j):
	t = x[i]
	del x[i]
	x.insert(j,t)
	
#distances
def dk(x,y):
	return ninver(dot(inv(y),x))
	
def du(x,y):
	return n-llis(dot(inv(y),x))
	
def dc(x,y):
	return n-ncycles(dot(inv(y),x))
	
#inversion function
def ninver(x):
	return len(inver(x))

def inver(x):
	return set([(i,j) for i,j in itertools.combinations(range(n),2) if x[i]>x[j]])
	
def ainver(x):
	return [(i,i+1) for i in range(n-1) if x[i]>x[i+1]]

#lis functions
def lis(x):
	#see http://rosettacode.org/wiki/Longest_increasing_subsequence#Python
	X = x[:]
	N = len(X)
	P = [0 for i in range(N)]
	M = [0 for i in range(N+1)]
	L = 0
	for i in range(N):
		lo = 1
		hi = L
		while lo <= hi:
			mid = (lo+hi)//2
			if (X[M[mid]] < X[i]):
				lo = mid+1
			else:
				hi = mid-1
		newL = lo
		P[i] = M[newL-1]
		M[newL] = i
		if (newL > L):
			L = newL
	S = []
	k = M[L]
	for i in range(L-1, -1, -1):
		S.append(X[k])
		k = P[k]       
	return S[::-1]
	
def llis(x):
	return len(lis(x))
	
def alis(x):
	#see http://stackoverflow.com/questions/9554266/finding-all-possible-longest-increasing-subsequence?rq=1
	count = [1]*n
	def longestIncreaseSubsequence(seq):
		n = len(seq)
		for i in range(1,n):
			maxi = 0
			for j in range(i-1,-1,-1):
				if seq[j]<seq[i]:
					maxi = max(maxi,count[j])
			count[i] = maxi + 1
		maxi = 0
		for i in range(len(count)):
			if count[i]>maxi:
				maxi = count[i]
		return maxi
	def allLIS(a,k,count,arr,maxi,result):
		if k==maxi:
			lista = []
			for i in range(maxi,0,-1):
				lista.append(arr[a[i]])
			result.append(lista)
		else:
			k = k+1
			candidates = [None]*len(arr)
			ncandidate = 0
			for i in range(a[k-1],-1,-1):
				if count[i]==maxi-k+1 and (arr[i]<arr[a[k-1]] or count[i]==maxi):
					candidates[ncandidate] = i
					ncandidate = ncandidate + 1
			for i in range(0,ncandidate):
				a[k] = candidates[i]
				allLIS(a,k,count,arr,maxi,result) 
	maxi = longestIncreaseSubsequence(x)
	a = [None]*(maxi+1)
	a[0] = len(x)-1
	result = []
	allLIS(a,0,count,x,maxi,result)
	return result
	
def urlis(x):
	return random.choice(alis(x))
	
def ind(l,x):
	return [x.index(v) for v in l]
	
#cycles functions
def cycles(perm):
	#see https://gist.github.com/begriffs/2211881
	remain = set(perm)
	result = []
	while len(remain) > 0:
		n = remain.pop()
		cycle = [n]
		while True:
			n = perm[n]
			if n not in remain:
				break
			remain.remove(n)
			cycle.append(n)
		result.append(cycle)
	return result
	
def ncycles(x):
	return len(cycles(x))

def cycleToPerm(c):
	z = range(n)
	for k in range(len(c)-1):
		i = c[k]
		j = c[k+1]
		z[i] = j
	z[c[-1]] = c[0]
	return z
	
#lis reduction functions
def prandUnredLis():
	while True:
		x = prand()
		lx = llis(x)
		flag = True
		for pins in pINS:
			y = dot(x,pins)
			ly = llis(y)
			if ly<lx:
				flag = False
		if flag:
			return x
			
def hasUnredLis(x):
	lx = llis(x)
	for pins in pINS:
		y = dot(x,pins)
		if llis(y)<lx:
			return False
	return True
	
def lisldsRedCases(x,verbose=True):
	r = { "<<":0, "<=":0, "<>":0, "=<":0, "==":0, "=>":0, "><":0, ">=":0, ">>":0, "other":0 }
	l1x,l2x = llis(x),llds(x)
	for g in pINS:
		y = dot(x,g)
		l1y,l2y = llis(y),llds(y)
		if l1y==l1x-1 and l2y==l2x-1:
			r["<<"] += 1
		elif l1y==l1x-1 and l2y==l2x:
			r["<="] += 1
		elif l1y==l1x-1 and l2y==l2x+1:
			r["<>"] += 1
		elif l1y==l1x and l2y==l2x-1:
			r["=<"] += 1
		elif l1y==l1x and l2y==l2x:
			r["=="] += 1
		elif l1y==l1x and l2y==l2x+1:
			r["=>"] += 1
		elif l1y==l1x+1 and l2y==l2x-1:
			r["><"] += 1
		elif l1y==l1x+1 and l2y==l2x:
			r[">="] += 1
		elif l1y==l1x+1 and l2y==l2x+1:
			r[">>"] += 1
		else:
			r["other"] += 1
	if verbose:
		printLisLdsRedCases(r)
	return r
	
def printLisLdsRedCases(d):
	print("ID")
	print("<<    : "+(str(d["<<"] if "<<" in d else 0)))
	print("<=    : "+(str(d["<="] if "<=" in d else 0)))
	print("<>    : "+(str(d["<>"] if "<>" in d else 0)))
	print("=<    : "+(str(d["=<"] if "=<" in d else 0)))
	print("==    : "+(str(d["=="] if "==" in d else 0)))
	print("=>    : "+(str(d["=>"] if "=>" in d else 0)))
	print("><    : "+(str(d["><"] if "><" in d else 0)))
	print(">=    : "+(str(d[">="] if ">=" in d else 0)))
	print(">>    : "+(str(d[">>"] if ">>" in d else 0)))
	print("other : "+(str(d["other"] if "other" in d else 0)))
	
def bothLisLdsIncrease():
	while True:
		x = prand()
		l1x,l2x = llis(x),llds(x)
		for g in pINS:
			y = dot(x,g)
			l1y,l2y = llis(y),llds(y)
			if l1y>l1x and l2y>l2x:
				return [x,g]
	
#random permutation with only 1 cycle
def prand1c():
	x = [None]*n
	c = range(1,n)
	i = 0
	while c:
		j = random.choice(c)
		c.remove(j)
		x[i] = j
		i = j
	x[i] = 0
	return x
	
#cycle distribution functions
def stirling1u(n,k):
	#stirling number of the 1st kind unsigned
	if n==0 and k==0:
		return 1
	if n==0 or k==0:
		return 0
	return (n-1)*stirling1u(n-1,k) + stirling1u(n-1,k-1)

def npermWithCycles(k):
	return stirling1u(n,k)
	
#lds functions
def lds(x):
	return compl(lis(compl(x)))
	
def llds(x):
	return len(lds(x))
	
def alds(x):
	return [compl(l) for l in alis(compl(x))]
	
def urlds(x):
	return compl(urlis(compl(x)))

#inversion distribution functions
def mahonian(n,k):
	#see http://stackoverflow.com/questions/19372991/number-of-n-element-permutations-with-exactly-k-inversions
	def mahonian_row(n):
		i = 1
		result = [1]
		while i < n:
			prev = result[:]
			result = [0] * int(1 + ((i + 1) * 0.5) * (i))
			m = [1] * (i + 1)
			for j in range(len(m)):
				for k in range(len(prev)):
					result[k+j] += m[j] * prev[k]
			i = i + 1
		return result
	return mahonian_row(n)[k]
	
def npermWithInvers(k):
	return mahonian(n,k)
	
#lis length distribution function
def seqA047874(n,k):
	#see https://oeis.org/A047874 and https://oeis.org/A047874/b047874.txt
	if n>60:
		print("Impossible to compute this value for n greater than 60")
	lineno = n*(n-1)/2 + k
	fn = "b047874.txt"
	if os.path.exists(fn):
		with open(fn,"r") as f:
			for line in f:
				if int(line.split()[0])==lineno:
					return int(line.split()[1])
	else:
		print "Trying to read the file from web https://oeis.org/A047874/b047874.txt"
		un = "https://oeis.org/A047874/b047874.txt"
		txt = urllib2.urlopen(un).read().split("\n")
		for line in txt:
			if int(line.split()[0])==lineno:
					return int(line.split()[1])
	return -1

def npermWithLisLength(k):
	return seqA047874(n,k)
	
#randomized sorting algorithms
def applySeq(x,s):
	z = x[:]
	for y in s:
		z = dot(z,y)
	return z
	
def composeSeq(s):
	return applySeq(e,s)
	
def mapAswSeq(s):
	return map(lambda p : asw(p), s)
	
def mapExcSeq(s):
	return map(lambda p : exc(p), s)
	
def mapInsSeq(s):
	return map(lambda p : ins(p), s)

def randbs(x):
	y = x[:]
	s = []
	ai = ainver(y)
	while len(ai)>0:
		sw = random.choice(ai)
		swap(y,sw[0],sw[1])
		ai = ainver(y)
		s.append(sw)
	return s
	
def randDecAsw(x):
	return randbs(inv(x))
	
def randss(x):
	y = x[:]
	s = []
	cyc = cycles(y)
	while len(cyc)<n:
		cyc = filter(lambda c : len(c)>1,cyc)
		q = list(numpy.cumsum([len(c)*(len(c)-1)/2 for c in cyc]))
		tot = q[-1]
		r = random.randint(0,tot-1)
		for i in range(len(cyc)):
			if r<q[i]:
				c = i
		c = cyc[c]
		i = random.choice(c)
		c.remove(i)
		j = random.choice(c)
		s.append(exc_nf((i,j)))
		swap(y,i,j)
		cyc = cycles(y)
	return s
	
def randmergess(x):
	y = x[:]
	s = []
	cyc = cycles(y)
	while len(cyc)>1:
		w = list(numpy.cumsum([len(cyc[k])*(n-len(cyc[k])) for k in range(len(cyc))]))
		r = random.randint(0,w[-1]-1)
		for c1 in range(len(cyc)):
			if r<w[c1]:
				break
		i = random.choice(cyc[c1])
		del cyc[c1]
		w = list(numpy.cumsum(map(lambda c : len(c),cyc)))
		r = random.randint(0,w[-1]-1)
		for c2 in range(len(cyc)):
			if r<w[c2]:
				break
		j = random.choice(cyc[c2])
		s.append(exc_nf((i,j)))
		swap(y,i,j)
		cyc = cycles(y)
	return s

def randDecExc(x):
	return randss(inv(x))
	
def randis(x,randlis=urlis):
	y = x[:]
	s = []
	lis = randlis(y)
	while len(lis)<n:
		u = [i for i in range(n) if i not in lis]
		i = random.choice(ind(u,y))
		ival = y[i]
		for b in range(len(lis)):
			if lis[b]>ival:
				break
		else:
			b = len(lis)
		if b==0:
			a,b = 0,y.index(lis[0])
		elif b==len(lis):
			a,b = y.index(lis[-1]),n-1
		else:
			a,b = y.index(lis[b-1]),y.index(lis[b])
		if a==b:
			j = a
		elif i<a:
			j = random.randint(a,b-1)
		elif i>b:
			j = random.randint(a+1,b)
		else:
			j = None
			print("Problem with randis")
		s.append(ins_nf((i,j)))
		lis.append(ival)
		lis = sorted(lis)
		insert(y,i,j)
		if lis not in alis(y):
			print("BIG PROBLEM")
	return s
	
def decInsSeq(x,randlis=urlis):
	return randis(inv(x),randlis)
			
def checkAllInsDiamRev():
	#return true if for all permutations x the Ulam distance between x and rev(x) equals the Ulam diameter
	#return false otherwise
	for p in itertools.permutations(e):
		x = list(p)
		r = rev(x)
		if du(x,r)!=Du:
			return False
	return True
		
def ssort(x):
	y = x[:]
	s = []
	for j in range(0,n-1):
		imin = j
		for i in range(j+1,n):
			if y[i]<y[imin]:
				imin = i
		if imin!=j:
			t = y[j]
			y[j] = y[imin]
			y[imin] = t
			s.append(exc_nf((j,imin)))
	return s

def expInertia(nexp=1000):
	anomalies = 0
	for i in xrange(nexp):
		x = prand()
		dx = randDecAsw(x)
		y = dot(x,asw(dx[0]))
		wx = ninver(x)
		wy = ninver(y)
		if wy!=wx+1:
			anomalies += 1
	print "Anomalies: " + str(anomalies) + " / " + str(nexp)
	
	


#init the environment and print usage
init()
#test()
help()