from sage.all import *
import sys
from PrimesList import *
from Legendre import *
from EuclideanAlgorithm import *

def CFRAC(N): # initialize variables
	k = 1
	a = 2 * floor(sqrt(N))
	y = floor(sqrt(N))
	z = 1
	Epair = (1, 0)
	Fpair = (0, 1)
	return kloop(N, k, a, y, z, Epair, Fpair)

def kloop(N, k, a, y, z, Epair, Fpair): # find l + 1 sets of Ak, Ak**2 pairs that are base-smooth
	D = k * N # moke base using D = k * N so that we can make a longer period if k = 1 is too short
	baselist, q = collectbase(D) # find base elements
	l = len(baselist) # find how many base elements there are
	looplen = calclooplen(D) # find how many elements before the continued fraction loops
	counter = 1 # counts how many elements have been calculated
	Q = [] # will store the pairs [Ak, abs(Ak**2)] such that Ak**2 can be factored over baselist (all mod N)
	R = [] # will store the exponents of the factorization of each Ak2 in Q
	while (len(Q) < l + 1) and (counter <= looplen): # need l + 1 equations to ensure linear dependence
		a, y, z, Epair, Fpair, Ak, Bk, Ak2 = nextAk2(D, a, y, z, Epair, Fpair) # calculate next set of numbers
		if bsmoothcheck(Ak2, q) == 1:
			Q.append([Ak, abs(Ak2)]) # add set of numbers to list
			explist = basefactorize(Ak2, baselist) # calculate set of exponents
			R.append(([len(R)], explist)) # add set of exponents to list
		counter += 1
	dep = lindep(N, Q, R) # factorization or None if we cannot yet find a factorization
	if dep != None:
		return dep # the factorization
	else: # try again with different multiplier
		k = nextk(k)
		return kloop(N, k, (2 * floor(sqrt(N))), floor(sqrt(N)), 1, (1, 0), (0, 1))

def nextk(k): # calculates next smallest squarefree multiplier k
	k += 1
	while moebius(k) == 0:
			k += 1
	return k

def calclooplen(D): # calculates the length before the continued fraction loops
	F = QuadraticField(D)
	sqrtDformatted = F(sqrt(D))
	cf = continued_fraction(sqrtDformatted)
	preperiod = cf.preperiod()
	period = cf.period()
	return len(preperiod) + len(period)

def nextAk2(D, a, y, z, Epair, Fpair): # finds next set of numbers needed using the continued fraction expansion
	rootD = sqrt(D)
	a, y, z, Epair, Fpair, Ak, Bk = nextfraction(D, a, y, z, Epair, Fpair)
	Ak2 = (Ak**2) % D
	if Ak2 > 2 * rootD:
		Ak2 -= D
	elif Ak2 < - 2 * rootD:
		Ak2 += D
	return a, y, z, Epair, Fpair, Ak, Bk, Ak2

def collectbase(D): # find base for given value of D
	b = floor(e**(sqrt(ln(D) * ln(ln(D))) / 2)) # find bound for base primes
	if b < 31: # never have too small a base
		b = 31
	plist = [2] # collect base primes N <= b such that Kroencker(N, p) = 1
	pindex = 1 # index 1 will start the loop at the second prime, 3
	p = primelist[pindex]
	q = 2 # the product of the primes in the base
	while p <= b:
		if Kronecker(D, p) == 1: # only use primes for which D is a square mod p
			plist.append(p)
			q *= p
		pindex += 1
		p = primelist[pindex] # go to next prime
	return plist, q

def nextfraction(N, a, y, z, Epair, Fpair): # find next continued fraction
	yk = a * z - y
	zk = floor((N - yk**2) / z)
	ak = floor((floor(sqrt(N)) + yk) / zk)
	Epair = (Epair[1], (ak * Epair[1] + Epair[0]) % N)
	Fpair = (Fpair[1], (ak * Fpair[1] + Fpair[0]) % N)
	Ak = (Epair[0] + floor(sqrt(N)) * Fpair[0]) % N
	Bk = Fpair[0]
	return ak, yk, zk, Epair, Fpair, Ak, Bk

def bsmoothcheck(t, q): # returns 1 if t is bsmooth over base and 0 otherwise
	while True:
		g = gcd(t, q)
		if g == 1:
			return 0
		else:
			while t % g == 0: # divide out all multiples of g
				t = t / g
			if abs(t) == 1:
				return 1

def basefactorize(t, baselist): # factorize over baselist and return exponents
	l = len(baselist)
	explist = [0] * l
	for x in range(0, l):
		if t == 1:
			break
		while (t % baselist[x] == 0):
			t = t / baselist[x]
			explist[x] += 1
	return explist

def lindep(N, Q, R): # find dependency, return either factorization or None
	if R == []:
		return None
	S = R.copy()
	width = len(S[0][1])
	height = len(S)
	finishedheight = 0 # how many rows have been put at top of the matrix and are therefore done
	for item in S: # take all vectors mod 2
		for s in range(0, width):
			item[1][s] = item[1][s] % 2
	for w in range(0, width): # iterate through width
		for h in range(finishedheight, height): # iterate through height below finished rows
			if S[h][1][w] == 1: # find row that contains 1 in the current column
				row = S[h] # save row
				S.remove(S[h]) # remove row
				S.insert(finishedheight, row) # insert row at top of unfinished rows, bottom of finished rows
				for lower in range(finishedheight + 1, height): # add current row to each lower row
					if S[lower][1][w] == 1:
						S[lower] = (S[lower][0] + S[finishedheight][0], [(S[lower][1][j] + S[finishedheight][1][j]) % 2 for j in range(0, width)])
						if S[lower][1] == [0] * width: # if any row has all 0s (meaning all rows multiplied together give all even exponents), check if it factors
							check = checknumbers(N, Q, R, S[lower][0])
							if check != None:
								return check # return a factorization
				finishedheight += 1
				break
	return None

def checknumbers(N, Q, R, rows): # multiply the appropriate numbers together, see if they will give a nontrivial solution
	x2 = 1
	y2 = 1
	for row in rows:
		Ak, Ak2 = Q[row]
		x2 = (x2 * Ak) % N # multiply all the Aks together mod N
		y2 = (y2 * Ak2) # multiply all the Ak**2s together
	y = sqrt(y2)
	if ((x2 % N) != (y % N)) and ((x2 % N) != (- y % N)): # check that x2 != +- y mod N
		return factorize(N, x2, y)
	else: # if x2 != +- y mod N we would get a trivial solution, so try new set of numbers
		return None

def factorize(N, x2, y): # find factorization of N
	r0 = x2 + y
	f0 = gcd(r0, N)
	return f0, int(N / f0)