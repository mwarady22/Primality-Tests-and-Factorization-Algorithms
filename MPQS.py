from sage.all import *
import sys
from PrimesList import *
from Legendre import *
from SqrtModp import *

def MPQS(N, f, M):
	baselist, baseproduct, largestbaseindex = collectbase(N, f) # list of base elements, product of all base elements, index in primelist of largest element
	sqrts = [] # will hold the modular square roots of each odd prime in the baselist, with None at index 0
	logs = []# will hold the base 2 log of each prime in the baselist
	for p in baselist: # fill sqrts and logs
		if p == 2:
			sqrts.append(1)
			logs.append(1)
		else:
			sqrts.append(SqrtModp(N, p))
			logs.append(log(p, 2).n())
	qindex, q = chooseq(N, M, largestbaseindex) # find an optimal q and its index in primelist
	if q == None:
		return 'MPQS will not work on ' + str(N)
	plusindex = qindex # index radiating out from q in the positive direction
	minusindex = qindex - 1 # index radiating out from q in the negative direction
	last = 'minus' # indicates which direction the last index was taken from
	while True:
		if minusindex > largestbaseindex: # make sure Kronecker(N, q) == 1
			if Kronecker(N, primelist[minusindex]) == 1:
				break
			else:
				minusindex -= 1
		else:
			last = 'stayplus'
			break
	return newpoly(N, M, baselist, baseproduct, largestbaseindex, sqrts, logs, plusindex, minusindex, last)

def newpoly(N, M, baselist,baseproduct, largestbaseindex, sqrts, logs, plusindex, minusindex, last):
	if last == 'minus': # last q was taken in the negative direction, now use positive directions
		qindex = plusindex
		last = 'plus'
		plusindex += 1
		while plusindex < len(primelist) and Kronecker(N, primelist[plusindex]) != 1:
			plusindex += 1
		if plusindex >= len(primelist):
			last = 'stayminus'
	elif last == 'plus': # last q was taken in the positive direction, now use negative directions
		qindex = minusindex
		last = 'minus'
		minusindex -= 1
		while minusindex > largestbaseindex and Kronecker(N, primelist[minusindex]) != 1:
			minusindex -= 1
		if minusindex < 0:
			last = 'stayplus'
	elif last == 'stayplus': # all q from now on taken in the positive direction
		qindex = plusindex
		plusindex += 1
		while plusindex < len(primelist) and Kronecker(N, primelist[plusindex]) != 1:
			plusindex += 1
	elif last == 'stayminus': # all q from now on taken in the negative direction
		qindex = minusindex
		minusindex -= 1
		while minusindex > largestbaseindex and Kronecker(N, primelist[minusindex]) != 1:
			minusindex -= 1
	if (last == 'stayplus') and (qindex >= len(primelist)): # we have used all available qs
		return 'cannot be solved with given bounds'
	if (last == 'stayminus') and (qindex < 0): # we have used all available qs
		return 'cannot be solved with given bounds'
	q = primelist[qindex]
	a = q**2 # a in a * x**2 + b * x + c
	b = int(mod(N, a).sqrt()) # b in a * x**2 + b * x + c
	c = int((b**2 - N) / a) # c in a * x**2 + b * x + c
	sollist = [] # will hold solutions to a * x**2 + b * x + c mod p for each p in baselist
	for pindex in range(0, len(baselist)):
		p = baselist[pindex]
		sqrtp = sqrts[pindex]
		ainv = inverse_mod(a, p)
		soln0 = (ainv * (sqrtp - b)) % p
		if p != 2:
			soln1 = (ainv * (- sqrtp - b)) % p
			sollist.append([soln0, soln1])
		else: # p == 2 has only one solution
			sollist.append([soln0])
	return sieve(N, M, baselist, baseproduct, largestbaseindex, sqrts, logs, plusindex, minusindex, last, a, b, c, sollist)

def sieve(N, M, baselist, baseproduct, largestbaseindex, sqrts, logs, plusindex, minusindex, last, a, b, c, sollist):
	spos = [0] * (M + 1) # holds 0 and positive indices for sieving
	sneg = [0] * (M + 1) # holds 0 and negative indices for sieving
	for pindex in range(0, len(baselist)):
		p = baselist[pindex]
		sol0 = sollist[pindex][0] # first solution mod p
		indpos = sol0
		if sol0 == 0:
			indneg = 0
		else:
			indneg = - (sol0 - p)
		while indpos < M + 1: # while in range, add logs[pindex]
			spos[indpos] = spos[indpos] + logs[pindex]
			indpos += p # do this every p indices
		while indneg < M + 1: # while in range, add logs[pindex]
			sneg[indneg] = sneg[indneg] + logs[pindex]
			indneg += p # do this every p indices
		if p != 2: # repeat with other solution if p != 2
			sol1 = sollist[pindex][1]
			indpos = sol1
			indneg = - (sol1 - p)
			while indpos < M + 1: # while in range, add logs[pindex]
				if indpos >= 0:
					spos[indpos] = spos[indpos] + logs[pindex]
					indpos += p # do this every p indices
			while indneg < M + 1: # while in range, add logs[pindex]
				if indneg >= 0:
					sneg[indneg] = sneg[indneg] + logs[pindex]
					indneg += p # do this every p indices
	return checksieve(N, M, baselist, baseproduct, largestbaseindex, sqrts, logs, plusindex, minusindex, last, a, b, c, sollist, spos, sneg)

def checksieve(N, M, baselist, baseproduct, largestbaseindex, sqrts, logs, plusindex, minusindex, last, a, b, c, sollist, spos, sneg):
	check = log((M * sqrt(N)), 2).n() - 10 # number over which we will perform bsmoothchecks
	Q = [] # will hold index, evaluated polynomial pairs
	R = [] # will hold the base factorization along with its index in R
	x = ZZ['x'].gen()
	g = a * x**2 + 2 * b * x + c
	for indpos in range(0, len(spos)):
		if spos[indpos] >= check: # if it is likely that gx is bsmooth check if it is
			gx = g(indpos)
			smooth = bsmoothcheck(gx, baseproduct)
			if smooth == 1: # if it is bsmooth, put it in R, Q
				R.append([[len(R)], basefactorize(gx, baselist)])
				Q.append([indpos, gx])
	for indneg in range(1, len(sneg)):
		if sneg[indneg] >= check: # if it is likely that gx is bsmooth check if it is
			gx = g(indneg)
			smooth = bsmoothcheck(gx, baseproduct)
			if smooth == 1: # if it is bsmooth, put it in R, Q
				R.append([[len(R)], basefactorize(gx, baselist)])
				Q.append([- indneg, gx])
	dep = lindep(N, Q, R) # factorization or None if we cannot yet find a factorization
	if dep != None:
		return dep # the factorization
	else:
		return newpoly(N, M, baselist, baseproduct, largestbaseindex, sqrts, logs, plusindex, minusindex, last) # didn't work, try new polynomial

def bsmoothcheck(t, q): # returns 1 if t is bsmooth over base and 0 otherwise
	while True:
		if t == 0:
			return 0
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


def collectbase(D, f): # find base for given value of D with values at most f
	baseproduct = 2
	plist = [2] # collect base primes N <= b such that Kroencker(N, p) = 1
	pindex = 1 # index 1 will start the loop at the second prime, 3
	p = primelist[pindex]
	largestbaseindex = 0
	while p <= f:
		if Kronecker(D, p) == 1: # only use primes for which D is a square mod p
			plist.append(p)
			largestbaseindex = pindex
			baseproduct *= p
		pindex += 1
		p = primelist[pindex] # go to next prime
	return plist, baseproduct, largestbaseindex


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
		xi, yi = Q[row]
		x2 = (x2 * xi) % N # multiply all the inputs to Q that generate the desired congruences
		y2 = (y2 * yi) # multiply all the outputs of Q in the desired congruences
	y = sqrt(y2)
	if ((x2 % N) != (y % N)) and ((x2 % N) != (- y % N)): # check that x2 != +- y mod N
		return factorize(N, x2, y)
	else: # if x2 != +- y mod N we would get a trivial solution, so try new set of numbers
		return None

def factorize(N, x2, y): # find factorization of N
	r0 = x2 + y
	f0 = gcd(r0, N)
	return f0, int(N / f0)

def chooseq(N, M, largestbaseindex):
	close = int(sqrt((sqrt(2 * N)) / M)) # find prime q close to this number which is a quadratic residue mod N
	if close <= primelist[largestbaseindex]: # find q above highest p in baselist
		pindex = largestbaseindex + 1
		kronindex = None
		while True: # find first prime larger than close for which Kronecker gives 1
			if Kronecker(N, primelist[pindex]) == 1:
				kronindex = pindex
				break
			pindex += 1
		return pindex, primelist[pindex]
	else:
		pindex = 0
		while primelist[pindex] < close: # find q larger than close
			pindex += 1
		if Kronecker(N, primelist[pindex]): # if this q has a Kronecker value of 1, return and use it
			return pindex, primelist[pindex]
		else:
			pindexminus = pindex - 1
			pindexplus = pindex + 1
			while True: # find closest q to close that has a Kronecker value of 1
				if (pindexminus < largestbaseindex) and (pindexplus >= len(primelist)):
					return 'cannot find usable q', None
				if pindexminus > largestbaseindex:
					if Kronecker(N, primelist[pindexminus]) == 1:
						return pindexminus, primelist[pindexminus]
					else:
						pindexminus -= 1
				if pindexplus < len(primelist):
					if Kronecker(N, primelist[pindexplus]) == 1:
						return pindexplus, primelist[pindexplus]
					else:
						pindexplus += 1


# print(MPQS(1050703, 100, 1000))
# print(MPQS(15347, 35, 10000))
# print(MPQS(15440, 50, 50))
# print(MPQS(64775585))
# print(MPQS(539873))
# print(MPQS(33569887))