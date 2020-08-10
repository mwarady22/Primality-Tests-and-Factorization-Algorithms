from sage.all import *
import sys
from PrimesList import *
from Legendre import *
from SqrtModp import *

def MPQS(N, f, M):
	print('MPQS')
	baselist, k = collectbase(N, f) # list of base elements, product of all base elements
	print('past baselist')
	sqrts = [] # will hold the modular square roots of each odd prime in the baselist, with None at index 0
	logs = [1]# will hold the base 2 log of each prime in the baselist
	mid = sqrt((sqrt(2 * N) / M))
	print('baselist : ' + str(baselist))
	for p in baselist:
		if p == 2:
			sqrts.append(1)
			logs.append(1)
		else:
			sqrts.append(SqrtModp(N, p))
			logs.append(log(p, 2).n())
	qindex, q = chooseq(N, M, baselist[len(baselist) - 1])
	print('q : ' + str(q))
	plusindex = qindex
	minusindex = qindex - 1
	last = 'minus'
	return newpoly(N, M, baselist, sqrts, logs, plusindex, minusindex, last)

def newpoly(N, M, baselist, sqrts, logs, plusindex, minusindex, last):
	print('newpoly')
	if last == 'minus':
		qindex = plusindex
		last = 'plus'
		plusindex += 1
		if plusindex >= len(primelist):
			last = 'stayminus'
	elif last == 'plus':
		qindex = minusindex
		last = 'minus'
		minusindex -= 1
		if minusindex < 0:
			last = 'stayplus'
	elif last == 'stayplus':
		qindex = plusindex
		plusindex += 1
	elif last == 'stayminus':
		qindex = minusindex
		minusindex -= 1
	if (last == 'stayplus') and (qindex >= len(primelist)):
		return 'cannot be solved with given bounds'
	if (last == 'stayminus') and (qindex < 0):
		return 'cannot be solved with given bounds'
	q = primelist[qindex]
	print('q : ' + str(q))
	a = q**2
	print('a : ' + str(a))
	b0 = mod(N, a).sqrt()
	print('b0 : ' + str(b0))
	b = int(mod(N, a).sqrt())
	print('b : ' + str(b))
	c = int((b**2 - N) / a)
	print('c : ' + str(c))
	sollist = []
	for pindex in range(0, len(baselist)):
		print('p : ' + str(baselist[pindex]))
		p = baselist[pindex]
		sqrtp = sqrts[pindex]
		ainv = inverse_mod(a, p)
		soln0 = (ainv * (sqrtp - b)) % p
		if p != 2:
			soln1 = (ainv * (- sqrtp - b)) % p
			sollist.append([soln0, soln1])
		else:
			sollist.append([soln0])
	return sieve(N, M, baselist, sqrts, logs, indexlist, indexlistindex, a, b, c, sollist)

def sieve(N, M, baselist, sqrts, logs, indexlist, indexlistindex, a, b, c, sollist):
	print('sieve')
	spos = [0] * (M + 1) # holds 0 and positive indices for sieving
	sneg = [0] * (M + 1) # holds 0 and negative indices for sieving
	for pindex in range(0, len(baselist)):
		p = baselist[pindex]
		sol0 = sollist[pindex][0]
		indpos = sol0
		indneg = - (sol0 - p)
		while indpos < M + 1:
			spos[indpos] = spos[indpos] + logs[pindex]
			indpos += p
		while ind < M + 1:
			sneg[indneg] = sneg[indneg] + logs[pindex]
			indneg += p
		if p != 2:
			sol1 = sollist[pindex][1]
			indpos = sol1
			indneg = sol1
			while ind < M + 1:
				if ind >= 0:
					spos[indpos] = spos[indpos] + logs[pindex]
					indpos += p
			while ind < M + 1:
				if ind >= 0:
					sneg[indneg] = sneg[indneg] + logs[pindex]
					indneg += p
	return checksieve(N, M, baselist, sqrts, logs, plusindex, minusindex, last, a, b, c, sollist, spos, sneg)

def checksieve(N, M, baselist, sqrts, logs, plusindex, minusindex, last, a, b, c, sollist, spos, sneg):
	print('checksieve')
	check = log((M * sqrt(N)), 2).n() - 1
	# smoothindices = []
	Q = []
	R = []
	x = ZZ['x'].gen()
	g = a**2 * x**2 + 2 * b * x + c 
	for indpos in range(0, len(spos)):
		if spos[indpos] >= check:
			gx = g(indpos)
			smooth = bsmoothcheck(gx, baselist)
			if smooth == 1:
				# smoothindices.append(indpos)
				R.append([[len(R)], basefactorize(gx, baselist)])
				Q.append(ind, gx)
	for indneg in range(1, len(sneg)):
		if sneg[indneg] >= check:
			gx = g(indneg)
			smooth = bsmoothcheck(sneg[indpos], baselist)
			if smooth == 1:
				# smoothindices.append(- indneg)
				R.append([[len(R)], basefactorize(gx, baselist)])
				Q.append(- ind, gx)
	if len(Q) > len(baselist):
		dep = lindep(N, Q, R) # factorization or None if we cannot yet find a factorization
		if dep != None:
			return dep # the factorization
	else:
		return newpoly(N, M, baselist, sqrts, logs, plusindex, minusindex, last) # didn't work, try new polynomial

def bsmoothcheck(t, q): # returns 1 if t is bsmooth over base and 0 otherwise
	print('bsmoothcheck')
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
	print('basefactorize')
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
	print('collectbase')
	b = floor(e**(sqrt(ln(D) * ln(ln(D))) / 2)) # find bound for base primes
	if b < 31: # never have too small a base
		b = 31
	if b > f:
		b = f
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


def lindep(N, Q, R): # find dependency, return either factorization or None
	print('lindep')
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
	print('checknumbers')
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
	print('factorize')
	r0 = x2 + y
	f0 = gcd(r0, N)
	return f0, int(N / f0)

def chooseq(N, M, largestbase):
	print('chooseq')
	close = int(sqrt((sqrt(2 * N)) / M)) # find prime q close to this number which is a quadratic residue mod N
	print('close : ' + str(close))
	pindex = 0
	p = primelist[pindex]
	while p <= largestbase: # q should not be in base
		pindex += 1
		p = primelist[pindex]
	pindex += 1
	p = primelist[pindex]
	primelen = len(primelist)
	while (p < close) and (p < primelen - 2): # if p < the goal and there are remaing primes available
		pindex += 1
		p = primelist[pindex] # use the next prime
	if p < primelen - 2: # we did not run out of primes
		lowdiff = close - primelist[pindex - 1]
		highdiff = p - close
		if highdiff <= lowdiff: # see if lower or higer prime is closer to goal
			qindex = pindex
			q = primelist[qindex]
		else:
			qindex = pindex - 1
			q = primelist[qindex]
	else:
		qindex = primelen - 1
		q = primelist[qindex]
	return qindex, q


print(MPQS(1050703, 100, 10))
# print(MPQS(15347))
# print(MPQS(15440))
# print(MPQS(64775585))
# print(MPQS(539873))
# print(MPQS(33569887))
# print(MPQS(6558422578784))