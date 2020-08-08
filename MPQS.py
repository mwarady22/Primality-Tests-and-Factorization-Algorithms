from sage.all import *
import sys
from PrimesList import *
from Legendre import *
from SqrtModp import *

def MPQS(N, f, M):
	baselist, k = collectbase(N, f) # list of base elements, product of all base elements
	sqrts = [None] # will hold the modular square roots of each odd prime in the baselist, with None at index 0
	logs = [1]# will hold the base 2 log of each prime in the baselist
	mid = sqrt((sqrt(2 * N) / M))
	qind = None
	q = None
	counter = 0
	for p in baselist:
		sqrts.append(SqrtModp(N, p))
		logs.append(log(p, 2).n())
		if p <= mid:
			q = p
			qind = counter
		counter += 1
	indexlist = [q]
	plus = q + 1
	minus = q - 1
	while (plus < len(sqrts)) and (minus >= 0):
		indexlist.append(plus)
		plus += 1
		indexlist.append(minus)
		minus -= 1
	while (plus < len(sqrts)):
		indexlist.append(plus)
		plus += 1
	while (minus >= 0):
		indexlist.append(minus)
		minus -= 1

	return newpoly() #

def newpoly(N, M, baselist, sqrts, logs, indexlist, indexlistindex):
	qind = indexlist[indexlistindex]
	q = baselist[qind]
	a = q**2
	ainv = inverse_mod(a, 2)
	b = mod(N, a).sqrt()
	sollist = []
	for pindex in range(0, len(baselist)):
		p = baselist[pindex]
		sqrtp = sqrts[pindex]
		ainv = inverse_mod(a, p)
		soln0 = (ainv * (sqrtp - b)) % p
		if p != 2:
			soln1 = (ainv * (- sqrtp - b)) % p
			sollist.append([soln0, soln1])
		else:
			sollist.append([soln0])
	return sieve()

def sieve(N, M, baselist, logs, sollist):
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
				spos[indpos] = spos[indpos] + logs[pindex]
				indpos += p
			while ind < M + 1:
				sneg[indneg] = sneg[indneg] + logs[pindex]
				indneg += p
	return checksieve()

def checksieve(N, M, spos, sneg):
	check = log((M * sqrt(N)), 2).n() - 1
	checkindices = []
	






def QS(N):
	x = ZZ['x'].gen()
	Qa = (x + ceil(sqrt(N)))**2 - N # create a polynomial to evaluate at different inputs
	l = [] # will hold Qa(z) at index z
	q = [] # will hold [z + ceil(sqrt(N)), Qa(z)] at index z
	ran = ceil(sqrt(N)) # the range over which Qa will be evaluated
	for z in range(0, ran):
		Qaz = Qa(z)
		l.append(int(Qaz))
		q.append([z + ceil(sqrt(N)), Qaz])
	m = l.copy() # save copy of l so that l can be modified without loss of information
	for p in baselist: # sieve
		w = Integers(p)['w'].gen() # create integer mod p ring in variable w
		Ya = (w + ceil(sqrt(N)))**2 - N # create the polynomial Qa mod p
		rootlist = Ya.roots() # find roots of this polynomial
		for root, mult in rootlist: # divide out each prime completely from root + n * p
			index = int(root)
			while index < len(l):
				while l[index] / p == int(l[index] / p): # while there is still a factor of p, remove it
					l[index] = l[index] / p
				index += p
	R = []
	Q = []
	for ind in range(0, len(l)): # find which indices hold numbers factorable in the prime base
		if l[ind] == 1:
			R.append([[len(R)], basefactorize(m[ind], baselist)])
			Q.append(q[ind])
	dep = lindep(N, Q, R) # factorization or None if we cannot yet find a factorization
	if dep != None:
		return dep # the factorization
	else:
		return 'Cannot be solved with Qudratic Sieve, try Multiple Polynomial Quadratic Sieve'


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



# print(QS(15347))
# print(QS(15440))
# print(QS(64775585))
# print(QS(539873))
# print(QS(33569887))
# print(QS(6558422578784))