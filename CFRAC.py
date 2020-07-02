from sage.all import *
import sys
from PrimesList import *
from Legendre import *
from EuclideanAlgorithm import *

def CFRAC(N): #
	print('CFRAC')
	rootN = sqrt(N)
	baselist, k = collectbase(N)
	l = len(baselist)
	a = 2 * floor(rootN)
	y = floor(rootN)
	z = 1
	E0 = 1
	E1 = 0
	F0 = 0
	F1 = 1
	counter = 0
	Q = [] # will store the pairs [Ak, abs(Ak**2)] such that Ak**2 can be factored over baselist (all mod N)
	while len(Q) < l + 1: # need l + 1 equations to ensure linear dependence
		print('make Q : ' + str(counter))
		Ak = None
		Bk = None
		# print('old ak, yk, zk, E2, E3, F2, F3, Ak, Bk')
		# print(a, y, z, E0, E1, F0, F1, Ak, Bk)
		a, y, z, E0, E1, F0, F1, Ak, Bk = nextfraction(N, a, y, z, E0, E1, F0, F1)
		# print('new ak, yk, zk, E2, E3, F2, F3, Ak, Bk')
		# print(a, y, z, E0, E1, F0, F1, Ak, Bk)
		Ak2 = (Ak**2) % N
		if Ak2 > 2 * rootN:
			Ak2 -= N
		elif Ak2 < - 2 * rootN:
			Ak2 += N
		if bsmoothcheck(Ak2, k) == 1:
			Q.append([Ak, abs(Ak2)])
		print('l + 1')
		print(l + 1)
		print('len(Q)')
		print(len(Q))
		counter += 1
	R = [] # will store the exponents of the factorization of each Ak2 in Q
	for pair in Q:
		print('add to R')
		Ak, Ak2 = pair
		explist = basefactorize(Ak2, baselist)
		R.append(explist)
	xlist = lindep(R, l)
	x2 = 1
	y2 = 1
	for w in xlist:
		print('iterate through xlist')
		Ak, Ak2 = Q[w]
		x2 = (x2 * Ak) % N
		y2 = (y2 * Ak2) % N
	# yfactors = basefactorize(y2)
	# y = 1 # sqareroot of y2
	# for f in range(0, len(basefactorize)):
	# 	exp = yfactors[f]
	# 	factor = baselist[f]
	# 	y *= factor**(exp / 2)
	y = sqrt(y2)

	# factor as x2 - y, x2 + y gcd with N
	r0 = x2 + y
	r1 = x2 - y
	f0 = gcd(r0, N)
	f1 = gcd(r1, N)
	return f0, f1

	
def collectbase(N): #
	print('collectbase')
	b = floor(e**(sqrt(ln(N) * ln(ln(N))) / 2)) # find bound for base primes
	plist = [2] # collect base primes N <= b such that Kroencker(N, p) = 1
	pindex = 1 # index 1 will start the loop at the second prime, 3
	p = primelist[pindex]
	k = 1 # the product of the primes in the base
	while p <= b:
		print('build plist')
		if Kronecker(N, p) == 1:
			plist.append(p)
			k *= p
		pindex += 1
		p = primelist[pindex]
	return plist, k #

def nextfraction(N, a, y, z, E0, E1, F0, F1): #
	print('nextfraction')
	yk = a * z - y
	zk = floor((N - yk**2) / z)
	ak = floor((floor(sqrt(N)) + yk) / zk)
	E2 = E1
	E3 = (ak * E1 + E0) % N
	F2 = F1
	F3 = (ak * F1 + F0) % N
	Ak = (E1 + floor(sqrt(N)) * F1) % N
	Bk = F1
	return ak, yk, zk, E2, E3, F2, F3, Ak, Bk #

def bsmoothcheck(t, k): # returns 1 if t is bsmooth over base and 0 otherwise
	print('bsmoothcheck')
	while True:
		g = gcd(t, k)
		if g == 1:
			return 0
		else:
			while t % g == 0:
				t = t / g
			if abs(t) == 1:
				return 1

def basefactorize(t, baselist): #
	print('basefactorize')
	l = len(baselist)
	explist = [0] * l
	for x in range(0, l):
		if t == 1:
			break
		while (t % baselist[x] == 0):
			t = t / baselist[x]
			explist[x] += 1
	return explist #

def lindep(R, l): #
	print('lindep')
	ZZmod2ZZ = ZZ(2)
	print('R')
	print(R)
	S = []
	for lst in R:
		for item in lst:
			S.append(item)
	M = matrix(ZZmod2ZZ, len(R[0]), S) #
	print(M)
	M = matrix(ZZmod2ZZ, R).transpose()
	redM = A.echelon_form()
	v = matrix([0] * l)
	xlist = [] # stores which exponent vectors are linearly dependent
	for x in range(0, l + 1):
		if redM[:,x] != v:
			xlist.append(x)
	return xlist

print(CFRAC(3053))