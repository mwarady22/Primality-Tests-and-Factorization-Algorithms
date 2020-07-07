from sage.all import *
import sys
from PrimesList import *
from Legendre import *
from EuclideanAlgorithm import *

def CFRAC(N): #
	print('CFRAC')
	k = 1
	a = 2 * floor(sqrt(N))
	y = floor(sqrt(N))
	z = 1
	Epair = (1, 0)
	Fpair = (0, 1)
	return kloop(N, k, a, y, z, Epair, Fpair)

def nextk(k): # calculates next smallest squarefree multiplier k
	k += 1
	while moebius(k) == 0:
			k += 1
	return k

def calclooplen(D):
	F = QuadraticField(D)
	sqrtDformatted = F(sqrt(D))
	cf = continued_fraction(sqrtDformatted)
	preperiod = cf.preperiod()
	period = cf.period()
	return len(preperiod) + len(period)

def kloop(N, k, a, y, z, Epair, Fpair): # 
	print('kloop')
	D = k * N
	baselist, q = collectbase(D)
	print('baselist')
	print(baselist)
	l = len(baselist)
	looplen = calclooplen(N)
	counter = 0
	Q = [] # will store the pairs [Ak, abs(Ak**2)] such that Ak**2 can be factored over baselist (all mod N)
	R = [] # will store the exponents of the factorization of each Ak2 in Q
	while (len(Q) < l + 1) and (counter < looplen): # need l + 1 equations to ensure linear dependence
		print('make Q : ' + str(counter))
		rootN = sqrt(N)
		a, y, z, Epair, Fpair, Ak, Bk = nextfraction(N, a, y, z, Epair, Fpair)
		Ak2 = (Ak**2) % N
		if Ak2 > 2 * rootN:
			Ak2 -= N
		elif Ak2 < - 2 * rootN:
			Ak2 += N
		if bsmoothcheck(Ak2, q) == 1:
			explist = basefactorize(Ak2, baselist)
			zerovector = True
			for r in range(0, len(explist)):
				explist[r] = explist[r] % 2
				if explist[r] != 0:
					zerovector = False
			alreadyin = False
			for item in R:
				if explist == item:
					alreadyin = True
					break
			if (alreadyin == False) and (zerovector == False):
				Q.append([Ak, abs(Ak2)])
				print('added : ' + str([Ak, abs(Ak2)]) + '****************')
				R.append(explist)
			else:
				print('already in : ' + str([Ak, abs(Ak2)]) + '****************')

		else:
			print('did not add : ' + str([Ak, abs(Ak2)]))
		counter += 1
	if not len(Q) == l + 1:
		k = nextk(k)
		print('k : ' + str(k))
		return kloop(N, k, (2 * floor(sqrt(N))), floor(sqrt(N)), 1, (1, 0), (0, 1))
	print('Q')
	print(Q)
	xlist = lindep(R, l)
	x2 = 1
	y2 = 1
	for w in xlist:
		print('iterate through xlist')
		Ak, Ak2 = Q[w]
		print('Ak, Ak2 : ' + str(Q[w]))
		x2 = (x2 * Ak) % N
		y2 = (y2 * Ak2)
	y = sqrt(y2)
	print('y2 : ' + str(y2))
	print('y : ' + str(y))
	print('x2 : ' + str(x2))
	print('k : ' + str(k))
	return factorize(N, x2, y)

	
def collectbase(N): #
	print('collectbase')
	b = floor(e**(sqrt(ln(N) * ln(ln(N))) / 2)) # find bound for base primes
	print('b')
	print(b)
	plist = [2] # collect base primes N <= b such that Kroencker(N, p) = 1
	pindex = 1 # index 1 will start the loop at the second prime, 3
	p = primelist[pindex]
	q = 2 # the product of the primes in the base
	while p <= b:
		print('build plist')
		if Kronecker(N, p) == 1:
			plist.append(p)
			q *= p
		pindex += 1
		p = primelist[pindex]
	return plist, q #

def nextfraction(N, a, y, z, Epair, Fpair): #
	print('nextfraction')
	yk = a * z - y
	zk = floor((N - yk**2) / z)
	ak = floor((floor(sqrt(N)) + yk) / zk)
	Epair = (Epair[1], (ak * Epair[1] + Epair[0]) % N)
	Fpair = (Fpair[1], (ak * Fpair[1] + Fpair[0]) % N)
	Ak = (Epair[0] + floor(sqrt(N)) * Fpair[0]) % N
	Bk = Fpair[0]
	return ak, yk, zk, Epair, Fpair, Ak, Bk #

def bsmoothcheck(t, q): # returns 1 if t is bsmooth over base and 0 otherwise
	print('bsmoothcheck')
	while True:
		g = gcd(t, q)
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
	print('R')
	print(R)
	M = matrix(Integers(2), R).transpose() #
	print(M)
	redM = M.echelon_form()
	print('\n')
	print(redM)
	v = matrix([0] * l).transpose()
	print(v)
	xlist = [] # stores which exponent vectors are linearly dependent
	for x in range(0, l + 1):
		if redM[:,x] != v:
			xlist.append(x)
			print('added col ' + str(x) + ' to xlist')
	print(xlist)
	return xlist

def factorize(N, x2, y): # factor as x2 - y, x2 + y gcd with N
	print('factorize')
	r0 = x2 + y
	print('r0 : ' + str(r0))
	r1 = x2 - y
	print('r1 : ' + str(r1))
	g = gcd(r0, r1)
	r0 = int(r0 / g)
	print('r0 : ' + str(r0))
	f0 = gcd(r0, N)
	f1 = gcd(r1, N)

	return f0, f1

# N = 10403 # 101 * 103 period 4
# N = 3053
# N = 3599 # period 2
# N = 3054
# N = 3052
# N = 2257
N = 9271 #
# N = 418679 #
# N = 6426786497 #
# N = 665059573553159 #
# N = 729461128210276840421 #
# N = 194551432662383450877400470563 #
c = CFRAC(N)
print('factorization : ' + str(c))
print(str(c[0]) + '*' + str(c[1]) + '=' + str(c[0] * c[1]))
print('N : ' + str(N))


# D = kN
