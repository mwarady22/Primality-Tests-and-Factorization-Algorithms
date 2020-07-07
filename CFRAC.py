from sage.all import *
import sys
from PrimesList import *
from Legendre import *
from EuclideanAlgorithm import *

def CFRAC(N): #
	if N % 2 == 0:
		return 2, int(N / 2)
	s = sqrt(N - 2)
	if s == int(s):
		return 'Cannot use CFRAC on numbers of the form n**2 + 2.'
	sqrtN = sqrt(N)
	if floor(sqrtN)**2 == N:
		return sqrtN, sqrtN
	# cf = continued_fraction(sqrt(N))
	F = QuadraticField(N)
	sqrtNformatted = F(sqrt(N))
	cf = continued_fraction(sqrtNformatted)
	print(cf)
	period = cf.period()
	print(period)
	k = 1
	factor0, factor1 = body(N, k, len(period))
	while isinstance(factor0, str) or (factor0 == k) or (factor1 == k):
		# change k
		k += 1
		while moebius(k) == 0:
			k += 1
		D = k * N
		G = QuadraticField(D)
		sqrtDformatted = G(sqrt(D))
		cf1 = continued_fraction(sqrtDformatted)
		period = cf1.period()
		factor0, factor1 = body(D, k, len(period))
	if (factor0 % k == 0):
		factor0 = int(factor0 / k)
	else:
		factor1 = int(factor1 / k)
	return factor0, factor1

def body(N, k, length): #
	print('body')
	rootN = sqrt(N)
	baselist, q = collectbase(N) # either returns a list of primes and the product of all those primes or a tuple containing a factorization
	if not isinstance(baselist, list): # in this case a tuple containing a factorization was returned
		return baselist, q # return the factorization
	print('baselist')
	print(baselist)
	l = len(baselist)
	a = 2 * floor(rootN)
	y = floor(rootN)
	z = 1
	Epair = (1, 0)
	Fpair = (0, 1)
	counter = 0 #
	Q = [] # will store the pairs [Ak, abs(Ak**2)] such that Ak**2 can be factored over baselist (all mod N)
	R = [] # will store the exponents of the factorization of each Ak2 in Q
	while (len(Q) < l + 1) and counter < length: # need l + 1 equations to ensure linear dependence
		print('make Q : ' + str(counter))
		Ak = None
		Bk = None
		a, y, z, Epair, Fpair, Ak, Bk = nextfraction(N, a, y, z, Epair, Fpair)
		Ak2 = (Ak**2) % N
		if Ak2 > 2 * rootN:
			Ak2 -= N
		elif Ak2 < - 2 * rootN:
			Ak2 += N
		if bsmoothcheck(Ak2, q) == 1:
			Q.append([Ak, abs(Ak2)])
			explist = basefactorize(Ak2, baselist)
			evenexps = True
			print('evenexps for : ' + str(explist))
			for exp in explist:
				if exp % 2 == 1:
					print('odd : '+ str(exp))
					evenexps = False
					break
			if evenexps == True:
				x2 = Ak
				yy = sqrt(Ak2)
				ret = factor(N, x2, yy)
				if ret[0] != N and ret[1] != N:
					return ret
			R.append(explist)
			print('added : ' + str([Ak, abs(Ak2)]) + '****************')
		else:
			print('did not add : ' + str([Ak, abs(Ak2)]))
		counter += 1 #
	if counter == length:
		return 'period too small', None
	print('Q')
	print(Q)
	for pair in Q:
		print('add to R')
		Ak, Ak2 = pair
		
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
	return factor(N, x2, y)

	
def collectbase(N): #
	print('collectbase')
	b = floor(e**(sqrt(ln(N) * ln(ln(N))) / 2)) # find bound for base primes
	print('b : ' + str(b))
	plist = [2] # collect base primes N <= b such that Kroencker(N, p) = 1
	pindex = 1 # index 1 will start the loop at the second prime, 3
	p = primelist[pindex]
	q = 2 # the product of the primes in the base
	while p <= b:
		print('build plist')
		K = Kronecker(N, p)
		if K == 1:
			plist.append(p)
			q *= p
		elif K == 0:
			return (p, int(N / p))
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
	xlist = [] # stores which exponent vectors are linearly dependent
	for x in range(0, l + 1):
		if redM[:,x] != v:
			xlist.append(x)
			print('added col ' + str(x) + ' to xlist')
	print(xlist)
	return xlist

def factor(N, x2, y): # factor as x2 - y, x2 + y gcd with N
	print('x2 : ' + str(x2))
	print('y : ' + str(y))
	r0 = x2 + y
	print('r0 : ' + str(r0))
	r1 = x2 - y
	print('r1 : ' + str(r1))
	f0 = gcd(r0, N)
	f1 = gcd(r1, N)
	return f0, f1

# print(CFRAC(10403)) # 101 * 103
# print(CFRAC(3053))
# print(CFRAC(3599)) #
# print(CFRAC(3054))
# print(CFRAC(3052))
# print(CFRAC(2257))
# print(CFRAC(9271)) #
print(CFRAC(418679)) #
# print(CFRAC(6426786497)) #
# print(CFRAC(665059573553159)) #
# print(CFRAC(729461128210276840421)) #
# print(CFRAC(194551432662383450877400470563)) #