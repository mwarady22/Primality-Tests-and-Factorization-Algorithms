from sage.all import *
import sys
from EuclideanAlgorithm import *
from PrimesList import *
from PrimitiveRoot import *
from PoweringAlgorithms import *
UCF = UniversalCyclotomicField()

def ecalc(t):
	if (t % 2 == 1):
		return 2
	ret = 2
	h = 0
	q = primelist[h]
	while q - 1 <= t:
		if t % (q - 1) == 0:
			v = vpcalc(q, t)
			ret *= q**(v + 1)
		h += 1
		q = primelist[h]
	return ret

def vpcalc(q, t): # v_q(t) denotes the highest exponent k such that q**k | t
	v = 0
	while (t % (q**v)) == 0:
		v += 1
	v -= 1
	return v


def etablecalclarge(tablebound): # if you want a full table from t = 1 to t = tablebound >= 1000
	table = []
	index = 0
	while index <= tablebound:
		table.append(ecalc(index))
		index += 1000
	return table

def etablecalcsmall(tablebound): # if you want a full table from t = 1 to t = tablebound < 1000
	table = []
	index = 0
	while index <= tablebound:
		table.append(ecalc(index))
		index += 1
	return table


def Precomps(B): # use for B less than 1000000

	if B < 1000:
		tab = etablecalcsmall(B)
		ind = 0
		while tab[ind]**2 <= B:
			ind += 1
		et = tab[ind]
		t = ind * 1000
	else:
		tab = etablecalclarge(B)
		ind = 0
		while tab[ind]**2 <= B:
			ind += 1
		et = tab[ind]
		t = ind * 1000

	h = 1
	q =  primelist[h]
	table = []
	while q <= floor(sqrt(et)):
		if (et % q == 0): # finding primes q >= 3 that divide et
			g, ftable = Jacobi1(q) #
			l = 1
			p = primelist[l]
			while p <= (q - 1):
				if ((q - 1) % p == 0):
					k, chilist = Jacobi2(p, q)
					Jpq, j, J3, J2 =  Jacobi3(p, q, k, chilist, ftable)
					table.append((q, g, p, ftable, Jpq, j, J3, J2))
				l += 1
				p = primelist[l]
		h += 1
		q = primelist[h]
	return t, et, table


def findt(B):
	for item in etable:
		if item > B:
			return item
		elif item**2 > B:
			return item
	return 'table contains no item large enough for use with B'


def Jacobi1(q):
	g = PrimRoot(q)
	ftable = []
	for x in range(1, q - 1):
		gf = (1 - PowerModm(g, x, q)) % q
		for y in range(1, q - 1):
			gfprime = PowerModm(g, y, q)
			if gf == gfprime:
				ftable.append(y)
				break
	return g, ftable


def Jacobi2(p, q):
	k = vpcalc(p, q - 1)
	pk = LRbin(p, k, calculate_e(k))
	chilist = []
	for x in range(1, q - 1):
		chilist.append(UCF.gen(pk, (x % pk)))
	return k, chilist

def Jacobi3(p, q, k, chilist, ftable):
	if (p >= 3) or ((p == 2) and (k == 3)):
		Jpq = 0
		for x in range(1, q - 1):
			Jpq += chilist[x - 1] * chilist[ftable[x - 1] - 1]
		return Jpq, None, None, None
	elif (p == 2) and (k >= 3):
		Jpq = 0
		for x in range(1, q - 1):
			Jpq += chilist[x - 1] * chilist[ftable[x - 1] - 1]
		
		j = 0
		for x in range(1, q - 1):
			fx = ftable[x - 1]
			zeta = UCF.gen(LRbin(2, k, calculate_e(k)), (2 * x + fx))
			j += zeta

		J3 = Jpq * j

		J2 = 0
		for x in range(1, q - 1):
			J2 += UCF.gen(8, (3 * x + ftable[x - 1]))
		J2 = J2 * J2

		return Jpq, j, J3, J2


def JacobiTest(N, t, et, table): # N <= B is a strong pseudo-prime in 20 randomly chosen bases
	# 1 gcd
	gcd = checkgcd(N, t, et)
	if gcd > 1:
		return 'composite'

	# 2 init
	# 3 loop on char
	# 4a check *b p >= 3
	# 4b check *b p = 2 and k >= 3
	# 4c check *b p = 2 and k = 2
	# 4d check *b p = 2 and k = 1
	# 5 check Lp
	# 6 induction

def checkgcd(N, t, et):
	return max(gcd(t, N), gcd(et, N))

def initJacobi(N, t): #
	l_p = [0, 0] # for p the jth prime dividing t insert 1 if p >= 3 and N**(p - 1) != 1 (mod p**2) else 0
	h = 2
	curprime = primelist[h]
	while curprime <= floor(sqrt(t)):
		print(curprime)
		if (t % curprime == 0) and (PowerModm(N, curprime - 1, curprime**2) != 1):
			l_p.append([curprime, 1])
		else:
			l_p.append([curprime, 0])
		h += 1
		curprime = primelist[h]
	return l_p

def chiloop(): #
	# figure out p, q, k
	if (p >= 3): # case 1
		# 4a
	elif (p == 2) and (k >= 3): # case 2
		# 4b
	elif (p == 2) and (k == 2): # case 3
		# 4c
	else: # (p == 2) and (k == 1) # case 4
		# 4d

def case1(p, k): #
	E = []
	pk = LRbin(p, k, calculate_e(k))
	r = N % pk
	for n in range(0, pk):
		if (n % p != 0):
			E.append(n)
	theta = 0
	
	####



def case2(k): # p = 2

def case3(k): # p = 2

def case4(): # p = 2, k = 1
	




# print(ecalc(5040))
# print(etable)
# print(Jacobi1(101))
# print(Jacobi2(17, 35))
# print(Precomps(7921))
# print(initJacobi(1000, 1700))
# print(JacobiTest(N, Precomps(B)))

# pages 463 - 464
# sigma 307 article