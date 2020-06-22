from sage.all import *
import sys
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
	return table


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






# print(ecalc(5040))
# print(etable)
# print(Jacobi1(101))
# print(Jacobi2(17, 35))
# print(Precomps(7921))

# pages 463 - 464