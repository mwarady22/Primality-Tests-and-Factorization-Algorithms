from sage.all import *
import sys
from PrimesList import *
from PrimitiveRoot import *
from PoweringAlgorithms import *
UCF = UniversalCyclotomicField()

# contains e(5040) (can use for numbers with 104 digits)
# and e(720720) (can be used for numbers with 474 digits)
etable = [15321986788854443284662612735663611380010431225771200, 2599265289938045790285087430718636500803510968991503583932015415681339876886274198699629634541272594255929473856752109528189646768583309952911477951292254958722748988899531120282182905307875781283119196897221941548215066799861837732382400]

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


def etablecalc(tablebound): # if you want a full table from t = 1 to t = tablebound
	table = []
	index = 0
	while index <= tablebound:
		table.append(ecalc(index))
	return table


def Precomps(B):
	et = findt(B)
	if et == etable[0]:
		t = 5040
	else:
		t = 720720

	h = 1
	q =  primelist[h]
	table
	while q <= et:
		if (et % q == 0): # finding primes q >= 3 that divide et
			g, ftable = Jacobi1(q) #
			k, chilist = Jacobi2(p, q) #
			Jpq, j, J3, J2 =  Jacobi3(p, k, chilist, ftable)
			table.append((q, g, ftable, Jpq, j, J3, J2))
		h += 1
		q = primelist[h]


def findt(B):
	for item in etable:
		if item > B:
			return item
		elif item**2 > B:
			return item
	return 'table contains no item large enough for use with B'


def Jacobi1(q): #
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


def Jacobi2(p, q): #
	k = vpcalc(p, q - 1)
	print(k)
	pk = LRbin(p, k, calculate_e(k))
	chilist = []
	for x in range(1, q - 1):
		chilist.append(UCF.gen(pk, x))
	return k, chilist

def Jacobi3(p, k, chilist, ftable): # # figure out return
	if (p >= 3) or ((p == 2) and (k == 3)):
		Jpq = 0
		for x in range(1, q - 1):
			Jpq += chilist[x - 1] * chilist[ftable[x - 1]]
		return Jpq, None, None, None
	elif (p == 2) and (k >= 3):
		
		Jpq = 0
		for x in range(1, q - 1):
			Jpq += chilist[x - 1] * chilist[ftable[x - 1]]
		
		j = 0
		for x in range(1, q - 1):
			fx = ftable[x - 1]
			zeta = UCF.gen(LRbin(2, k, calculate_e(k)), (2 * x + fx))
			# summand = LRbin(zeta, 2 * x, calculate_e(2 * x)) * LRbin(zeta, fx, calculate_e(fx))
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

# pages 463 - 464