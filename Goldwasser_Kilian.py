from sage.all import *
import sys
import random
from TrialDivisionFactoring import *
from Rabin_Miller import *
from Legendre import *
from SqrtModp import *


def GKTest(N): #
	print('GKTest')
	h = 0
	Ni = N
	return Nismall(Ni, h) #


# def Nismall(Ni, h): # 2
# 	print('Nismall')
# 	if Ni < 2**30:
# 		t = TrialFactor(Ni, 2**15)
# 		if not isinstance(t, str): # then it is not the case that Ni is prime or has only factors larger than 2**15
# 			return backtrack(Ni, h)
# 	return choosecurve(Ni, h) #

def Nismall(Ni, h): # 2
	print('Nismall')
	if Ni < 2**30:
		t = TrialFactor(Ni, 2**15)
		if not isinstance(t, str): # then it is not the case that Ni is prime or has only factors larger than 2**15
			return backtrack(Ni, h)
		return 'prime' #
	return choosecurve(Ni, h) #

def choosecurve(Ni, h): # 3
	print('choosecurve')
	a = randint(0, Ni - 1)
	b = randint(0, Ni - 1)
	z = 4 * a**3 + 27 * b**2
	if z % Ni == 0:
		a, b = choosecurve(Ni, h) #
	return schoof(Ni, h, a, b) #


def schoof(Ni, h, a, b): # 4 a, b, Ni
	print('schoof')
	ZNiZ = Integers(Ni)
	E = EllipticCurve(ZNiZ, [a, b])
	try:
		m = E.cardinality()
	except AttributeError:
		return backtrack(Ni, h)
	return mcheck(Ni, h, E, a, b, m) #


def mcheck(Ni, h, E, a, b, m): # 5
	print('mcheck')
	if (m % 2 == 1):
		return choosecurve(Ni, h) #
	else:
		q = m / 2
	print('q : ' + str(q))
	rabinmiller = RM(q)
	if rabinmiller == 'composite':
		return choosecurve(Ni, h) #
	else:
		return findp(Ni, h, E, a, b, m, q) #


def findp(Ni, h, E, a, b, m, q): # 6
	print('findp')
	while True:
		x = randint(0, Ni - 1)
		K = Kronecker((x**3 + a * x + b), Ni)
		if (K == 0) or (K == 1):
			break
	y2 = x**3 + a * x + b
	y = SqrtModp(y2, Ni)
	if isinstance(y, str):
		return backtrack(Ni, h) #
	else:
		return checkp(Ni, h, E, a, b, m, q, x, y) #


def checkp(Ni, h, E, a, b, m, q, x, y): # 7
	print('checkp')
	P = E(x, y)
	P1 = m * P
	OE = E(0)
	try:
		z = int(m / q)
		P2 = z * P
	except ZeroDivisionError:
		return backtrack(Ni, h) #

	if P1 != OE:
		return backtrack(Ni, h) #
	elif P2 == OE:
		return findp(Ni, h, E, a, b, m, q) #
	else:
		return recurse(Ni, h, q) #


def recurse(Ni, h, q): # 8
	print('recurse')
	h += 1
	Ni = q
	return Nismall(Ni, h) #


def backtrack(Ni, h): # 9
	print('backtrack')
	if h == 0:
		return 'composite'
	else:
		h -= 1
		return choosecurve(Ni, h) #

print(GKTest(11037271757))