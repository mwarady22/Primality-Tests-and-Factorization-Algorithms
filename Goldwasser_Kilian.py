from sage.all import *
import sys
import random
from TrialDivisionFactoring import *
from Rabin_Miller import *
from Legendre import *
from SqrtModp import *


def GKTest(N):
	h = 0
	Ni = N
	return Nismall(Ni, h)

def Nismall(Ni, h):
	if Ni < 2**20: # if Ni is small enough just trial divide
		t = TrialFactor(Ni, 2**10)
		if not isinstance(t, str): # then it is not the case that Ni is prime or has only factors larger than 2**15
			return backtrack(Ni, h)
		return 'prime'
	return choosecurve(Ni, h)

def choosecurve(Ni, h): # find coefficients a, b for an elliptic curve y**2 = x**3 + a * x + b
	a = randint(0, Ni - 1)
	b = randint(0, Ni - 1)
	z = 4 * a**3 + 27 * b**2
	if z % Ni == 0: # if discriminant is a multiple of Ni, choose new curve
		a, b = choosecurve(Ni, h)
	return schoof(Ni, h, a, b)


def schoof(Ni, h, a, b): # find the number of integer points on the chosen curve
	ZNiZ = Integers(Ni)
	E = EllipticCurve(ZNiZ, [a, b]) # construct curve
	try:
		m = E.cardinality() # throws error if fails
	except AttributeError:
		return backtrack(Ni, h)
	return mcheck(Ni, h, E, a, b, m)


def mcheck(Ni, h, E, a, b, m):
	if (m % 2 == 1):
		return choosecurve(Ni, h)
	else:
		q = m / 2
	rabinmiller = RM(q) # check that q passes Rabin-Miller
	if rabinmiller == 'composite':
		return choosecurve(Ni, h)
	else:
		return findp(Ni, h, E, a, b, m, q)


def findp(Ni, h, E, a, b, m, q):
	while True: # choose random x until curve gives an integer point
		x = randint(0, Ni - 1)
		y2 = (x**3 + a * x + b) % Ni
		K = Kronecker(y2, Ni)
		if (K == 0) or (K == 1):
			break
	y = SqrtModp(y2, Ni)
	if isinstance(y, str):
		return backtrack(Ni, h)
	else:
		return checkp(Ni, h, E, a, b, m, q, x, y)


def checkp(Ni, h, E, a, b, m, q, x, y):
	P = E(x, y) # make point P
	P1 = m * P
	OE = E(0) # make point at infinity
	try:
		z = int(m / q)
		P2 = z * P
	except ZeroDivisionError:
		return backtrack(Ni, h)
	if P1 != OE: # check if P1 is not the identity element
		return backtrack(Ni, h)
	elif P2 == OE: # check if P2 is the identity element 
		return findp(Ni, h, E, a, b, m, q)
	else:
		return recurse(Ni, h, q)


def recurse(Ni, h, q):
	h += 1
	Ni = q
	return Nismall(Ni, h)


def backtrack(Ni, h):
	if h == 0:
		return 'composite'
	else:
		h -= 1
		return choosecurve(Ni, h)