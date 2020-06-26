from sage.all import *
import sys
import random
from TrialDivisionFactoring import *
from Rabin_Miller import *
from Legendre import *
from SqrtModp import *

def GKTest(N): #
	h = 0
	Ni = N
	return # 2

def Nismall(): # 2
	if Ni < 2**30:
		t = TrialFactor(Ni, 2**15)
		if not isinstance(t, str):
			return # 9
	return # 3

def choosecurve(): # 3
	a = randint(0, Ni - 1)
	b = randint(0, Ni - 1)
	z = 4 * a**3 + 27 * b**2
	if z % Ni == 0:
		a, b = choosecurve() #
	return a, b # 4




def schoof(a, b, Ni): # 4 a, b, Ni
	ZNiZ = Integers(Ni)
	E = EllipticCurve(ZNiZ, [a, b])
	m = E.cardinality()
	return # 5



def mcheck(): # 5
	if (m % 2 == 1):
		return # 3
	else:
		q = m / 2

	rabinmiller = RM(q)
	if rabinmiller == 'composite':
		return # 3
	else:
		return # 6


def findp(): # 6
	while True:
		x = randint(0, Ni - 1)
		K = Kronecker((x**3 + a * x + b), Ni)
		if (K == 0) or (K == 1):
			break
	y2 = x**3 + a * x + b
	y = SqrtModp(y2, Ni)
	if isinstance(y, str):
		return # 9
	else:
		return # 7


def checkp(): # 7
	P = E(x, y)
	P1 = m * P
	P2 = (m / q) * P
	# http://sage-doc.sis.uta.fi/reference/curves/sage/schemes/elliptic_curves/ell_point.html
	# 'arithmetic over'











	


def recurse(): # 8
	h += 1
	Ni = q
	return Nismall() #

def backtrack(): # 9
	if h == 0:
		return 'composite'
	else:
		h -= 1
		return # 3

# schoof(4, 5, 6)