from sage.all import *
import sys
import random
from EuclideanAlgorithm import *
sys.setrecursionlimit(10**4)


def ECM(N, dmax=128):
	print('ECM')
	g = gcd(N, 6)
	if g > 1:
		return g, int(N / g)
	r = N
	power = 2
	while r >= 2:
		print('while')
		r = N**(1 / power)
		if r == floor(r):
			return r, int(N / r)
		power += 1
	return choosecurve(N, dmax) #

def choosecurve(N, dmax): #
	print('choosecurve')
	b = randint(0, N - 1)
	x1 = randint(0, N - 1)
	y1 = randint(0, N - 1)
	return initpoint(N, dmax, b, x1, y1) #

def initpoint(N, dmax, b, x1, y1): #
	print('initpoint')
	c = (y1**2 - x1**3 - b * x1) % N
	ZNZ = Integers(N)
	E = EllipticCurve(ZNZ, [b, c])
	P = E(x1, y1)
	return loopd(N, dmax, b, x1, y1, P, 2) # d = 2, P

def loopd(N, dmax, b, x1, y1, P, d): # d #
	print('loopd : ' + str(d))
	try:
		print('try')
		Q = d * P
	except ZeroDivisionError as z:
		print('ZeroDivisionError')
		noinverse = z.args[0].split(' ')[2]
		print('noinverse')
		return failure(N, dmax, noinverse) #
	P = Q
	d += 1
	if d <= dmax:
		return loopd(N, dmax, b, x1, y1, P, d) #
	else:
		return choosecurve(N, dmax) #

def failure(N, dmax, noinverse): #
	print('failure')
	g = gcd(N, noinverse)
	if g < N:
		return g, int(N / g)
	else:
		return choosecurve(N, dmax) #

# num = 492876847 * 982451653
# print(ECM(num, 128))