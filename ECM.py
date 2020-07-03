from sage.all import *
import sys
import random
from EuclideanAlgorithm import *
sys.setrecursionlimit(10**4)


def ECM(N, dmax=128): # N is composite, dmax (optional) is the maximum search depth
	g = gcd(N, 6) # check that 2, 3 are not factors
	if g > 1:
		return g, int(N / g)
	r = N
	power = 2
	while r >= 2: # check that N is not a perfect power
		r = N**(1 / power)
		if r == floor(r):
			return r, int(N / r)
		power += 1
	return choosecurve(N, dmax)

def choosecurve(N, dmax): # chose a random curve to use for computations
	b = randint(0, N - 1)
	x1 = randint(0, N - 1) # x coordinate of initial point
	y1 = randint(0, N - 1) # y coordinate of initial point
	return initpointcurve(N, dmax, b, x1, y1)

def initpointcurve(N, dmax, b, x1, y1): # initialize point and curve
	c = (y1**2 - x1**3 - b * x1) % N # find constant term of curve such that P = (x1, y1) is on the curve
	ZNZ = Integers(N)
	E = EllipticCurve(ZNZ, [b, c]) # create curve y**2 = x**3 + b * x + c (mod N)
	P = E(x1, y1) # create initial point
	return loopd(N, dmax, b, x1, y1, P, 2)

def loopd(N, dmax, b, x1, y1, P, d): # perform multiplications of P until one fails or the maximum search depth is reached
	try:
		Q = d * P
	except ZeroDivisionError as z:
		noinverse = z.args[0].split(' ')[2] # take number that failed to have inverse from the error message
		return failure(N, dmax, noinverse)
	P = Q
	d += 1
	if d <= dmax:
		return loopd(N, dmax, b, x1, y1, P, d)
	else: # give up on this curve, choose another and try again
		return choosecurve(N, dmax)

def failure(N, dmax, noinverse): # multiplication failed, have found a factor of N
	g = gcd(N, noinverse)
	if g < N: # have found a nontrivial factor of N, return
		return g, int(N / g)
	else: # factor of N is N, try again with new curve
		return choosecurve(N, dmax) #

# num = 492876847 * 982451653
# print(ECM(num, 128))