from sage.all import *
import sys
from EuclideanAlgorithm import *

def Pollard(N): # N is a composite integer, will find a non-trivial factor of N
	return initialize(N)

def initialize(N):
	return accumulate(N, 2, 2, 2, 1, 1, 1, 0)

def accumulate(N, y, x, x1, k, l, P, c):
	g = 1
	x = (x**2 + 1) % N
	P = (P * (x1 - x)) % N
	c += 1
	if c == 20:
		g = gcd(P, N)
		if g > 1:
			return backtrack(N, y, g, x1) # have found factor, maybe N
		else:
			y = x
			c = 0
	return advance(N, y, x, x1, k, l, P, c, g)

def advance(N, y, x, x1, k, l, P, c, g):
	k = k - 1
	if k != 0:
		return accumulate(N, y, x, x1, k, l, P, c)
	else:
		g = gcd(P, N)
		if g > 1:
			return backtrack(N, y, g, x1) # have found factor, maybe N
		else:
			x1 = x
			k = l
			l = 2 * l
			for j in range(0, k):
				x = (x**2 + 1) % N
			y = x
			c = 0
			return accumulate(N, y, x, x1, k, l, P, c)

def backtrack(N, y, g, x1):
	while not (g > 1):
		y = (y**2 + 1) % N
		g = gcd(x1 - y, N)
	if g < N:
		return g
	else:
		return 'algorithm fails'