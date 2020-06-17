from sage.all import *
import sys
from PrimesList import *
from EuclideanAlgorithm import *


def FirstStage(N, B):
	# find k
	x = 2
	y = x
	c = 0
	h = 0
	j = h
	return nextprime(N, B, k, x, y, c, h, j)

def nextprime(N, B, k, x, y, c, h, j):
	h += 1
	if h > k:
		g = gcd(x - 1, N)
		if g = 1:
			return 'the algorithm did not succeed in spliting N'
		else:
			i = j
			x = y
			return # 5
	else:
		q = PrimesList[i]
		q1 = q
		l = floor(B / q)
		return computepower() #

def computepower(): #
	while q1 <= l:
		q1 = q * q1
	x = 