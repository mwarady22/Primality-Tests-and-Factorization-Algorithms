from sage.all import *
import sys
from PrimesList import *
from EuclideanAlgorithm import *
from PoweringAlgorithms import *


def FirstStage(N, B):
	z = 0
	print(B)
	while primelist[z] < B:
		print(primelist[z])
		z += 1
	k = z - 1
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
		if g == 1:
			return 'the algorithm did not succeed in spliting N'
		else:
			h = j
			x = y
			return backtrack(N, B, k, x, y, c, h, j) #
	else:
		q = primelist[h]
		q1 = q
		l = floor(B / q)
		return computepower(N, B, k, x, y, c, h, j, l, q, q1) #

def computepower(N, B, k, x, y, c, h, j, l, q, q1): #
	while q1 <= l:
		q1 = q * q1
	x = PowerModm(x, q1, N)
	c += 1
	if c < 20:
		return nextprime(N, B, k, x, y, c, h, j)
	else:
		return computegcd(N, B, k, x, y, c, h, j) #

def computegcd(N, B, k, x, y, c, h, j): #
	g = gcd(x - 1, N)
	if g == 1:
		c = 0
		j = h
		y = x
		return nextprime(N, B, k, x, y, c, h, j)
	else:
		h = j
		x = y
		return backtrack(N, B, k, x, y, c, h, j) #

def backtrack(N, B, k, x, y, c, h, j): #
	h += 1
	q = primelist[h]
	q1 = q
	return finshedcheck(N, B, k, x, y, c, h, j) #

def finishedcheck(N, B, k, x, y, c, h, j): #
	x = PowerModm(x, q1, N)
	g = gcd(x - 1, N)
	if g == 1:
		q1 = q * q1
		if q1 <= B:
			return finishedcheck(N, B, k, x, y, c, h, j) #
		else:
			return backtrack(N, B, k, x, y, c, h, j) #
	else:
		if g < N:
			return g
		else:
			return 'the algorithm has failed'

# print(FirstStage(191657, 10000))