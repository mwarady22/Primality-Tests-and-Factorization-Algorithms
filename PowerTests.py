from sage.all import *
import sys

def IntSqrt(N): # finds the floor of the square root of N
	if (N < 0): # if negative, return the answer of the positive version multiplied by i
		return (loop(- N, - N))*i
	else:
		return loop(N, N)

def loop(N, x):
	y = floor((x + floor(N / x)) / 2)
	if (y < x):
		return loop(N, y)
	else:
		return x

def precomp(): # only needed once
	q11 = [0] * 11
	for k in range(0, 6):
		q11[(k**2 % 11)] = 1
	q63 = [0] * 63
	for l in  range(0, 32):
		q63[(l**2 % 63)] = 1
	q64 = [0] * 64
	for m in  range(0, 32):
		q64[(m**2 % 64)] = 1
	q65 = [0] * 65
	for p in  range(0, 33):
		q65[(p**2 % 65)] = 1
	return q11, q63, q64, q65

def SquareTest(N, q11, q63, q64, q65): # tests if N is a perfect square
	neg = False
	if N < 0: # if N is negative, return the answer of the positive versin multiplied by i
		neg = True
		N = - N
	t = N % 64
	if q64[t] == 0:
		return 'not square'
	r = N % 45045
	if q63[r % 63] == 0:
		return 'not square'
	elif q65[r % 65] == 0:
		return 'not square'
	elif q11[r % 11] == 0:
		return 'not square'
	q = IntSqrt(N) # calculate the floor of the square root of N
	if q**2 != N:
		return 'not square'
	else: # in this case the floor of the square root of N is equal to N so that is the square root
		if neg:
			return q * i
		else:
			return q