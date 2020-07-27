from sage.all import *
import sys
from PrimesList import *

def SL(N): # 1

	B = floor((L(N))**.5)
	plist = getplist(B)
	K = 1
	d = floor(ln(B))
	return Kinit(N, B, plist, K, D) # 2



def getplist(B): # precompute a list of primes up to L(N)**.5
	l = []
	c = 0
	while primelist[c] <= B:
		l.append(primelist[c])
		c += 1
	return l

def L(x):
	return e**(sqrt(ln(x) * ln(ln(x))))

def Kinit(N, B, plist, K, D): # 2
	if ((K * N % 4) == 3):
		D = - K * N
	else:
		D = - 4 * K * N
	return # 3

def chooseform(): # 3 #
	########
	x = fp
	c = 0
	h = 1
	return # 4

def nextprime(): # 4 #
	h += 1
	if h > k:
		K += 1
		return # 2
	else:
		q = plist[h - 1]
		q1 = q
		l = floor(B / q)
		return # 5

def comppower(): # 5 #
	while q1 <= l:
		q1 *= 1
	x = x**q1 #### power in class group ####
	c += 1
	if c < 20:
		return # 4
	else:
		return # 6

def checksuccess(): # 6 #
	d1 = 0
	while x # not ambiguous form and d1 < d: #
		x *= x
		d1 += 1
	if x # not ambiguous form: #
		c = 0
		return # 4
	else:
		return # 7

def checkfinished(): # 7 #
	####
	# find factorization of KN corresponding to x
	# if factorization does not split N:
		return # 3
	# else:
		return # nontrivial factor of N