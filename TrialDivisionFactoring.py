from sage.all import *
import sys
from PrimesList import *

t = [6, 4, 2, 4, 2, 4, 6, 2] # unknown why these numbers
k = lastindex
lastelem = primes[lastindex] # the last element of the list of primes
lastpmod = lastelem % 30

if lastpmod == 1:
	j = 0
elif lastpmod == 7:
	j = 1
elif lastpmod == 11:
	j = 2
elif lastpmod == 13:
	j = 3
elif lastpmod == 17:
	j = 4
elif lastpmod == 19:
	j = 5
elif lastpmod == 23:
	j = 6
elif lastpmod == 29:
	j = 7


def TrialFactor(N, B=lastelem): # Finds a factor of N of size up to B
	h, m, l, d, r = 0, 0, 0, 0, 0
	return initialize(N, k, j, B, h, m, l, d, r) # branches to step 1


def initialize(N, k, j, B, h, m, l, d, r): # step 1
	if N <= 5:
		if N == 1:
			return [1]
		elif N == 2:
			return [2]
		elif N == 3:
			return [3]
		elif N == 4:
			return [2, 2]
		elif N == 5:
			return [5]
	else:
		i = - 1
		m = - 1
		l = floor(sqrt(N))
		h = - 1
		return nextprime(N, k, j, B, h, m, l, d, r) # branches to step 2


def nextprime(N, k, j, B, h, m, l, d, r): # step 2
	m += 1
	if m > k:
		h = j - 1
		return nextdivisor(N, k, j, B, h, m, l, d, r) # branches to step 5
	else:
		d = primes[m]
		return trialdivide(N, k, j, B, h, m, l, d, r) # branches to step 3
	

def trialdivide(N, k, j, B, h, m, l, d, r): # step 3
	r = N % d
	if r == 0:
		return [d]
	return primecheck(N, k, j, B, h, m, l, d, r) # branches to step 4
	

def primecheck(N, k, j, B, h, m, l, d, r): # step 4
	if d >= l:
		if N > 1:
			return "prime"
	if h < 0:
		return nextprime(N, k, j, B, h, m, l, d, r) # branches to step 2
	return nextdivisor(N, k, j, B, h, m, l, d, r) # branches to step 5


def nextdivisor(N, k, j, B, h, m, l, d, r): # step 5
	h = (h + 1) % 8
	d = d + t[h]
	if d > B:
		return "remaining divisors are greater than " + B
	else:
		return trialdivide(N, k, j, B, h, m, l, d, r) # branches to step 3




def TrialCompleteFactorization(N, B=lastelem):
	factorlist = []
	factor = Trial_Factor(N, B)
	while not isinstance(factor, str):
		factorlist = factorlist + factor
		for f in factor:
			N = int(N / f)
		factor = Trial_Factor(N, B)
		if N == 1:
			break
	if factor == 'prime':
		factorlist.append(N)
		return factorlist
	elif isinstance(factor, str):
		notprime = factor[len(factor) - 1]
		return factorlist, notprime



# print(Trial_Complete_Factorization(2065762, 15))
# print(Trial_Factor(2065762, 15))