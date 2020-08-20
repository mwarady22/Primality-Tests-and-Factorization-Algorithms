from sage.all import *
import sys
import random
from Legendre import *
from PoweringAlgorithms import *

def twokcalc(r): # twokcalc(r) denotes the highest exponent k such that 2**k | r
	v = 0
	while (r % (2**v)) == 0:
		v += 1
	v -= 1
	return v

def SqrtModp(a, p): # finds square root mod p for odd prime p
	a = a % p # reduce a so as to make calculations faster
	k = twokcalc(p - 1)
	q = int((p - 1) / (2**k)) # p - 1 = 2**k * q
	while True: # find n such that (n / p) = - 1
		n = randint(0, p - 1)
		kr = Kronecker(n, p)
		if kr == - 1:
			break
	z = PowerModm(n, q, p)
	y = z
	r = k
	x = PowerModm(a, int((q - 1) / 2), p)
	b = (a * x**2) % p
	x = (a * x) % p
	return findexp(b, p, x, r, a, y)

def findexp(b, p, x, r, a, y):
	if (b % p) == 1:
		return x # found square root
	m = 1
	twom = 2
	while True: # find smallest b**(2**m) % p = 1
		b2m = PowerModm(b, twom, p)
		if b2m == 1:
			break
		m += 1
		twom *= 2
	if m == r:
		return str(a) + ' is not a quadratic residue mod ' + str(p)
	return redexp(m, y, b, p, x, r, a)

def redexp(m, y, b, p, x, r, a):
	expexp = r - m - 1
	exp = LRbin(2, expexp, calculate_e(expexp))
	t = PowerModm(y, exp, p) # t = y**(2**(r - m - 1))
	y = (t**2) % p
	r = m % p
	x = (x * t) % p
	b = (b * y) % p
	return findexp(b, p, x, r, a, y)