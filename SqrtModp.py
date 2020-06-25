from sage.all import *
import sys
import random
from Legendre import *
from PoweringAlgorithms import *

def twokcalc(r): # 2kcalc(p) denotes the highest exponent k such that 2**k | t
	v = 0
	while (r % (2**v)) == 0:
		v += 1
	v -= 1
	return v

def SqrtModp(a, p): # finds square root mod p for odd prime p
	k = twokcalc(p - 1)
	q = int((p - 1) / (2**k))
	while True:
		print('rand')
		n = randint(0, p - 1)
		kr = Kronecker(n, p)
		if kr == - 1:
			break
	z = PowerModm(n, q, p)
	y = z
	r = k
	print('r : ' + str(r))
	x = PowerModm(a, (q - 1) / 2, p)
	b = a * PowerModm(x, 2, p)
	x = (a * x) % p
	return findexp(b, p, x, r, a, y) #

def findexp(b, p, x, r, a, y): #
	if (b % p) == 1:
		return x
	m = 1
	twom = 2
	x = 0
	while True:
		print('True : ' + str(x))
		b2m = PowerModm(b, twom, p)
		if b2m == 1:
			break
		m += 1
		twom *= 2
		x += 1
	if m == r:
		return str(a) + ' is not a quadratic residue mod ' + str(p)
	return redexp(m, y, b, p, x, r, a) #

def redexp(m, y, b, p, x, r, a): #
	exp = r - m - 1
	print('exp : ' + str(exp))
	t = LRbin(y, exp, calculate_e(exp))
	y = t**2
	r = m
	x = x * t
	b = b * y
	return findexp(b, p, x, r, a, y) 

print(SqrtModp(8, 17))