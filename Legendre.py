from sage.all import *
import sys

def Kronecker(a, b): # computes the Kronecker symbol for a, b
	if b == 0: # two cases when b == 0
		if (abs(a) == 1):
			return 1
		else:
			return 0

	if (a % 2 == 0) and (b % 2 == 0): # if both a or b is even, return 0
		return 0

	v = 0
	while (b % 2 == 0): # if b is even but a is not, divide the highest power of 2 out of b
		v += 1
		b = b / 2
	if (v % 2 == 0): # if the highest power of 2 divided out of b is even, k = 1
		k = 1
	else: # otherwise k = 1 or k = - 1 depending on a
		exp = (a**2 - 1) / 8
		k = (- 1)**exp
	if b < 0: # make b positive
		b = - b
		if a < 0: # if both a and b were negative, negate k
			k = - k
	return finishedcheck(a, b, k)

def finishedcheck(a, b, k):
	if (a == 0):
		if b > 1:
			return 0
		if b == 1:
			return int(k)
	v = 0
	while (a % 2 == 0): # divide the highest power of 2 out of a
		v += 1
		a = int(a / 2)
	if (v % 2 == 1): # k = k or k = - k depending on b
		exp = (b**2 - 1) / 8
		k = (- 1)**exp * k
	return reciprocity(a, b, k)

def reciprocity(a, b, k):
	exp = (a - 1) * (b - 1) / 4
	flexp = floor(exp)
	dec = exp - flexp
	expmod2 = flexp % 2
	exp = expmod2 + dec
	k = (- 1)**exp * k
	r = abs(a)
	a = b % r
	b = r
	return finishedcheck(a, b, k)