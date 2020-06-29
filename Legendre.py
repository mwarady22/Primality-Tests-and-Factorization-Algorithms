from sage.all import *
import sys

def Kronecker(a, b): # computes the Kronecker symbol for a, b
	if b == 0:
		if (abs(a) == 1):
			return 1
		else:
			return 0

	if (a % 2 == 0) and (b % 2 == 0):
		return 0

	v = 0
	while (b % 2 == 0):
		v += 1
		b = b / 2
	if (v % 2 == 0):
		k = 1
	else:
		exp = (a**2 - 1) / 8
		k = (- 1)**exp
	if b < 0:
		b = - b
		if a < 0:
			k = - k
	return finishedcheck(a, b, k)

def finishedcheck(a, b, k):
	if (a == 0):
		if b > 1:
			return 0
		if b == 1:
			return int(k)
	v = 0
	while (a % 2 == 0):
		v += 1
		a = int(a / 2)
	if (v % 2 == 1):
		exp = (b**2 - 1) / 8
		k = (- 1)**exp * k
	return reciprocity(a, b, k)

def reciprocity(a, b, k):
	exp = (a - 1) * (b - 1) / 4
	k = (- 1)**exp * k
	r = abs(a)
	a = b % r
	b = r
	return finishedcheck(a, b, k)

# print(Kronecker(-1, 5))