from sage.all import *
import sys

# Computes the greatest common divisor of two integers using the Euclidean Algorithm.
# The gcd does not exist if both integers are 0.
def gcd(a, b):

	a = abs(int(a)) # gcd of negatives is the same as that of their absolute values
	b = abs(int(b))

	if (a == 0) and (b == 0): # gcd of 0 and 0 does not exist
		return None

	if b > a: # order does not matter for gcd
		a, b = b, a

	div = a # the gcd that we will eventually return
	while b != 0:
		q = 0
		div = b
		r = a % b
		a, b = b, r

	return div


def Bezout(a, b):
	a = int(a)
	b = int(b)

	if (a == 0) and (b == 0): # gcd of 0 and 0 does not exist
		return None

	aneg = 1
	bneg = 1

	if a < 0: # we only work with a, b positive and reintroduce negatives to their coefficients later if necessary
		aneg = - 1
		a = abs(a)
	if b < 0:
		bneg = - 1
		b = abs(b)

	if b > a: # order of a, b does not matter so we put the larger one first
		a, b = b, a
		aneg, bneg = bneg, aneg

	# a1 = a
	# b1 = b

	lst = []

	while b != 0:
		q = 0
		while (a - q * b >= b):
			q += 1
		lst.append([a, q, b, a - q * b])
		a, b = b, a - q * b
	lst.pop() # get rid of last equation ending in 0

	# p = m * x + n * y
	m, a, n, b, p = 0, 0, 0, 0, 0

	loop = 0 # even or odd iteration of the while loop, 0 indexed

	while lst != []:
		eqn = lst.pop()
		if loop == 0:
			m = 1
			a = eqn[0]
			n = - eqn[1]
			b = eqn[2]
			p = eqn[3]
		elif loop % 2 == 0:
			# m = m
			a = eqn[0]
			n = n + (- eqn[1] * m)
			# b = b
		else:
			m = m + (- eqn[1] * n)
			# a = a
			# n = n
			b = eqn[0]

		loop += 1

	if loop % 2 == 0:
		return [bneg * m, bneg * a, aneg * n, aneg * b, p]
	else:
		return [aneg * m, aneg * a, bneg * n, bneg * b, p]


# Control

# print(gcd(*sys.argv[1 : ]))

# print(Bezout(*sys.argv[1 : ]))