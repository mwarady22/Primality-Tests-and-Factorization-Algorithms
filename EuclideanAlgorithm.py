from sage.all import *
import sys

def gcd(a, b):
	a = int(a)
	b = int(b)

	a = abs(a)
	b = abs(b)

	if b > a:
		a, b = b, a

	div = a
	while b != 0:
		q = 0
		div = b
		while (a - q * b >= b):
			q += 1

		a, b = b, a - q * b

	print(div)
	return div


def Bezout(a, b):
	a = int(a)
	b = int(b)

	aneg = 1
	bneg = 1

	if a < 0:
		aneg = - 1
		a = abs(a)
	if b < 0:
		bneg = - 1
		b = abs(b)

	if b > a:
		a, b = b, a
		aneg, bneg = bneg, aneg

	a1 = a
	b1 = b

	lst = []

	while b1 != 0:
		q = 0
		while (a1 - q * b1 >= b1):
			q += 1
		lst.append([a1, q, b1, a1 - q * b1])
		a1, b1 = b1, a1 - q * b1
	lst.pop() # get rid of last equation ending in 0

	# p = m * x + n * y
	m, a1, n, b1, p = 0, 0, 0, 0, 0

	loop = 0 # is this an even or odd iteration of the while loop, 0 indexed

	while lst != []:
		eqn = lst.pop()
		if loop == 0:
			m = 1
			a1 = eqn[0]
			n = - eqn[1]
			b1 = eqn[2]
			p = eqn[3]
		elif loop % 2 == 0:
			# m = m
			a1 = eqn[0]
			n = n + (- eqn[1] * m)
			# b1 = b1
		else:
			m = m + (- eqn[1] * n)
			# a1 = a1
			# n = n
			b1 = eqn[0]

		loop += 1

	if loop % 2 == 0:
		print([bneg * m, bneg * a1, aneg * n, aneg * b1, p])
		return [bneg * m, bneg * a1, aneg * n, aneg * b1, p]
	else:
		print([aneg * m, aneg * a1, bneg * n, bneg * b1, p])
		return [aneg * m, aneg * a1, bneg * n, bneg * b1, p]


# Control

# gcd(*sys.argv[1 : ])

Bezout(*sys.argv[1 : ])