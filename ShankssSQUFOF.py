from sage.all import *
import sys
from Rabin_Miller import *
from PowerTests import *

def SQUFOF(N):
	print('SQUFOF')
	if N % 2 == 0:
		return 2
	rm = RM(N)
	if rm == 'proably prime':
		return rm
	else:
		a, b, c, d = precomp()
		square = SquareTest(N, a, b, c, d)
		if square == 'not square':
			return initialization(N)
		else:
			return square

def initialization(N):
	print('init')
	if N % 4 == 1:
		print('case 1')
		D = N
		d = IntSqrt(D) # the floor of the square root of D
		b = 2 * floor((d - 1) / 2) + 1
	else:
		print('case 2')
		D = 4 * N
		d = IntSqrt(D) # the floor of the square root of D
		b = 2 * floor(d / 2)
	f = (1, b, (b**2 - D) / 4)
	print(b)
	print(D)
	Q = [] # the queue
	h = 0
	L = IntSqrt(d)
	return applyrho(N, D, d, b, f, Q, h, L)

def applyrho(N, D, d, b, f, Q, h, L):
	print('applyrho')
	a, b, c = f
	f = rho(a, b, c, D)
	h += 1
	if h % 2 == 1:
		return fillq(N, D, d, b, f, Q, h, L)
	else:
		return squareform(N, D, d, b, f, Q, h, L)

def rho(a, b, c, D):
	print('rho')
	c0 = c
	if c < 0:
		c = - c
	r = - b % (2 * c)
	if (abs(c) > sqrt(D)):
		while (r <= - abs(c)):
			r += 2 * c
		while (r > abs(c)):
			r -= 2 * c
	elif (abs(c) < sqrt(D)):
		while (r <= sqrt(D) - 2 * abs(c)):
			r += 2 * c
		while (r >= sqrt(D)):
			r -= 2 * c
	rh = (c0, r, (r**2 - D) / (4 * c0))
	return rh


def squareform(N, D, d, b, f, Q, h, L):
	print('squareform')
	A, B, C = f
	q11, q63, q64, q65 = precomp()
	square = SquareTest(N, q11, q63, q64, q65)
	if isinstance(square, int):
		a = square
		if a in Q:
			return shortperiod(N, D, d, b, f, Q, h, L)
		else:
			return initbackcycle(N, D, d, b, f, Q, h, L, a)
	return shortperiod(N, D, d, b, f, Q, h, L)

def shortperiod(N, D, d, b, f, Q, h, L):
	print('shortperiod')
	A, B, C = f
	if A == 1:
		return 'ran through principle cycle, did not find non-trivial squareform'
	return fillq(N, D, d, b, f, Q, h, L)

def fillq(N, D, d, b, f, Q, h, L):
	print('fillq')
	A, B, C = f
	if abs(A) <= L:
		Q.append(abs(A))
	return applyrho(N, D, d, b, f, Q, h, L)

def initbackcycle(N, D, d, b, f, Q, h, L, a):
	print('initbackcycle')
	A, B, C = f
	s = gcd(a, B, D)
	if s > 1:
		return s**2
	g = (a, - B, a * C)
	a, b, c = g
	while not (abs(sqrt(D) - 2 * abs(c)) < b and b < sqrt(D)):
		u, v, w = g
		g = rho(u, v, w, D)
	return backcycle(N, D, d, b, f, Q, h, L, a, g)

def backcycle(N, D, d, b, f, Q, h, L, a, g):
	print('backcycle')
	b1 = b
	a, b, c = g
	g = rho(a, b, c, D)
	a, b, c = g
	if b1 != b:
		return backcycle(N, D, d, b, f, Q, h, L, a, g)
	else:
		if a % 2 == 1:
			return abs(a)
		else:
			return abs(a / 2)



print(SQUFOF(105))