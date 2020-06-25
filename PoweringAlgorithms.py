from sage.all import *
import sys

# Right-Left Binary powering algorithm.  Computes g**n efficiently. 
def RLbin(g, n):

	if n == 0:
		return 1
	elif n < 0:
		N = - n
		z = g**-1
	else:
		N = n
		z = g

	y = 1

	while True:
		if N % 2 == 1:
			y = z * y

		N = floor(N / 2)
		if N == 0:
			return y
		z *= z

	

# Left-Right Binary powering algorithm.  Computes g**n efficiently.  Utilizes the unique integer
# e such that 2**e <= |n| < 2**e+1.
def LRbin(g, n, e):

	if n == 0:
		return 1
	elif n < 0:
		N = - n
		z = g**-1
	else:
		N = n
		z = g

	y = z
	E = 2**e
	N = N - E

	while True:
		if E == 1:
			return y
		else:
			E = E / 2
		y *= y
		if N >= E:
			N = N - E
			y *= z


def calculate_e(n):
	if n == 0:
		return 0
	return floor(log(abs(n), 2))


def PowerModm(g, n, m):
	if n == 0:
		return 1 % m
	if n == 1:
		return g % m

	nbin = bin(n)[2:]
	nbinlen = len(nbin)

	g1 = g % m
	l = [g1]
	for j in range(1, nbinlen):
		l.insert(0, (l[0]**2 % m))

	product = 1
	for k in range(0, nbinlen):
		product = product * (l[k]**int(nbin[k])) % m

	return product



# Control

# print(RLbin(2, -2))
# print(LRbin(2, -3, 1))
# print(calculate_e(65))
# print(PowerModm(32254, 548522, 172522))