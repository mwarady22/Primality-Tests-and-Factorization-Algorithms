from sage.all import *
import sys

# Right-Left Binary powering algorithm.  Computes g**n efficiently. 
def RLbin(g, n):

	if n == 0:
		# print(1)
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
			# print(y)
			return y
		z *= z

	

# Left-Right Binary powering algorithm.  Computes g**n efficiently.  Utilizes the unique integer
# e such that 2**e <= |n| < 2**e+1.
def LRbin(g, n, e):

	if n == 0:
		# print(1)
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
			# print(y)
			return y
		else:
			E = E / 2
		y *= y
		if N >= E:
			N = N - E
			y *= z



# Control

# RLbin(2, -2)
# LRbin(2, -3, 1)