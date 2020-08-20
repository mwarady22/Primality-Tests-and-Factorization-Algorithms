from sage.all import *
import sys

# Right-Left Binary powering algorithm.  Computes g**n efficiently. 
def RLbin(g, n):

	if n == 0: # g**0 = 1
		return 1
	elif n < 0:
		N = - n # make N positive
		z = g**-1 # use the inverse of g as z
	else:
		N = n
		z = g

	y = 1

	while True:
		if N % 2 == 1: # if N is odd, multiply it by z
			y = z * y

		N = floor(N / 2) # subtracts one if N is odd and divides out a factor of two either way
		if N == 0:
			return y
		z *= z # square z

# Left-Right Binary powering algorithm.  Computes g**n efficiently.  Utilizes the unique integer e such that 2**e <= |n| < 2**e+1.
def LRbin(g, n, e):

	if n == 0: # g**0 = 1
		return 1
	elif n < 0:
		N = - n # make N positive
		z = g**-1 # use the inverse of g as z
	else:
		N = n
		z = g

	y = z
	E = 2**e # largest factor of 2 less than n
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

def calculate_e(n): # returns e such that 2**e <= |n| < 2**e+1
	if n == 0:
		return 0
	return floor(log(abs(n), 2))

def PowerModm(g, n, m):
	if n == 0: # g**0 = 1
		return 1 % m
	if n == 1: # g**1 = g
		return g % m

	nbin = bin(n)[2:] # gets integer part of binary n
	nbinlen = len(nbin) # length of nbin

	g1 = g % m
	l = [g1] # list holding powers of g1 from 1 to nbinlen - 1
	for j in range(1, nbinlen): # fill l
		l.insert(0, (l[0]**2 % m))

	product = 1
	for k in range(0, nbinlen): # multiply g to various powers chosen based on binary n using the list made in l
		product = product * (l[k]**int(nbin[k])) % m

	return product