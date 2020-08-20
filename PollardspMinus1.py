from sage.all import *
import sys
from PrimesList import *
from EuclideanAlgorithm import *
from PoweringAlgorithms import *

def FirstStage(N, B): # runs independently, finds a factor of N only if there is one that is B-powersmooth
	x = 2
	g = gcd(x, N)
	if g == N:
		return 'must use x < N'
	elif g > 1:
		return g

	for h in range(2, B + 1):
		x = PowerModm(x, h, N)
		g = gcd(x - 1, N)
		if g == N:
			return 'the algorithm has failed' # x does not work, can try another value
		elif g > 1:
			return g
		h += 1
	return 'the algorithm did not succeed in spliting N' # B not high enough


def FirstStageforTwoStage(N, B, x=2): # same as FirstStage, but also outputs x for use in SecondStage

	if x >= N:
		x = 2
	g = gcd(x, N)
	if g == N:
		return 'must use x < N'
	elif g > 1:
		return g

	for h in range(2, B + 1):
		x = PowerModm(x, h, N)
		g = gcd(x - 1, N)
		if g == N:
			return 'the algorithm has failed : x = ' + str(x) # x does not work
		elif g > 1:
			return g
		h += 1
	return 'the algorithm did not succeed in spliting N : x = ' + str(x) # B not high enough

def secondstageprecomps(B1, B2): # precomputes 
	z = 0
	while primelist[z] < B1: # find index of first prime larger than B1
		z += 1
	k1 = z
	while primelist[z] < B2: # find index of last prime smaller than B2
		z += 1
	k2 = z - 1

	dlist = [] # list that will hold the differences of consecutive primes between B1, B2
	y = k1
	while y + 1 < k2: # fill dlist
		dlist.append(primelist[y + 1] - primelist[y])
		y += 1
	return k1, k2, dlist

def SecondStage(N, B1, B2):
	if B1 > B2: # make sure that B1 is smaller or equal to B2
		B1, B2 = B2, B1

	first = FirstStageforTwoStage(N, B1) # run FirstStage
	if isinstance(first, (int, Integer)): # if it returns an answer, return that
		return first
	elif first[14] == 'h': # get x from the result
		x = int(first[30 :])
	else: # get x from the result
		x = int(first[49 :])
	# initialize variables
	b = x
	c = 0
	P = 1
	h = 0
	j = h
	y = x

	k1, k2, prec = secondstageprecomps(B1, B2) # unpack precomputed variables
	bdlist = [] # list to hold b**d for every difference d in dlist
	for elem in prec:
		bdlist.append(b**elem)
	x = x**primelist[k1]
	return advance(N, k1, k2, bdlist, b, c, P, h, j, x, y)

def advance(N, k1, k2, bdlist, b, c, P, h, j, x, y): # go to next prime
	h += 1
	if (h >= k2 - k1): # if h exceeds possible indices of bdlist, return fail
		return fail(N, k1, k2, bdlist, b, c, P, h, j, x, y)
	x = x * bdlist[h - 1]
	P = P * (x - 1)
	c += 1
	if h >= k2:
		return fail(N, k1, k2, bdlist, b, c, P, h, j, x, y)
	elif c < 20:
		return advance(N, k1, k2, bdlist, b, c, P, h, j, x, y)
	return compgcd(N, k1, k2, bdlist, b, c, P, h, j, x, y)

def compgcd(N, k1, k2, bdlist, b, c, P, h, j, x, y): # compute gcd
	g = gcd(P, N)
	if g == 1:
		c = 0
		j = h
		y = x
		return advance(N, k1, k2, bdlist, b, c, P, h, j, x, y)
	return backtrack(N, k1, k2, bdlist, b, c, P, h, j, x, y, g)

def backtrack(N, k1, k2, bdlist, b, c, P, h, j, x, y, g): # backtrack and see if we have found a factor
	h = j
	x = y
	while not g > 1:
		x = x * bdlist[h - 1]
		h += 1
		g = gcd(x - 1, N)
	if g < N:
		return g
	else:
		return 'the algorithm has failed' # x does not work starting in the first stage

def fail(N, k1, k2, bdlist, b, c, P, h, j, x, y): # either fail or go to backtrack
	g = gcd(P, N)
	if g == 1:
		return 'the algorithm has failed'
	else:
		return backtrack(N, k1, k2, bdlist, b, c, P, h, j, x, y, g)