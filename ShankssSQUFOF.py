from sage.all import *
import sys
from Rabin_Miller import *
from PowerTests import *
from EuclideanAlgorithm import *

def SQUFOF(N, k):
	return init0(N, k) #

def init0(N, k):
	q11, q63, q64, q65 = precomp()
	P0 = floor(sqrt(k * N))
	Pj = P0
	Q0 = 1
	Qj = Q0
	Q1 = k * N - P0**2
	Qh = Q1
	square = 1 # not square yet
	Pj, Qj, Qh = loop0(P0, Pj, Qj, Qh)
	while square == 1:
		Pj, Qj, Qh = loop0(P0, Pj, Qj, Qh)
		Pj, Qj, Qh = loop0(P0, Pj, Qj, Qh)
		square = squaretest(Qh, q11, q63, q64, q65)
	return init1(N, k, P0, Pj, Q0, Qh) #

def loop0(P0, Pj, Qj, Qh): # h = j + 1, k = h + 1
	bh = floor((P0 + Pj) / Qh)
	Ph = bh * Qh - Pj
	Qk = Qj + bh * (Pj - Ph)
	return Ph, Qh, Qk

def squaretest(Qk, q11, q63, q64, q65): # return 0 if square, 1 if not square
	sq = SquareTest(Qk, q11, q63, q64, q65)
	if isinstance(sq, str):
		return 1
	else:
		return 0

#############################################################################

def init1(N, k, P0, Pj, Q0, Qh):
	b0 = floor((P0 - Pj) / sqrt(Qh))
	P0 = b0 * sqrt(Qh) + Pj
	Pj = P0
	oldPj = None
	Q0 = sqrt(Qh)
	Qj = Q0
	Q1 = (k * N - P0**2) / Q0
	Qh = Q1
	while oldPj != Pj:
		oldPj = Pj
		Pj, Qj, Qh = loop1(P0, Pj, Qj, Qh)
	return outputcheck(N, k, Pj)

def loop1(P0, Pj, Qj, Qh):
	bh = floor((P0 + Pj) / Qh)
	Ph = bh * Qh - Pj
	Qk = Qj + bh * (Pj - Ph)
	return Ph, Qh, Qk

def outputcheck(N, k, Ph):
	f = gcd(N, Ph)
	if f != 1:
		return f
	# else:
		return SQUFOF(N, k + 1)


# print(SQUFOF(88526117, 1))