from sage.all import *
import sys
from TrialDivisionFactoring import *
from Legendre import *
from Cornacchia import *
from TrialDivisionFactoring import *
from PollardspMinus1 import *



# Dlist



def AtkinTest(N): # N != 1 is coprime to 6 and passes the Rabin-Miller test
	h = 0
	nn = 0
	Ni = N
	#







def Nismall(Ni): # 2
	if Ni < 2**20: # if Ni is small enough just trial divide
		t = TrialFactor(Ni, 2**10)
		if not isinstance(t, str): # then it is not the case that Ni is prime or has only factors larger than 2**15
			return # 14
		return 'prime' #
	return # 3

def nextD(): # 3 #
	nn += 1
	D = Dlist[nn]
	if Kronecker(D, N) != 1:
		return # 3
	else:
		sol = Cornacchia(4 * N, abs(D)) # the solution (x, y) to x**2 + abs(D) * y**2 = 4 * N or 'no  solution'
		if isinstance(sol, str):
			return # 3
		else:
			x = sol[0]
			y = sol[1]
			return # 4

def factorm(): # 4 #
	m0 = N + 1 + x
	m1 = N + 1 - x
	####
	m2 = N + 1 + 2 * y
	m3 = N + 1 - 2 * y
	####
	m2 = N + 1 + (x + 3 * y) / 2
	m3 = N + 1 - (x + 3 * y) / 2
	m4 = N + 1 + (x - 3 * y) / 2
	m5 = N + 1 - (x - 3 * y) / 2
	####
	