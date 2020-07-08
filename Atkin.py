from sage.all import *
import sys
from TrialDivisionFactoring import *
from Legendre import *



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
		return 'prime'

def nextD(): # 3 #
	nn += 1
	D = Dlist[nn]
	if Kronecker(D, N) != 1:
		return # 3
	else:
		# Cornacchia