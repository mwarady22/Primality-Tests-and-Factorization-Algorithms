from sage.all import *
import sys
from PoweringAlgorithms import *
from EuclideanAlgorithm import *

def AKSTest(N): # deterministically outputs if N is prime or not

	for a in range(2, floor(sqrt(N))): # check if N is a perfect power
		d = 1
		while d < N:
			d *= a
			if d == N:
				return 'composite'

	logNlogN = (log(N)**2).n() # log**2(N)

	r = 2

	for R in range(ceil(log(N).n()), N): # finding smallest r such that ord(N) (mod r) > log**2(N)
		ZmodRZ = Integers(R)
		if gcd(R, N) == 1:
			if ZmodRZ(N).multiplicative_order() > logNlogN:
				r = R
				break

	for b in range(2, r): # check if gcd(b, N) is greater than 1
		if gcd(N, b) > 1:
			return 'composite'

	if N >= r:
		return 'prime'

	phir = euler_phi(r) # phi(r)
	floorphirlogN = floor(phir * (log(N).n())) # floor(phi(r) * log(N))

	x = var('x')
	R = Integers(N)[x]
	f = x**N - 1
	g = N
	J = R.ideal(f)
	S = R.quotient(J, 'x') # S is ZZ / (N, x**r - 1)

	for c in range(1, floorphirlogN):
		print(c)
		LHS = ((x + c)**N).expand() # calculate left hand side polynomial
		RHS = x**N + c # calculate right hand side polynomial
		LHS1 = S(LHS) # LHS % (N, x**r - 1)
		RHS1 = S(RHS) # RHS % (N, x**r - 1)
		if LHS1 != RHS1: # check if they are the same
			return 'composite'
	return 'prime'


# print(AKSTest(2099))