from sage.all import *
import sys
from PoweringAlgorithms import *
from EuclideanAlgorithm import *

def AKSTest(N): # deterministically outputs if N is prime or not

	for a in range(2, floor(sqrt(N))):
		d = 1
		while d < N:
			d *= a
			if d == N:
				print('1')
				return 'composite'

	logNlogN = (log(N)**2).n()

	r = 2

	for R in range(ceil(log(N).n()), N):
		ZmodRZ = Integers(R)
		if gcd(R, N) == 1:
			if ZmodRZ(N).multiplicative_order() > logNlogN:
				r = R
				break

	for b in range(2, r):
		if gcd(N, b) > 1:
			print('2')
			return 'composite'

	# if N >= r:
	# 	print('3')
	# 	return 'prime'

	phir = euler_phi(r)
	floorphirlogN = floor(phir * (log(N).n()))

	# R = PolynomialRing(ZZ, 'X')
	# Y = R.gen()
	# S = R.quotient(Y**r - 1)
	# X = S.gen()

	R = PolynomialRing(ZZ, 'X')
	x = R.gen()
	f = x**N - 1
	g = N
	J = R.ideal(f, g)
	S = R.quotient(J)
	print('S')
	print(S)
	X = S.gen()

	for c in range(1, floorphirlogN):
		LHS = (X + c)**N
		RHS = X**N + c
		if LHS != RHS:
			print('LHS')
			print(LHS)
			print('RHS')
			print(RHS)
			print('4')
			return 'composite'

	return 'prime'


print(AKSTest(31))