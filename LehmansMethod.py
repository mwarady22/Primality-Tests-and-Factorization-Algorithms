from sage.all import *
import sys
from TrialDivisionFactoring import *
from Squares import *
from EuclideanAlgorithm import *

def Lehman(N): # finds a factor of N

	B = floor(N**(1 / 3)) # first checks by trial division that there is not a small factor ....
	trial = TrialFactor(N, B) # .... this ensures that there are at most two factors of N
	if isinstance(trial, list):
		return trial[0]
	else:
		return loopk(N, B)

def loopk(N, B, k=0):
	k += 1
	if k > B:
		return 'prime'
	else:
		if k % 2 == 0:
			r = 1
			m = 2
		else:
			r = k + N
			m = 4
		return loopa(N, B, k, r, m)

def loopa(N, B, k, r, m):
	low = 4 * k * N
	a = IntSqrt(low)
	if a**2 < low: # find first integer a such that 4 * k * N <= a**2 < 4 * k * N + B**2
		a += 1
	while (a**2 < (4 * k * N + B**2)) and ((a % m) == (r % m)):
		c = a**2 - 4 * k * N
		b = IntSqrt(c)
		if isinstance(b, int):
			return gcd(a + b, N)
	return loopk(N, B, k)


# print(Lehman(25))