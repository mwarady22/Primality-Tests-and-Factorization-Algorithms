from sage.all import *
import sys
from Legendre import *
from SqrtModp import *
from PowerTests import *

def Cornacchia(p, d): # p prime, 0 < d < p, outputs solution to x**2 + d * y**2 = p or says there is none
	k = Kronecker(- d, p) # determine if there can by a solution
	if k == - 1:
		return 'no solution'
	x0 = SqrtModp(- d, p)
	while x0 < p / 2: # find number congruent to +-x0 between p / 2, p
		x0 += p
	while x0 > p:
		x0 -= p
	a = p
	b = x0
	l = floor(sqrt(p))
	b, l, a = euclidalg(b, l, a)
	c = int((p - b**2) / d)
	q11, q63, q64, q65 = precomp()
	if (not ((p - b**2)  % d == 0)) or isinstance(SquareTest(c, q11, q63, q64, q65), str):
		return 'no solution'
	else:
		return b, sqrt(c)

def euclidalg(b, l, a):
	if b > l:
		r = a % b
		a = b
		b = r
		return euclidalg(b, l, a)
	return b, l, a

# print(Cornacchia(5, 2))