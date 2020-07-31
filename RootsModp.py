from sage.all import *
import sys
from SqrtModp import *

def RootsModp(p, poly):
	R = PolynomialRing(GF(p), 'x')
	x = R.gen()
	poly = x**4 + 10 * x + 4

print(poly.roots())