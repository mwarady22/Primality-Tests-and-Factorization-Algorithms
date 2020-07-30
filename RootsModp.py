from sage.all import *
import sys
from SqrtModp import *

def RootsModp(p, poly, field):
	retlist = []
	return #

def isolateroots(p, poly, field=0, retlist=[]):
	R = ZZ['x']
	S = R.quotient(poly, 'x')
	x = S.gen()
	xp = x**p
	xpx = xp - x
	print(xpx)
	print('R')
	print(R)
	print('S')
	print(S)
	l = xpx.number_of_terms()
	xpx0 = R.gen()
	print(xpx0)
	# for term in range(0, l):



	return xpx0





# 	S = ZZ['a'].quotient(poly, 'x')
# 	print('S')
# 	print(S)
# 	x = S.gen()
# 	print('x')
# 	print(x)
# 	XpminusX = x**p - x
# 	print('XpminusX')
# 	print(XpminusX)
# 	print('poly')
# 	print(poly)
# 	# g = XpminusX.gcd(poly)
# 	R = parent(poly)
# 	print('R')
# 	print(R)
# 	print('XpminusX')
# 	print(XpminusX)

# 	print('parent(poly)')
# 	print(parent(poly))
# 	print('parent(XpminusX)')
# 	print(parent(XpminusX))

# 	R(XpminusX)
# 	g = poly.gcd(XpminusX)
# 	g0 = g(0)
# 	if g0 == 0:
# 		retlist.append(0)
# 		g = g / x
# 	return g
	# return smalldeg() #

def smalldeg(p, poly, field, retlist, g): #
	if degree(g) == 0:
		return retlist
	elif degree(g) == 1:
		a1, a0 = g.coefficient(1, 0)
		retlist.append(- a1 / a0)
		return retlist
	elif degree(g) == 2:
		a2, a1, a0 = g.coefficient(2, 1, 0)
		d = a1**2 - 4 * a0 * a2
		f = SqrtModp(d, p)
		retlist.append((- a1 + f) / (2 * a2))
		retlist.append((- a1 - f) / (2 * a2))
		return retlist
	return randomsplitting() #

def randomsplitting(p, poly, field, retlist, g): #
	# a = ####
	gcdelem1 = (x + a)**((p - 1) / 2) - 1
	B = gcd(gcdelem1, g)
	if (degree(B) == 0) or (degree(B) == degree(g)):
		return randomsplitting(p, poly, field, retlist, g)
	return recurse() #

def recurse(): #
	listB = RootsModp(p, B, field)
	listgbyB = RootsModp(p, g / B, field)
	for item in listB:
		retlist.append(item)
	for item in listgbyB:
		retlist.append(item)
	return retlist




R = ZZ['x']
x = R.gen()
poly = x**4 + 10 * x + 4

print(isolateroots(3, poly))