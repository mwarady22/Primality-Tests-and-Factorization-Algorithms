from sage.all import *
import sys
from EuclideanAlgorithm import *

# Takes in simultaneous congruences in coprime moduli and outputs their unique solution.
# If the moduli are not coprime, it will output 'None'.
def CRT_basic(xlist, mlist): # takes in 2 lists, the first containing the integers and the second containing the moduli

	if len(xlist) != len(mlist): # list lengths must match
		return None

	M = 1 # the value m_1 * m_2 * .... * m_n, starting as an empty product
	for m in mlist:
		M *= m

	alist = [] # this list will hold the inverses of each M_i = M / m_i mod m_i
	for m_i in mlist:
		M_i = M / m_i
		B = Bezout(M_i, m_i) # returns a list [u, x, v, y, d] such that ux + vy = d
		if B[4] != 1: # in this case the m_i are not all coprime
			# print(None)
			return None
		if B[1] == M_i:
			alist.append(B[0])
		else:
			alist.append(B[2])

	x = 0 # the return value, starting as an empty sum
	for j in range(len(mlist)):
		x += alist[j] * (M / mlist[j]) * xlist[j]
	x = x % M

	return int(x)

# Computes the Clist needed to compute multiple congruences mod the same mlist.
def CRT_precomp(mlist):

	Clist = [1]
	j = 2

	while j <= len(mlist):
		p = 1 # an empty product that will contain m_1 * m_2 * .... * m_j-1 mod m_j
		m_j = mlist[j - 1]
		for l in range(0, j - 1):
			p *= mlist[l]
		p = p % m_j
		B = Bezout(p, m_j) # returns a list [u, x, v, y, d] such that ux + vy = d
		if B[4] != 1: # in this case the m_i are not all coprime
			# print(None)
			return None
		if B[1] == p:
			Clist.append(B[0])
		else:
			Clist.append(B[2])
		j += 1

	return Clist

# Single congruence computation using the advanced method.
def CRT_comp(xlist, mlist, Clist):

	if len(xlist) != len(mlist): # list lengths must match
		return None

	if Clist == None:
		return None

	ylist = [xlist[0] % mlist[0]]
	for j in range(1, len(mlist)):
		y_j = ylist[len(ylist) - 1]
		for k in range(j - 2, - 1, - 1):
			y_j *= mlist[k]
			y_j += ylist[k]
		y_j = (xlist[j] - y_j) * Clist[j]
		ylist.append(y_j % mlist[j])

	x = ylist[len(ylist) - 1]
	for l in range(len(ylist) - 2, - 1, - 1):
		x *= mlist[l]
		x += ylist[l]

	return x

# Runs multiple congruence computatins using the advanced method.
def CRT_multiple(xlistlist, mlist):

	Clist = CRT_precomp(mlist)

	if Clist == None:
		return None

	outlist = []

	for xlist in xlistlist:
		outlist.append(CRT_comp(xlist, mlist, Clist))

	return outlist

def CRT_ind(xlist, mlist):

	if len(xlist) != len(mlist): # list lengths must match
		return None

	x = xlist[0]
	m = mlist[0]
	j = 1

	while j < len(xlist):
		B = Bezout(m, mlist[j])
		if B[4] != 1:
			return None
		if B[1] == m:
			u = B[0]
			v = B[2]
		else:
			u = B[2]
			v = B[0]
		x = u * m * xlist[j] + v * mlist[j] * x
		m = m * mlist[j]
		x = x % m
		j += 1

	return x