from sage.all import *
import sys
from EuclideanAlgorithm import *

# Takes in simultaneous congruences in coprime moduli and outputs their unique solution.
# If the moduli are not coprime, it will output 'None'.
def CRT_basic(xlist, mlist): # takes in 2 lists, the first containing the integers and the second containing the moduli

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

	# print(int(x))
	return int(x)






# Control

# CRT_basic(*sys.argv[1 : ])

CRT_basic([6, 3, 2], [7, 6, 12])