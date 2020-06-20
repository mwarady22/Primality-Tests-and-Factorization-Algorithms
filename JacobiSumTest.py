from sage.all import *
import sys
from PrimesList import *


def ecalc(t):
	# v_q(t) denotes the highest exponent k such that q**k | t
	if (t % 2 == 1):
		return 2
	ret = 2
	q = primelist[0]
	while q - 1 < t:
		if t % (q - 1) == 0:
			r = 0
			while (t % (q**r)) == 0:
				r *= q
			ret *= r
	return ret


def etable(tablebound): ######### maybe just make list of like 3 see page 464 ########
	table = []
	index = 0
	while index < tablebound:
		table.append(ecalc(index))
	return table


def Precomps(B):


def findt(B):
	table = etable(sqrt(B))
	index = 0
	while table[index]**2 <= B:
		index += 1
	return index