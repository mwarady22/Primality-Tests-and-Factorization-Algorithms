from sage.all import *
import sys
from PoweringAlgorithms import *

def PrimRoot(p): # finds number g that, taken to various powers, gives each residue mod p
	a = 1
	q = p - 1
	plist = []
	for j in range(2, floor(sqrt(q)) + 1): # make list of prime factors
		if (q % j == 0):
			while (q % j == 0):
				q = int(q / j)
			plist.append(j)
	return initcheck(p, a, plist)

def initcheck(p, a, plist): # check new possible primitive root
	a += 1
	h = 0
	return checkph(p, a, h, plist)

def checkph(p, a, h, plist): # check new exponent
	exp = int((p - 1) / plist[h])
	aexp = PowerModm(a, exp, p)
	if aexp == 1:
		return initcheck(p, a, plist) # not a primitive root, check new base
	h += 1
	return finishedcheck(p, a, h, plist)

def finishedcheck(p, a, h, plist): # check if we have found a primitive root
	if h > len(plist) - 1:
		return a # done, return
	else:
		return checkph(p, a, h, plist) # keep checking exponents

# print(PrimRoot(7))