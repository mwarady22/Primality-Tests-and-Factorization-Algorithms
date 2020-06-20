from sage.all import *
import sys
from PoweringAlgorithms import *

def PrimRoot(p):
	print('primroot')
	a = 1
	q = p - 1
	plist = []
	for j in range(2, floor(sqrt(q)) + 1):
		if (q % j == 0):
			while (q % j == 0):
				q = int(q / j)
			plist.append(j)
	return initcheck(p, a, plist) #

def initcheck(p, a, plist): #
	print('initcheck')
	a += 1
	h = 0
	return checkph(p, a, h, plist) #

def checkph(p, a, h, plist): #
	print('checkph')
	exp = int((p - 1) / plist[h])
	aexp = PowerModm(a, exp, p)
	if aexp == 1:
		return initcheck(p, a, plist) # 2
	h += 1
	return finishedcheck(p, a, h, plist)

def finishedcheck(p, a, h, plist): #
	print('finishedcheck')
	if h > len(plist) - 1:
		return a
	else:
		return checkph(p, a, h, plist) # 3

print(PrimRoot(7))