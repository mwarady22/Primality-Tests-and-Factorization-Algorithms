from sage.all import *
import sys
from PrimesList import *

def SL(N): #

	B = floor((L(N))**.5)
	plist = getplist(B)
	return SLinit(N) #
	# D = Kinit(K, N)
	# x, c, h = chooseform() #
	# q, q1, l = nextprime() #



def getplist(B):
	l = []
	c = 0
	while primelist[c] <= B:
		l.append(primelist[c])
	return l

def L(x):
	return e**(sqrt(ln(x) * ln(ln(x))))

def SLinit(N): # 1
	B = floor((L(N))**.5)
	K = 1
	d = floor(ln(B))
	return # 2

def Kinit(K, N): # 2
	if ((K * N % 4) == 3):
		D = - K * N
	else:
		D = - 4 * K * N
	return # 3

def chooseform(): # 3 #
	####
	return # 4

def nextprime(): # 4 #
	h += 1
	if h > k:
		K += 1
		return Kinit() #
	else:
		q = plist[h]
		q1 = q
		l = floor(B / q)
		return # 5

def comppower(): # 5 #
	while q1 <= l:
		q1 *= 1
	x = x**q1 #### power in class group ####
	c += 1
	if c < 20:
		return # 4
	else:
		return # 6

def checksuccess(): # 6 #
	d1 = 0
	while x # not ambiguous form and d1 < d: #
		x *= x
		d1 += 1
	if x # not ambiguous form: #
		c = 0
		return # 4
	else:
		return # 7

def checkfinished(): # 7 #
	####