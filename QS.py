from sage.all import *
import sys
from PrimesList import *
from Legendre import *


def QS(N, mult=1):
	baselist, q = collectbase(N) # list of base elements, product of all base elements
	print('baselist')
	print(baselist)
	x = ZZ['x'].gen()
	Qa = (x + ceil(sqrt(N)))**2 - N
	l = []
	q = []
	for z in range(0, ceil(sqrt(N) * mult)):
		Qaz = Qa(z)
		l.append(Qaz)
		q.append([z + 124, Qaz])
	m = l.copy()
	for p in baselist:
		print('for ' + str(p))
		w = Integers(p)['w'].gen()
		Ya = (w + ceil(sqrt(N)))**2 - N
		print('Ya')
		print(Ya)
		rootlist = Ya.roots()
		print('rootlist')
		print(rootlist)
		
		for root, mult in rootlist:
			print('len(l)')
			print(len(l))
			index = int(root)
			print('index')
			print(index)
			counter = 0
			print('index < len(l)')
			print(index < len(l))
			while index < len(l):
				print('counter : ' + str(counter))
				print('l')
				print(l)
				print('l[index]')
				print(l[index])
				l[index] = l[index] / p
				print('l[index]')
				print(l[index])
				index += p
				print('index')
				print(index)
				counter += 1
				print('index < len(l)')
				print(index < len(l))
	R = []
	for ind in range(0, len(l)):
		if l[ind] == 1:
			R.append([[len(R)], basefactorize(m[ind], baselist)])
	print('R')
	print(R)
	dep = lindep(N, q, R) # factorization or None if we cannot yet find a factorization
	if dep != None:
		return dep # the factorization
	else:
		# return QS(N, mult * 2)
		None


def basefactorize(t, baselist): # factorize over baselist and return exponents
	l = len(baselist)
	explist = [0] * l
	for x in range(0, l):
		if t == 1:
			break
		while (t % baselist[x] == 0):
			t = t / baselist[x]
			explist[x] += 1
	return explist


def collectbase(D): # find base for given value of D
	b = floor(e**(sqrt(ln(D) * ln(ln(D))) / 2)) # find bound for base primes
	if b < 31: # never have too small a base
		b = 31
	plist = [2] # collect base primes N <= b such that Kroencker(N, p) = 1
	pindex = 1 # index 1 will start the loop at the second prime, 3
	p = primelist[pindex]
	q = 2 # the product of the primes in the base
	while p <= b:
		if Kronecker(D, p) == 1: # only use primes for which D is a square mod p
			plist.append(p)
			q *= p
		pindex += 1
		p = primelist[pindex] # go to next prime
	return plist, q


def lindep(N, Q, R): # find dependency, return either factorization or None
	print('lindep')
	if R == []:
		return None
	S = R.copy()
	width = len(S[0][1])
	height = len(S)
	finishedheight = 0 # how many rows have been put at top of the matrix and are therefore done
	for item in S: # take all vectors mod 2
		for s in range(0, width):
			item[1][s] = item[1][s] % 2
	for w in range(0, width): # iterate through width
		print('S')
		print(S)
		for h in range(finishedheight, height): # iterate through height below finished rows
			if S[h][1][w] == 1: # find row that contains 1 in the current column
				row = S[h] # save row
				S.remove(S[h]) # remove row
				S.insert(finishedheight, row) # insert row at top of unfinished rows, bottom of finished rows
				for lower in range(finishedheight + 1, height): # add current row to each lower row
					if S[lower][1][w] == 1:
						S[lower] = (S[lower][0] + S[finishedheight][0], [(S[lower][1][j] + S[finishedheight][1][j]) % 2 for j in range(0, width)])
						if S[lower][1] == [0] * width: # if any row has all 0s (meaning all rows multiplied together give all even exponents), check if it factors
							check = checknumbers(N, Q, R, S[lower][0])
							if check != None:
								return check # return a factorization
				finishedheight += 1
				break
	return None

def checknumbers(N, Q, R, rows): # multiply the appropriate numbers together, see if they will give a nontrivial solution
	print('checknumbers')
	print('rows')
	print(rows)
	x2 = 1
	y2 = 1
	for row in rows:
		Ak, Ak2 = Q[row]
		print('Ak')
		print(Ak)
		print('Ak2')
		print(Ak2)
		x2 = (x2 * Ak) % N # multiply all the Aks together mod N
		y2 = (y2 * Ak2) # multiply all the Ak**2s together
	y = sqrt(y2)
	if ((x2 % N) != (y % N)) and ((x2 % N) != (- y % N)): # check that x2 != +- y mod N
		return factorize(x2, y)
	else: # if x2 != +- y mod N we would get a trivial solution, so try new set of numbers
		return None

def factorize(x2, y): # find factorization of N
	r0 = x2 + y
	f0 = gcd(r0, N)
	return f0, int(N / f0)



print(QS(15347))