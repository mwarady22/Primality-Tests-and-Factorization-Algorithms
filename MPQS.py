from sage.all import *
import sys
from Primeslist import *
from Legendre import *
from SqrtModp import *


def MPSQ(N): #
	return #






def chooseM(): # choose seiving length of 2 * M #
	return #

def chooseA(): # choose A for polynomial Q(x) = A * x**2 + 2 * B * x + C #
	approx = sqrt(2 * N) / M
	index = 0
	while primelist[index] < approx:
		index += 1
	A = primelist[index]
	if Kronecker(N, A) == 1:
		return A
	else:
		indexminus = index
		indexplus = index
		while True:
			indexminus -= 1
			indexplus += 1
			if indexminus >= 0:
				Aminus = primelist[indexminus]
			if indexplus < len(primelist):
				Aplus = primelist[indexplus]
			if Kronecker(N, Aminus) == 1:
				A = Aminus
				break
			elif Kronecker(N, Aplus) == 1:
				A = Aplus
				break
			if (indexminus >= 0) and (indexplus < len(primelist)):
				# choose new M
	return A

def chooseB(): # choose B for polynomial Q(x) = A * x**2 + 2 * B * x + C #
	B = SqrtModp(N, A)
	return B


def chooseC(): # choose C for polynomial Q(x) = A * x**2 + 2 * B * x + C #
	C = (B**2 - N) / A
	return C