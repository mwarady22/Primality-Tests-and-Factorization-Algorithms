from sage.all import *
import sys

def Wilson(p):
	ans = 1
	for j in range(1, p): # compute (p - 1)! mod p
		ans *= j
		ans = ans % p
	if ans == (p - 1): # if (p - 1)! == - 1 mod p, p is prime
		return 'prime'
	else: # otherwise p is composite
		return 'composite'