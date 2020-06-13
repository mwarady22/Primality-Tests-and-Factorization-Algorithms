from sage.all import *
import sys

def Wilson(p):
	ans = 1
	for j in range(1, p): # compute (p - 1)! % p
		ans *= j
		ans = ans % p
	if ans == (p - 1):
		return 'prime'
	else:
		return 'composite'

# print(Wilson(101))