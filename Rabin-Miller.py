from sage.all import *
import sys
import random
from PoweringAlgorithms import *

def RM(N):
	if N % 2 == 0: # check if even
		return 'composite'
	return initialize(N) # branch to step 1


def initialize(N): # step 1
	q = N - 1
	t = 0
	while q % 2 == 0:
		q = int(q / 2)
		t += 1
	c = 20
	return new_a(N, q, t, c) # branch to step 2

def new_a(N, q, t, c): # step 2
	a = randint(1, N - 1)
	e = 0
	b = LRbin(a, q, calculate_e(q))
	b = b % N
	if b == 1:
		return repeat_test(N, q, t, c, b) # branch to step 4
	return squarings(N, q, t, c, b, e)


def squarings(N, q, t, c, b, e): # step 3
	while ((b % N) != 1) and ((b % N) != N - 1) and (e <= t - 2):
		b = (b**2) % N
		e += 1
	if (b != N - 1):
		return 'composite'
	return repeat_test(N, q, t, c, b) # branch to step 4


def repeat_test(N, q, t, c, b): # step 4
	c -= 1
	if c < 0:
		return new_a(N, q, t, c) # branch to step 2
	else:
		return 'probably prime'

# print(RM(2047))