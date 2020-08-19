#######################################################################

# Primality Tests and Factorization Algorithms

# A library of primality tests and factorization algorithms.  Tests 
# and algorithms sourced from Cohen's "A Course in Computational 
# Algebraic Number Theory"

#######################################################################

# Usage instructions

# File Name
	# Method Name
		# Inputs
		# Outputs
		# Notes

# AKS.py
	# AKSTest(N)
		# Inputs
			# N, the number which will be tested for primality
		# Outputs
			# 'prime' if N is prime
			# 'composite' if N is composite
		# Notes
			# None

# CFRAC.py
	# CFRAC(N)
		# Inputs
			# N, the composite number which we are trying to factor
		# Outputs
			# a tuple of two nontrivial factors of N which multiply to 
				# give N
		# Notes
			# None

# ChineseRemainderTheorem.py
	# CRT_basic(xlist, mlist)
		# Inputs
			# xlist, list of values which the output should be congruent 
				# to
			# mlist, list of coprime moduli corresponding to each x in 
				# xlist
		# Outputs
			# the unique (mod M = product of all m in mlist) solution that 
				# satisfies all of the input congruences
			# None if mlist not coprime
		# Notes
			# None

	# CRT_precomp(mlist)
		# Inputs
			# mlist, list of coprime moduli
		# Outputs
			# Clist, a list of numbers z such that p * u + m_j * v = 1 
				# where p is the product of all m in mlist except m_j
			# None if mlist not coprime
		# Notes
			# for use with CRT_comp(xlist, mlist, Clist) and in 
				# CRT_multiple(xlistlist, mlist)

	# CRT_comp(xlist, mlist, Clist)
		# Inputs
			# xlist, list of values which the output should be congruent 
				# to
			# mlist, list of coprime moduli corresponding to each x in 
				# xlist
			# Clist, a list of numbers z such that p * u + m_j * v = 1 
				# where p is the product of all m in mlist except m_j
		# Outputs
			# the unique (mod M = product of all m in mlist) solution that 
				# satisfies all of the input congruences
			# None if mlist not coprime or lenghts of xlist and mlist are 
				# not equal
		# Notes
			# a single congruence computation using an advanced method

	# CRT_multiple(xlistlist, mlist)
		# Inputs
			# xlistlist, a list containing multiple xlists, which are each  
				# a list of values which the output should be congruent to
			# mlist, list of coprime moduli corresponding to each x in 
				# each xlist
		# Outputs
			# a list continaing the unique (mod M = product of all m in 
				# mlist) solution that satisfies all of the input 
				# congruences for each xlist in xlistlist
			# None if mlist not coprime or lenghts of xlist and mlist are 
				# not equal
		# Notes
			# multiple congruence computations using an advanced method

	# CRT_ind(xlist, mlist)
		# Inputs
			# xlistlist, a list containing multiple xlists, which are each  
				# a list of values which the output should be congruent to
			# mlist, list of coprime moduli corresponding to each x in 
				# each xlist
		# Outputs
			# the unique (mod M = product of all m in mlist) solution that 
				# satisfies all of the input congruences
			# None if mlist not coprime or lenghts of xlist and mlist are 
				# not equal
		# Notes
			# inductive version of the Chinese Remainder Theorem

# Cornacchia.py
	# Cornacchia(p, d)
		# Inputs
			# p, a prime
			# d, an integer 0 < d < p
		# Outputs
			# x, y, a solution to x**2 + d * y**2 = p if one exists
			# 'no solution' if such a solution does not exist
		# Notes

# ECM.py
	# ECM(N, dmax=128)
		# Inputs
			# N, the composite number which we are trying to factor
			# d, the maximum recursion depth per curve used, if not 
				# specified set to 128
		# Outputs
			# a tuple of two nontrivial factors of N which multiply to 
				# give N
		# Notes
			# None

# EuclideanAlgorithm.py
	# gcd(a, b)
		# Inputs
			# a, b, two integers
		# Outputs
			# the greatest common divisor of a, b
		# Notes
			# None

	# Bezout(a, b)
		# Inputs
			# a, b, two integers
		# Outputs
			# list of five numbers w, x, y, z, g = gcd(a, b) such that 
				# w * x + y * z = g
		# Notes
			# None

# Golwasser_Kilian.py
	# GKTest(N)
		# Inputs
			# N, the number we want to determine to be prime or 
				# composite
		# Outputs
			# 'prime' if N is prime
			# 'composite' if N is composite
		# Notes
			# N must pass Rabin-Miller test before GKTest is used, might 
				# run forever if N is composite

# Legendre.py
	# Kronecker(a, b)
		# Inputs
			# a, the integer on the top of the Kronecker symbol
			# b, the integer on the bottom of the Kronecker symbol
		# Outputs
			# 0 if a even and b = 2 or if a != += 1 and b = 0
			# 1 if a = +- 1 (mod 8) and b = 2 or if a >= 0 and b = -1 or 
				# if a = +- 1 and b = 0
			# - 1 if a = +- 3 (mod 8) and b = 2 or if a < 0 and b = -1
		# Notes
			# when b is an odd prime, equal to the Legendre symbol

# LehmansMethod.py
	# Lehman(N)
		# Inputs
			# N, the number which we are trying to factor or prove to be 
			# prime
		# Outputs
			# a factor of N if N is composite
			# 'prime' if N is prime
		# Notes
			# None

# PollardspMinus1.py
	# FirstStage(N, B, x=2)
		# Inputs
			# N, the number which we are trying to factor
			# B, a bound such that for any prime factor p that can be 
				# found, p - 1 is B-powersmooth (all primes in its prime 
				# factorization are less than or equal to B)
		# Outputs
			# a factor of N if the algorithm succeeds
			# 'the algorithm has failed' if the algorithm cannot factor 
				# N
			# 'the algorithm did not succeed in spliting N' if the bound 
				# B is too low to factor N
		# Notes
			# None

	# SecondStage(N, B1, B2)
		# Inputs
			# N, the number which we are trying to factor
			# B1 and B2, bounds such that for any prime factor p that can 
				# be found, p - 1 is equal to a B1-powersmooth number 
				# times a prime less than or equal to B2
		# Outputs
			# a factor of N if the algorithm succeeds
			# 'the algorithm has failed' if the algorithm cannot factor 
				# N
		# Notes
			# None

# PollardsRhoMethod.py
	# Pollard(N)
		# Inputs
			# N, a composite number which we are trying to factor
		# Outputs
			# a factor of N if the algorithm succeeds
			# 'algorithm fails' if the algorithm fails
		# Notes
			# None

# PoweringAlgorithms.py
	# RLbin(g, n)
		# Inputs
			# g, n, with which we want to calculate g**n
		# Outputs
			# g**n
		# Notes
			# None

	# calculate_e(n)
		# Inputs
			# n, a number that will be used as the exponent in 
				# LRbin(g, n, e)
		# Outputs
			# e such that 2**e <= |n| < 2**e+1
		# Notes
			# None

	# LRbin(g, n, e)
		# Inputs
			# g, n, with which we want to calculate g**n
			# e such that 2**e <= |n| < 2**e+1
		# Outputs
			# g**n
		# Notes
			# None

	# PowerModm(g, n, m)
		# Inputs
			# g, n, m, with which we want to calculate g**n mod m
		# Outputs
			# g**n mod m
		# Notes
			# None

# PowerTests.py
	# IntSqrt(N)
		# Inputs
			# N, the number we want the floor square root of
		# Outputs
			# the floor square root of N
		# Notes
			# None

	# precomp()
		# Inputs
			# None
		# Outputs
			# q11, q63, q64, q65, arrays to be used in 
				# SquareTest(N, q11, q63, q64, q65)
		# Notes
			# None

	# SquareTest(N, q11, q63, q64, q65)
		# Inputs
			# N, the number we want to determine to be a perfect square or not
			# q11, q63, q64, q65, arrays created by precomp()
		# Outputs
			# q11, q63, q64, q65, arrays t be used in 
				# SquareTest(N, q11, q63, q64, q65)
		# Notes
			# None

# PrimitiveRoot.py
	# PrimRoot(p)
		# Inputs
			# p, a prime number
		# Outputs
			# g, a number which to various powers gives each number between 0 
				# and p (noninclusive) mod p
		# Notes
			# None

# QS.py
	# QS(N)
		# Inputs
			# N, the composite number which we are trying to factor
		# Outputs
			# a tuple of two nontrivial factors of N which multiply to N
			# 'Cannot be solved with Qudratic Sieve, try Multiple Polynomial 
				# Quadratic Sieve' if algorithm fails
		# Notes
			# None

# Rabin_Miller.py
	# RM(N)
		# Inputs
			# N, the composite number which we are trying to factor
		# Outputs
			# 'composite' meaning that N is certainly composite
			# 'probably prime' meaning that N is probably prime, but could be 
				# composite
		# Notes
			# if RM(N) returns 'probably prime' there is only a 4**(-20) chance 
				# that it is composite

# SqrtModp.py
	# twokcalc(r)
		# Inputs
			# r, for which we want to find the largest integer t such that 
				# 2**t | r
		# Outputs
			# t, the largest integer such that 2**t | r
		# Notes
			# None

	# SqrtModp(a, p)
		# Inputs
			# a, a number we want the square root of mod p
			# p, the modulus for this equation
		# Outputs
			# x, such that x**2 = a mod p if it exists
			# str(a) + ' is not a quadratic residue mod ' + str(p) if it does 
				# not exist
		# Notes
			# None

# TrialDivisionFactoring.py
	# TrialFactor(N, B=lastelem)
		# Inputs
			# N, the number which we are trying to factor
			# B, an optional argument which gives the largest factor we want to 
				# test for
		# Outputs
			# a list containing one or two factors if there are factors less 
				# than B
			# 'prime' if the algorithm shows N is prime
			# 'remaining divisors are greater than ' + str(B) if the algorithm 
				# finds no factors less than or equal to B
		# Notes
			# finds one or two factors, not necessarily entire prime 
				# factorization

	# TrialCompleteFactorization(N, B=lastelem)
		# Inputs
			# N, the number which we are trying to factor
			# B, an optional argument which gives the largest factor we want to 
				# test for
		# Outputs
			# a list containing the prime factorization if all prime factors 
				# are less than B
			# 'prime' if the algorithm shows N is prime
			# a list containing all factors less than B and the message 
				# 'remaining factors of ' + str(N) + ' are greater than ' 
				# + str(B)
		# Notes
			# finds all factors less than B, not just one or two

# WilsonsTheorem.py
	# Wilson(p)
		# Inputs
			# p, a number we want to know to be prime or composite
		# Outputs	
			# 'prime' if p is prime
			# 'composite' if p is composite
		# Notes
			# None

#######################################################################