from sage.all import *
import sys

#######################################################################

# Usage instructions for all other files

# File Name
	# Method Name
		# Inputs
		# Outputs
		# Notes

#######################################################################

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