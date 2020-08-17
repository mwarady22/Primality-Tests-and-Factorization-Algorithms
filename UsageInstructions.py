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

