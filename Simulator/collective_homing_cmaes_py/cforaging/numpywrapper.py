#!/usr/bin/env python

import numpy as np
from numpy import sqrt, real, diag, dot
from numpy.linalg import eig, norm
from numpy.random import randn

def seedwrapper(s):
	np.random.seed(s)

def randomwrapper():
	n = np.random.random()
	return n

def randnwrapper(rows,cols):
	M = randn(rows,cols)
	return M

def normwrapper(v):
	n = norm(v)
	return n

def eigwrapper(C):
	eigval, eigvec = eig(C)	   # eigen decomposition, B==normalized eigenvectors
	eigval = real(eigval)	   # enforce real value
	D = sqrt(eigval) # D contains standard deviations now
	eigvec = real(eigvec)        # enforce real value
	rows, cols = np.shape(C)
	size = (rows * cols)
	B = np.arange(size, dtype=np.float64)
	i = 0
	j = 0
	elem = 0
	while elem < size:
		B[elem] = eigvec[i,j]
		j += 1
		if (j == cols):
			i += 1
			j = 0
		elem += 1
	return D, B

def dotwrapper(a,b):
	# Get array ndims
	n_a = a.ndim
	n_b = b.ndim
	# Consistency check
	# Number of dimensions of arrays
	assert n_a >= 1 and n_a <= 2, "Invalid array <a>!!!"
	assert n_b >= 1 and n_b <= 2, "Invalid array <b>!!!"
	# Consistent number of cols (array <a>) and rows (array <b>)
	if n_a == 1 and n_b == 1:
		# If both <a> and <b> are one-dimensional arrays, their lengths must match
		assert len(a) == len(b), "Incosistent sizes for one-dimensional arrays!!!"
	else:
		if n_a == 1:
			# If <a> is a one-dimensional array, its length must be equal to number of rows of <b>
			b_r, b_c = np.shape(b)
			assert len(a) == b_r, "Inconsistent sizes: length of <a> do not match rows of <b>!!!"
		elif n_b == 1:
			# If <b> is a one-dimensional array, its length must be equal to number of cols of <a>
			a_r, a_c = np.shape(a)
			assert a_c == len(b), "Inconsistent sizes: cols of <a> do not match length of <b>!!!"
		else:
			# Number of cols of <a> must be equal to number of rows of <b>
			a_r, a_c = np.shape(a)
			b_r, b_c = np.shape(b)
			assert a_c == b_r, "Inconsistent sizes: cols of <a> do not match rows of <b>!!!"
	# Compute dot product
	c = dot(a,b)
	return c

