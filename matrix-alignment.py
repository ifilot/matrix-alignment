#!/usr/bin/env python

import numpy as np

# build two adjacency matrices
V1 = np.matrix('0 1 0 0 0 ; 1 0 1 0 0 ; 0 1 0 1 1 ; 0 0 1 0 0 ; 0 0 1 0 0')
V2 = np.matrix('0 0 1 0 0 ; 0 0 1 1 1 ; 1 1 0 0 0 ; 0 1 0 0 0 ; 0 1 0 0 0')

# define matrices
matsize = V1.shape[0]
Asize = matsize * matsize
A = np.zeros([Asize, Asize]) 		# stochastic transition matrix
E = np.zeros([Asize, 1])			# network alignment vector
b = np.ones([Asize,1])				# network overlap vector

for i in range(0, matsize):			# loop over V1
	for j in range(0, matsize):		# loop over V2
		if np.sum(V1[:,i]) == np.sum(V2[:,j]):
			E[i*5+j] = 1
E /= np.linalg.norm(E)

u = 0
for i in range(0, matsize):			# loop over V1
	for j in range(0, matsize):		# loop over V2
		v = 0
		for k in range(0, matsize):			# loop over V1
			for l in range(0, matsize):		# loop over V2
				if V1[i,k] == 1 and V2[j,l] == 1:
					A[u,v] = 1.0 / (np.sum(V1[:,k]) * np.sum(V2[:,l]))
				v += 1
		u += 1

# use both network information as well as "diffusion" information; use a power method
alpha = 0.6
for i in range(0,20):
	prod = alpha * np.dot(A,b) + (1 - alpha) * E
	b = prod / np.linalg.norm(prod)

# construct B in matrix form and output result
bmatrix = np.zeros([matsize, matsize])
for i in range(0,matsize):
	for j in range(0,matsize):
		bmatrix[i,j] = b[i * 5 + j]
print bmatrix

# construct index lists
rr = range(0, matsize)
cc = range(0, matsize)

# build empty permutation list
perm = np.zeros([matsize, 1])

# find all permutations
for i in range(0, matsize):		# find 5 permutations
	maxsize = -1
	ir = -1
	ic = -1
	for r in rr:
		for c in cc:
			if bmatrix[r,c] > maxsize:
				maxsize = bmatrix[r,c]
				ir = r
				ic = c
	# store index
	perm[ir] = ic
	# remove highest indices from list
	rr.remove(ir)
	cc.remove(ic)

print perm

# construct permutation matrix from these permutations
P = np.zeros([matsize, matsize])
for idx,val in enumerate(perm):
	P[idx,int(val)] = 1

# verify the alignment procedure
print V1
print P * V2 * P.transpose()