#!/usr/bin/env python

import numpy as np

# build two adjacency matrices
V1 = np.matrix('0 1 0 0 0 ; 1 0 1 0 0 ; 0 1 0 1 1 ; 0 0 1 0 0 ; 0 0 1 0 0')
V2 = np.matrix('0 1 0 0 0 ; 1 0 1 0 0 ; 0 1 0 1 1 ; 0 0 1 0 0 ; 0 0 1 0 0')

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

alpha = 0.6
for i in range(0,20):
	prod = alpha * np.dot(A,b) + (1 - alpha) * E
	b = prod / np.linalg.norm(prod)

for i in range(0,matsize):
	for j in range(0,matsize):
		print b[i * 5 + j],
	print