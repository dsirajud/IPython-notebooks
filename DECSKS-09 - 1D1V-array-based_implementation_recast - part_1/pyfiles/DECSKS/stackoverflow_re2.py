import numpy as np

N, M, L = 10, 100, 200

A = np.random.randn(N,M,M)
x = np.random.randn(M,L)

B = np.zeros((N,M,L))
for i in range(N):
    B[i,:,:] = np.dot(A[i,:,:], x)

def for_dot2(A,x):
    for i in range(A.shape[0]):
        np.dot(A[:,i,:], x)

def for_dot3(A,x):
    for i in range(A.shape[0]):
        np.dot(A[:,:,i], x)
