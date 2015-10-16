import numpy as np
import numpy.random as npr
import time

dim1, dim2, dim3 = 10, 1000,1200
A = npr.randn(dim1, dim2, dim2)
x = npr.randn(dim2, dim3)
Ax = np.zeros([A.shape[0], A.shape[1], x.shape[1]]) 

def for_dot1(A,x):
    Ax[:,:,:] = np.dot(A,x)
    return Ax

def for_dot2(A,x):
    for i in range(A.shape[0]):
        Ax[i,:,:] = A[i,:,:].dot(x)
    return Ax

def for_dot3(A,x):
    Ax = np.einsum("ijk, kl -> ijl", A, x)

def for_dot4(A,x):
    Ax = np.tensordot(A,x, (2,0))

def for_dot5(A,x):
    Ax[:,:,:] = np.dot(A.reshape(-1,dim2), x).reshape(dim1, dim2, dim3)
    return Ax

