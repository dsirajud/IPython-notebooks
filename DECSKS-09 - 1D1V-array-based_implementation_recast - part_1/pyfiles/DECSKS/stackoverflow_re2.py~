
import numpy as np



dim1, dim2, dim3 = 20, 20, 20
A = np.random.rand(dim1, dim2, dim2)
x = np.random.rand(dim2, dim3)

def for_dot1(A,x):
    for i in range(A.shape[0]):
        np.dot(A[i,:,:], x)

def for_dot2(A,x):
    for i in range(A.shape[0]):
        np.dot(A[:,i,:], x)

def for_dot3(A,x):
    for i in range(A.shape[0]):
        np.dot(A[:,:,i], x)
