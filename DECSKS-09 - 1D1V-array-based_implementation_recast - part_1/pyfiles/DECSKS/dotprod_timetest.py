import numpy as np
import numpy.random as npr
import time

dim1, dim2, dim3 = 10, 1000, 1200
A = npr.randn(dim1, dim2, dim2)
x = npr.randn(dim2, dim3)

t = time.time()
Ax1 = np.zeros([A.shape[0], x.shape[0], x.shape[1]]) 
Ax1[:,:,:] = np.dot(A,x)
td1 = time.time() - t
print "a single dot product of A [shape = (%d, %d, %d)] with x [shape = (%d, %d)] completes in %.3f s" \
  % (A.shape[0], A.shape[1], A.shape[2], x.shape[0], x.shape[1], td1)


Ax2 = np.zeros([A.shape[0], x.shape[0], x.shape[1]])
t = time.time()
for i in range(A.shape[0]):
    Ax2[i,:,:] = A[i,:,:].dot(x)
td2 = time.time() - t
print "taking %d dot products of 2D dot products A[i,:,:] [shape = (%d, %d)] with x [shape = (%d, %d)] completes in %.3f s" \
  % (A.shape[0], A.shape[1], A.shape[2], x.shape[0], x.shape[1], td2)

t = time.time()
Ax3 = np.einsum("ijk, kl -> ijl", A, x)
td3 = time.time() - t
print "using np.einsum, it completes in %.3f s" % td3


print "validate computation, is the difference less than 100*machine_eps?"
eps = np.finfo(np.float).eps
TOL = 100*eps
print np.all(np.abs(Ax1 - Ax2) < TOL)
print np.all(np.abs(Ax1 - Ax3) < TOL)
print np.all(np.abs(Ax2 - Ax3) < 10*TOL)
print Ax2.shape
