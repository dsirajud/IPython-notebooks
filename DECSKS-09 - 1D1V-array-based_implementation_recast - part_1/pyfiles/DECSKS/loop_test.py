import numpy as np

def setup():
    A = np.random.randn(200,200)
    B = np.random.randn(200,200)

    return A,B


def main(A,B):

    C = np.zeros(A.shape)
    for j in range(A.shape[1]):
        for i in range(B.shape[0]):
            for k in range(A.shape[1]):
                C[i,j] += A[i,k] * B[k,j]

    return C
