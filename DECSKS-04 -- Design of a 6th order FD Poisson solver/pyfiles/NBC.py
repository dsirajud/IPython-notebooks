import numpy as np
import numpy.linalg as LA

N = 10
a, b = 0.0, 1.0
L = float(b - a)
dx = L / (N - 1)

# assemble differencing matrix
D = np.zeros((N+1,N+1))
D[0,0] = -1/6.
D[0,1] = -77/60.
D[0,2] = 5/2.
D[0,3] = -5/3.
D[0,4] = 5/6.
D[0,5] = -1/4.
D[0,6] = 1/30.

for i in range(1,D.shape[0]-1):
    D[i,i-1] = 1
    D[i,i] = -2
    D[i,i+1] = 1


# Assemble FD coefficient matrix on n: B
B = np.zeros((N+1,N+1))
for i in range(B.shape[0]):
    if i == 1:
        B[i,i] = 318 / 240.
        B[i,i+1] = -133/120.
        B[i,i+2] = 187 / 120.
        B[i,i+3] = -23 / 20.
        B[i,i+4] = 109 / 240.
        B[i,i+5] = -3/40.

    elif i == 2:

        B[i, i-1] = 3 / 40.
        B[i, i] = 209 / 240.
        B[i,i+1] = 1 / 60.
        B[i,i+2] = 7 / 120.
        B[i,i+3] = -1 / 40.
        B[i,i+4] = 1 / 240.

    elif 3 <= i <= N-2:

        B[i,i-2] = -1/240.
        B[i,i-1] = 1/10.
        B[i,i] = 97/120.
        B[i,i+1] = 1/10.
        B[i,i+2] = -1/240.

    elif i == N-1:

        B[i,i+1] = 3 / 40.
        B[i,i] = 209 / 240.
        B[i,i-1] = 1 / 60.
        B[i,i-2] = 7 / 120.
        B[i,i-3] = -1 / 40.
        B[i,i-4] = 1 / 240.

    # else i == N: row of zeros

