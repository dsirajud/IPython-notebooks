import numpy as np

class anything:
    def __init__(self, N, a = 0., b = 1.):

        self.N = N
        self.a = 0.
        self.b = 1.
        self.L = self.b - self.a
        
z = anything(N = 256)

xi = np.zeros(z.N)
for r in range(z.N):
    if r <= z.N/2 :
        xi[r] = 2*np.pi*r / z.L
    else:
        xi[r] = 2*np.pi*(r - z.N) / z.L
