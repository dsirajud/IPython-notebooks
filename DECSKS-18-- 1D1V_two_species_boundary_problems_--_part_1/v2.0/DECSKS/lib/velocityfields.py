import numpy as np

class Velocity:
    """Generates constant velocity fields for various test problems

    inputs:
    field -- (str) type of velocity field to construct
    z1 -- (instance) phase space variable
    z2 -- (instance) phase space variable
    """
    def __init__(self, z1, z2 = None, field = 'const'):
        if z2 is None:
            self.N = z1.N
            self.gridpoints = np.array(range(self.N))

        else:
            self.N1 = z1.N
            self.N2 = z2.N

        if field.lower() == 'rotating':
            self.z1values = -self.rotating(z2)
            self.z2values = self.rotating(z1)

        else:
            self.cells = eval('self.' + field + '(z1.cells)')

    def rotating(self, z):
        return 2*np.pi*z.cells

    def linear(self, x):
        v_low = -0.25
        v_high = 0.75
        return np.where(x < 0.0, (v_low + 1*x/3.), (v_high - 2*x/3.))

    def step(self, x):
        v_low = -1
        v_high = 2
        return np.where(x < 0, v_low,v_high)

    def LCS(self, x):
        v = np.zeros(len(x))
        v_low = -1
        T      = 0.4
        for each in range(len(x)):
            if x[each] < 0:
                v[each] = v_low
            elif 0 <= x[each] < 0.10:
                v[each] = 1 + 2*x[each]/1.
            else:
                v[each] = 1.5 - 1.*np.sin(2*np.pi*x[each]/T)
        return v

    def sinusoid(self, x):
        return np.sin(2*np.pi*x / 1.0)

    def const(self, x):
        v_const = 2.1
        return v_const*np.ones(len(x))

 
