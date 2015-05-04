from coeffs import *
import numpy as np
import math
from periodic_matrix import *
import numpy as np
import matplotlib.pyplot as plt
import numpy

NumGrids = 10
Nx = 18

#Nmin = 3
#Nmax = 5 + 1
error_norm = np.zeros(NumGrids)
gridpoint = 0
p = 1
dn = 2
dn_max = 12
p = 12 # LTE on derivatives
abscissa = []
order = []

for q in range(NumGrids):
    print q
    error_norm[q] = main(dn, Nx, dn_max, p)

    #    if q > Nmin:
        #        order = math.log( error_norm[q-1] / error_norm[q], 2)
        #        print "Nx%d/Nx%d , order = %g" % (2**(q-1), 2 ** q, order)
        #        pass
    abscissa.append(Nx)
    Nx *= 2

abscissa = np.array(abscissa)
h = 2./abscissa
order = -numpy.polyfit(np.log2(abscissa), np.log2(error_norm), 1)[0]
print "error = "
print error_norm

print "spacings = "
print h

print "Nx = "
print abscissa


print "order = %g" % order


orders = np.zeros(NumGrids)

for s in range(1,NumGrids):

    orders[s] = np.log2(error_norm[s-1] / error_norm[s])

print "orders on each step"
print orders

offset = np.max(error_norm)
order_line = offset + 2 ** (-p) * np.ones(len(abscissa))

fig, ax = plt.subplots()
ax.loglog(abscissa, error_norm, 'o')
#ax.hold('on')
#ax.loglog(abscissa, order_line, '-b')
ax.set_xscale('log', basex=2)
ax.set_yscale('log', basey=2)
#ax.set_xlim(1, 1e4)
#ax.set_ylim(1e-15, 1)
plt.show()
