#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np

Nx = [9,
      18,
      36,
      72,
      144,
      288,
      576]



L2 = [0.161384101359,
      0.0665531537064,
      0.000896903347292,
      7.15887167286e-07,
      1.32551033016e-09,
      2.54590670652e-12,
      1.25901665399e-14]


 
order = -np.polyfit(np.log2(Nx[1:-1]), np.log2(L2[1:-1]), 1)[0]

print "Nx                L2 error                order"
print "--------------------------------------------------------"
for i in range(len(Nx)):
    if i == 0:
        print "%d                %2.5e" % (Nx[i], L2[i])
    else:
        order = np.log2(L2[i-1] / L2[i])
        print "%d                %2.5e                %g" % (Nx[i], L2[i], order)


fig, ax = plt.subplots()
ax.loglog(Nx, L2, '--o')
ax.set_xscale('log', basex=2)
ax.set_yscale('log', basey=2)
plt.ylabel('$L^2$ norm of $f_{num} - f_{exact}$', fontsize = 16)
plt.xlabel('number of gridpoints $N_x$', fontsize = 16)
plt.text(np.max(Nx), np.max(L2), r'order ($|slope|$) = %g' % order, horizontalalignment = 'right', verticalalignment = 'top', fontsize = 16)
plt.grid()
plt.show()
