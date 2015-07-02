#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np

Nx = [24,
      48,
      96,
      192,
      384,
      768,
      1536]

L2 = [0.0338980800679,
      0.000220313520433,
      1.18516971559e-06,
      9.18647824875e-09,
      7.09906039822e-11,
      5.50658380527e-13,
      4.31540389839e-15]




print "Nx                L2 error                order"
print "--------------------------------------------------------"
for i in range(len(Nx)):
    if i == 0:
        print "%d                %2.5e" % (Nx[i], L2[i])
    else:
        order = np.log2(L2[i-1] / L2[i])
        print "%d                %2.5e                %g" % (Nx[i], L2[i], order)



order = -np.polyfit(np.log2(Nx[2:]), np.log2(L2[2:]), 1)[0]
fig, ax = plt.subplots()
ax.loglog(Nx, L2, '--o')
ax.set_xscale('log', basex=2)
ax.set_yscale('log', basey=2)
plt.ylabel('$L^2$ norm of $f_{num} - f_{exact}$', fontsize = 16)
plt.xlabel('number of gridpoints $N_x$', fontsize = 16)
plt.text(np.max(Nx), np.max(L2), r'order ($|slope|$) = %g' % order, horizontalalignment = 'right', verticalalignment = 'top', fontsize = 16)
plt.grid()
plt.show()
