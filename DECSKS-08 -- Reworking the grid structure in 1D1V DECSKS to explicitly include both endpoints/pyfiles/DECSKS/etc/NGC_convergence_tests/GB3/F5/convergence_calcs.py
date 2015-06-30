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

L2 = [0.062051276717,
      0.00659756512187,
      0.000297852626006,
      6.22854832121e-06,
      1.02530277211e-07,
      1.60778470232e-09,
      2.50506401429e-11]

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
