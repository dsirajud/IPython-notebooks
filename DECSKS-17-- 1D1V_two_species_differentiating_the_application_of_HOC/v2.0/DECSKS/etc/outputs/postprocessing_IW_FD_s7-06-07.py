
#!/usr/bin/env python

import numpy as np
import pylab as plt

# PLOT param resets

params = {'legend.fontsize': 12,
          'legend.linewidth': 2}
plt.rcParams.update(params)

# PLOT TABLEAU

# plot tableau

tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),
         (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),
         (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),
         (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),
         (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]

# Scale the RGB values to the [0, 1] range, which is the format matplotlib accepts.
for i in range(len(tableau20)):
    r, g, b = tableau20[i]
    tableau20[i] = (r / 255., g / 255., b / 255.)


##### s13-06
infile = open('s13-06/out_IW','r')
lines = infile.readlines()

WE = []
for line in lines:
    WE.append(eval(line))

WE = np.array(WE)
WE = (WE - WE[0]) / WE[0] # normalize
WE = np.abs(WE)

Nt = 60
T = 60.

t = np.zeros(Nt+1) # includes t = 0
t = np.linspace(0,T, Nt+1)

plt.semilogy(t, WE, '--o', color = tableau20[6],linewidth = 3,
             label = 's13-06')


##### s13-07
infile = open('s13-07/out_IW','r')
lines = infile.readlines()

WE = []
for line in lines:
    WE.append(eval(line))

WE = np.array(WE)
WE = (WE - WE[0]) / WE[0] # normalize
WE = np.abs(WE)

Nt = 60
T = 60.

t = np.zeros(Nt+1) # includes t = 0
t = np.linspace(0,T, Nt+1)

plt.semilogy(t, WE, '--o', color = tableau20[2],linewidth = 3,
             label = 's13-07')


# decorate plot
plt.grid()
plt.xlabel('time', fontsize = 14)
plt.ylabel('Fractional error in total energy $I_W$', fontsize = 12)
plt.legend(loc = 'best')
plt.savefig('Landau_IW_FD_s7-07_and_s7-06.png')
plt.clf()
