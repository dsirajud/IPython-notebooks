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


# PLOT THE ELECTROSTATIC ENERGY IN A SYSTEM

##### s13-03
infile = open('s13-03/out_WE','r')
lines = infile.readlines()

WE = []
for line in lines:
    WE.append(eval(line))

WE = np.array(WE)
WE = WE / WE[0] # normalize

Nt = 60
dt = 1.0
T = 60.

t = np.zeros(Nt+1) # includes t = 0
t = np.linspace(0,T, Nt+1)

plt.semilogy(t, WE, '--o', color = tableau20[0], linewidth = 3,
             label = 's13-03')



##### s13-04
infile = open('s13-04/out_WE','r')
lines = infile.readlines()

WE = []
for line in lines:
    WE.append(eval(line))

WE = np.array(WE)
WE = WE / WE[0] # normalize

Nt = 200
T = 60.
dt = T / Nt

t = np.zeros(Nt+1) # includes t = 0
t = np.linspace(0,T, Nt+1)

plt.semilogy(t, WE, '--o', color = tableau20[4],linewidth = 3,
             label = 's13-04')

##### s13-05
infile = open('s13-05/out_WE','r')
lines = infile.readlines()

WE = []
for line in lines:
    WE.append(eval(line))

WE = np.array(WE)
WE = WE / WE[0] # normalize

Nt = 60
T = 60.

t = np.zeros(Nt+1) # includes t = 0
t = np.linspace(0,T, Nt+1)

plt.semilogy(t, WE, '--o', color = tableau20[5],linewidth = 3,
             label = 's13-05')



##### s13-06
infile = open('s13-06/out_WE','r')
lines = infile.readlines()

WE = []
for line in lines:
    WE.append(eval(line))

WE = np.array(WE)
WE = WE / WE[0] # normalize

Nt = 60
T = 60.

t = np.zeros(Nt+1) # includes t = 0
t = np.linspace(0,T, Nt+1)

plt.semilogy(t, WE, '--o', color = tableau20[6],linewidth = 3,
             label = 's13-06')


##### s13-06
infile = open('s13-07/out_WE','r')
lines = infile.readlines()

WE = []
for line in lines:
    WE.append(eval(line))

WE = np.array(WE)
WE = WE / WE[0] # normalize

Nt = 60
T = 60.

t = np.zeros(Nt+1) # includes t = 0
t = np.linspace(0,T, Nt+1)

plt.semilogy(t, WE, '--o', color = tableau20[0],linewidth = 3,
             label = 's13-07')


###### PLOT LEGACY DATA (v1.2), from simulations s7-##
infile = open('../../../../v1.2/DECSKS/etc/outputs/s7-03/out_WE','r')
lines = infile.readlines()

WE = []
for line in lines:
    WE.append(eval(line))

WE = np.array(WE)
WE = WE / WE[0] # normalize

Nt = 60
dt = 1.0
T = 60.

t = np.zeros(Nt+1) # includes t = 0
t = np.linspace(0,T, Nt+1)

plt.semilogy(t, WE, '--o', color = tableau20[7],linewidth = 3,
             label = 's7-03')



##### s7-04

infile = open('../../../../v1.2/DECSKS/etc/outputs/s7-03/out_WE','r')
lines = infile.readlines()

WE = []
for line in lines:
    WE.append(eval(line))

WE = np.array(WE)
WE = WE / WE[0] # normalize

Nt = 60
T = 60.

t = np.linspace(0,T, Nt+1)

plt.semilogy(t, WE, '--o', color = tableau20[11],linewidth = 3,
             label = 's7-04')



# LINEAR THEORY, WE = WE[0] exp(-2*gamma*t), if
#                E = E[0] exp(jw_r t) exp(-gamma t)
#                normalize: WE / WE[0] = exp(-gamma*t)

gamma = 0.153359 # damping constant

WE_linear = np.exp(-2*gamma*t)

plt.semilogy(t, WE_linear, '--m', linewidth = 3, label = 'Linear theory')

# decorate plot
plt.grid()
plt.xlabel('time', fontsize = 14)
plt.ylabel('Normalized electrostatic energy, $W_E/W_{E0}$', fontsize = 12)
plt.legend(loc = 'lower left')
plt.savefig('Landau_WE_FD.png')
plt.clf()
