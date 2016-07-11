import numpy as np
import matplotlib.pyplot as plt


def file_to_array(quantity = 'I1', notebook_number = 7, sim_number = 21, version = ''):

    # assemble relative path
    rel_path = './etc/outputs/'
    if version != '':
         rel_path += 'DECSKS-' + str(version) + '_outputs/'

    # assemble filepath
    sim_name = 'out_' + quantity + '_' + 's' + str(notebook_number) + '-' + str(sim_number) + 'c'
    filepath = rel_path + sim_name
    infile = open(filepath, 'r')
    lines = infile.readlines()

    data = []
    for line in lines:
        data.append(eval(line.strip()))

    data_ndarray = np.zeros_like(data)

    for _ in range(len(data)):
        data_ndarray[_] = data[_]

    return data_ndarray

T = 2.
Nt = 500 # number of time steps, does not include time zero
t = np.linspace(0,T, Nt+1)

# I1

# s17-03
I1 = file_to_array(quantity = 'I1', notebook_number = 18, sim_number = 21)

I1 = (I1 - I1[0]) / I1[0]

plt.plot(t[1:], I1, 'Db', markersize = 12, linewidth = 2, label = 'DECSKS-2.2: s18-21c')


plt.xlabel('time', fontsize = 14)
plt.ylabel('Abs. value of fractional error in the $L^1$ norm', fontsize = 12)
plt.title('Fractional change in the $L^1$ norm', fontsize = 14)

#plt.axis([0,T,1e-16, 1e-13])
plt.legend(loc = 'best')
plt.grid()
plt.figure()

# I2

# s18-21c
I2 = file_to_array(quantity = 'I2', notebook_number = 18, sim_number = 21)
I2 = np.abs((I2 - I2[0]) / I2[0])
plt.semilogy(t[1:], I2, '-b', markersize = 12, linewidth = 2, label = 'DECSKS-2.2: s18-21c')

plt.xlabel('time', fontsize = 14)
plt.ylabel('Abs. value of fractional error in the $L^2$ norm', fontsize = 12)
plt.title('Fractional change in the $L^2$ norm', fontsize = 14)
#plt.axis([0,T,1e-16, 1e-5])
plt.legend(loc = 'best')
plt.grid()
plt.figure()

# IW

# s18-21c
IW = file_to_array(quantity = 'IW', notebook_number = 18, sim_number = 21)
IW = np.abs((IW - IW[0]) / IW[0])
plt.semilogy(t[1:], IW, '-b', markersize = 12, linewidth = 2, label = 'DECSKS-2.2: s18-21c')

#plt.axis([0,T,1e-16, 1e-5])
plt.xlabel('time', fontsize = 14)
plt.ylabel('Abs. value of fractional error in total energy, $W$', fontsize = 12)
plt.title('Fractional error in the total energy invariant $W$', fontsize = 14)
plt.legend(loc = 'best')
plt.grid()
plt.figure()


# S

# s18-21c
S = file_to_array(quantity = 'S', notebook_number = 18, sim_number = 21)
S = (S - S[0]) / S[0]
plt.semilogy(t[1:], S, '-b', markersize = 12, linewidth = 2, label = 'DECSKS-2.2: s18-21c')

plt.xlabel('time', fontsize = 14)
plt.ylabel('Abs. value of fractional error in the entropy $S$ invariant', fontsize = 12)
plt.title('Fractional change in the entropy invariant $S$', fontsize = 14)

#plt.axis([0,T,1e-16, 1e-4])
plt.legend(loc = 'best')
plt.grid()
plt.figure()


# WE

# s18-21c
WE = file_to_array(quantity = 'WE', notebook_number = 18, sim_number = 21)

plt.semilogy(t[1:], WE, '-b', markersize = 12, linewidth = 2, label = 'DECSKS-2.2: s18-21c')


plt.xlabel('time', fontsize = 14)
plt.ylabel('Normalized electrostatic energy, $W_E/W_{E0}$', fontsize = 12)
plt.title('DECSKS-2.2: Normalized electrostatic energy, $W_E / W_{E0}$')
#plt.axis([0,T,1e-16, 1e-4])

plt.legend(loc = 'best')
plt.grid()
