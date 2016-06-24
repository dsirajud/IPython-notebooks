import numpy as np
import matplotlib.pyplot as plt


def file_to_array(quantity = 'I1', notebook_number = 7, sim_number = 3, version = ''):

    # assemble relative path
    rel_path = './etc/outputs/'
    if version != '':
         rel_path += 'DECSKS-' + str(version) + '_outputs/'

    # assemble filepath
    sim_name = 'out_' + quantity + '_' + 's' + str(notebook_number) + '-0' + str(sim_number)
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

T = 1000.
Nt = 2000 # number of time steps, does not include time zero
t = np.linspace(0,T, Nt+1)

# I1

# s17-03
I1 = file_to_array(quantity = 'I1', notebook_number = 17, sim_number = 3)
I1 = np.abs((I1 - I1[0]) / I1[0])
plt.semilogy(t, I1, '-b', markersize = 12, linewidth = 2, label = 'DECSKS-2.2: s17-03')



# s17-04
I1 = file_to_array(quantity = 'I1', notebook_number = 17, sim_number = 4)
I1 = np.abs((I1 - I1[0]) / I1[0])
plt.semilogy(t, I1, '--c', markersize = 12, linewidth = 2, label = 'DECSKS-2.2: s17-04')



plt.xlabel('time', fontsize = 14)
plt.ylabel('Abs. value of fractional error in the $L^1$ norm', fontsize = 12)
plt.title('Fractional error in the $L^1$ norm invariant (mass conservation)', fontsize = 14)

plt.axis([0,T,1e-16, 1e-13])
plt.legend(loc = 'best')
plt.grid()
plt.figure()

# I2

# s17-03
I2 = file_to_array(quantity = 'I2', notebook_number = 17, sim_number = 3)
I2 = np.abs((I2 - I2[0]) / I2[0])
plt.semilogy(t, I2, '-b', markersize = 12, linewidth = 2, label = 'DECSKS-2.2: s17-03')




# s17-04
I2 = file_to_array(quantity = 'I2', notebook_number = 17, sim_number = 4)
I2 = np.abs((I2 - I2[0]) / I2[0])
plt.semilogy(t, I2, '--c', markersize = 12, linewidth = 2, label = 'DECSKS-2.2: s17-04')


plt.xlabel('time', fontsize = 14)
plt.ylabel('Abs. value of fractional error in the $L^2$ norm', fontsize = 12)
plt.title('Fractional error in the $L^2$ norm invariant', fontsize = 14)
#plt.axis([0,T,1e-16, 1e-5])
plt.legend(loc = 'lower right')
plt.grid()
plt.figure()


# IW

# s17-03
IW = file_to_array(quantity = 'IW', notebook_number = 17, sim_number = 3)
IW = np.abs((IW - IW[0]) / IW[0])
plt.semilogy(t, IW, '-b', markersize = 12, linewidth = 2, label = 'DECSKS-2.2: s17-03')

# s17-04
IW = file_to_array(quantity = 'IW', notebook_number = 17, sim_number = 4)
IW = np.abs((IW - IW[0]) / IW[0])
plt.semilogy(t, IW, '--c', markersize = 12, linewidth = 2, label = 'DECSKS-2.2: s17-04')

#plt.axis([0,T,1e-16, 1e-5])
plt.xlabel('time', fontsize = 14)
plt.ylabel('Abs. value of fractional error in total energy, $W$', fontsize = 12)
plt.title('Fractional error in the total energy invariant $W$', fontsize = 14)
plt.legend(loc = 'lower right')
plt.grid()
plt.figure()


# S

# s17-03
S = file_to_array(quantity = 'S', notebook_number = 17, sim_number = 3)
S = np.abs((S - S[0]) / S[0])
plt.semilogy(t, S, '-b', markersize = 12, linewidth = 2, label = 'DECSKS-2.2: s17-03')

# s17-04
S = file_to_array(quantity = 'S', notebook_number = 17, sim_number = 4)
S = np.abs((S - S[0]) / S[0])
plt.semilogy(t, S, '--c', markersize = 12, linewidth = 2, label = 'DECSKS-2.2: s17-04')


plt.xlabel('time', fontsize = 14)
plt.ylabel('Abs. value of fractional error in the entropy $S$ invariant', fontsize = 12)
plt.title('Fractional error in the entropy invariant $S$', fontsize = 14)

#plt.axis([0,T,1e-16, 1e-4])
plt.legend(loc = 'lower right')
plt.grid()
plt.figure()


# WE

# s17-03
WE = file_to_array(quantity = 'WE', notebook_number = 17, sim_number = 3)
WE = WE / WE[0]
plt.semilogy(t, WE, '-b', markersize = 12, linewidth = 2, label = 'DECSKS-2.2: s17-03')

# s17-04
WE = file_to_array(quantity = 'WE', notebook_number = 17, sim_number = 4)
WE = WE / WE[0]
plt.semilogy(t, WE, '--c', markersize = 12, linewidth = 2, label = 'DECSKS-2.2: s17-04')



plt.xlabel('time', fontsize = 14)
plt.ylabel('Normalized electrostatic energy, $W_E/W_{E0}$', fontsize = 12)
plt.title('DECSKS-2.2: Normalized electrostatic energy, $W_E / W_{E0}$')
#plt.axis([0,T,1e-16, 1e-4])

plt.legend(loc = 'lower right')
plt.grid()
