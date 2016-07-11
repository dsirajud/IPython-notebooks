import numpy as np
import matplotlib.pyplot as plt
import numpy as np

def file_to_array(quantity = 'I1', notebook_number = 7, sim_number = 20, sim_number_version = 'b', version = ''):

    # assemble relative path
    rel_path = './etc/outputs/'
    if version != '':
         rel_path += 'DECSKS-' + str(version) + '_outputs/'

    # assemble filepath
    sim_name = 'out_' + quantity + '_' + 's' + str(notebook_number) + '-' + str(sim_number) + sim_number_version
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

T = 50.
Nt = 100 # number of time steps, does not include time zero
t = np.linspace(0,T, Nt+1)

# I1

I1 = file_to_array(quantity = 'I1', notebook_number = 18, sim_number = 20, sim_number_version = '')
I1 = np.abs((I1 - I1[0]) / I1[0])
plt.semilogy(t, I1, 'Db', markersize = 12, linewidth = 2, label = r'DECSKS-2.2: s18-20, $v_{th} = 0.5$')

I1 = file_to_array(quantity = 'I1', notebook_number = 18, sim_number = 20, sim_number_version = 'b')
I1 = np.abs((I1 - I1[0]) / I1[0])
plt.semilogy(t, I1, 'Dc', markersize = 12, linewidth = 2, label = r'DECSKS-2.2: s18-20b, $v_{th} = 0.0625$')

I1 = file_to_array(quantity = 'I1', notebook_number = 18, sim_number = 20, sim_number_version = 'c')
I1 = np.abs((I1 - I1[0]) / I1[0])
plt.semilogy(t, I1, 'Dm', markersize = 12, linewidth = 2, label = r'DECSKS-2.2: s18-20c, $v_{th} = 0.03125$')


plt.xlabel('time', fontsize = 14)
plt.ylabel('Abs. value of fractional error in the $L^1$ norm', fontsize = 12)
plt.title('Fractional error in the $L^1$ norm invariant (mass conservation)', fontsize = 14)

plt.axis([0,T,1e-16, 1e-13])
plt.legend(loc = 'best')
plt.grid()
plt.figure()

# I2

# s18-20
I2 = file_to_array(quantity = 'I2', notebook_number = 18, sim_number = 20, sim_number_version = '')
I2 = np.abs((I2 - I2[0]) / I2[0])
plt.semilogy(t, I2, '-b', markersize = 12, linewidth = 2, label = r's18-20, $v_{th} = 0.5$')

I2 = file_to_array(quantity = 'I2', notebook_number = 18, sim_number = 20, sim_number_version = 'b')
I2 = np.abs((I2 - I2[0]) / I2[0])
plt.semilogy(t, I2, '-c', markersize = 12, linewidth = 2, label = r's18-20b, $v_{th} = 0.0625$')

I2 = file_to_array(quantity = 'I2', notebook_number = 18, sim_number = 20, sim_number_version = 'c')
I2 = np.abs((I2 - I2[0]) / I2[0])
plt.semilogy(t, I2, '-m', markersize = 12, linewidth = 2, label = r's18-20c, $v_{th} = 0.03125$')


plt.xlabel('time', fontsize = 14)
plt.ylabel('Abs. value of fractional error in the $L^2$ norm', fontsize = 12)
plt.title('Fractional error in the $L^2$ norm invariant', fontsize = 14)
#plt.axis([0,T,1e-16, 1e-5])
plt.legend(loc = 'best')
plt.grid()
plt.figure()

# IW

# s18-20
IW = file_to_array(quantity = 'IW', notebook_number = 18, sim_number = 20, sim_number_version = '')
IW = np.abs((IW - IW[0]) / IW[0])
plt.semilogy(t, IW, '-b', markersize = 12, linewidth = 2, label = r'DECSKS-2.2: s18-20, $v_{th} = 0.5$')

# s18-20
IW = file_to_array(quantity = 'IW', notebook_number = 18, sim_number = 20, sim_number_version = 'b')
IW = np.abs((IW - IW[0]) / IW[0])
plt.semilogy(t, IW, '-c', markersize = 12, linewidth = 2, label = r'DECSKS-2.2: s18-20b, $v_{th} = 0.0625$')

# s18-20
IW = file_to_array(quantity = 'IW', notebook_number = 18, sim_number = 20, sim_number_version = 'c')
IW = np.abs((IW - IW[0]) / IW[0])
plt.semilogy(t, IW, '-m', markersize = 12, linewidth = 2, label = r'DECSKS-2.2: s18-20c, $v_{th} = 0.03125$')


#plt.axis([0,T,1e-16, 1e-5])
plt.xlabel('time', fontsize = 14)
plt.ylabel('Abs. value of fractional error in total energy, $W$', fontsize = 12)
plt.title('Fractional error in the total energy invariant $W$', fontsize = 14)
plt.legend(loc = 'best')
plt.grid()
plt.figure()


# S

# s18-20
S = file_to_array(quantity = 'S', notebook_number = 18, sim_number = 20, sim_number_version = '')
S = np.abs((S - S[0]) / S[0])
plt.semilogy(t, S, '-b', markersize = 12, linewidth = 2, label = r's18-20, $v_{th} = 0.5$')

S = file_to_array(quantity = 'S', notebook_number = 18, sim_number = 20, sim_number_version = 'b')
S = np.abs((S - S[0]) / S[0])
plt.semilogy(t, S, '-c', markersize = 12, linewidth = 2, label = r's18-20b, $v_{th} = 0.0625$')

S = file_to_array(quantity = 'S', notebook_number = 18, sim_number = 20, sim_number_version = 'c')
S = np.abs((S - S[0]) / S[0])
plt.semilogy(t, S, '-m', markersize = 12, linewidth = 2, label = r's18-20c, $v_{th} = 0.03125$')


plt.xlabel('time', fontsize = 14)
plt.ylabel('Abs. value of fractional error in the entropy $S$ invariant', fontsize = 12)
plt.title('Fractional error in the entropy invariant $S$', fontsize = 14)

#plt.axis([0,T,1e-16, 1e-4])
plt.legend(loc = 'best')
plt.grid()
plt.figure()


# WE

WE = file_to_array(quantity = 'WE', notebook_number = 18, sim_number = 20, sim_number_version = '')
WE = WE / WE[0]
plt.semilogy(t, WE, '--b', markersize = 12, linewidth = 2, label = r'DECSKS-2.2: s18-20, $v_{th} = 0.5$')

WE = file_to_array(quantity = 'WE', notebook_number = 18, sim_number = 20, sim_number_version = 'b')
WE = WE / WE[0]
plt.semilogy(t, WE, '--c', markersize = 12, linewidth = 2, label = r'DECSKS-2.2: s18-20b, $v_{th} = 0.0625$')

WE = file_to_array(quantity = 'WE', notebook_number = 18, sim_number = 20, sim_number_version = 'c')
WE = WE / WE[0]
plt.semilogy(t, WE, '--m', markersize = 12, linewidth = 2, label = r'DECSKS-2.2: s18-20c, $v_{th} = 0.03125$')

# Linear wave growth estimate, gamma = 1 / sqrt(8)
gamma = 1 / np.sqrt(8) # from linear theory for this distribution function (s18-20)
WE_linear = np.exp(2*gamma*t - 5.) # Linear Landau theory prediction for normalized W_E / W_E0
plt.semilogy(t, WE_linear, '-k', linewidth = 3, label = 'linear theory')

plt.xlabel('time', fontsize = 14)
plt.ylabel('Normalized electrostatic energy, $W_E/W_{E0}$', fontsize = 12)
plt.title('DECSKS-2.2: Normalized electrostatic energy, $W_E / W_{E0}$')
plt.axis([0,T,1e-5, 1e7])

plt.legend(loc = 'best')
plt.grid()
