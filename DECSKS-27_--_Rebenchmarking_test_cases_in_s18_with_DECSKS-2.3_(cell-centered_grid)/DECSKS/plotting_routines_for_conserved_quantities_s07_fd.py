import numpy as np
import matplotlib.pyplot as plt


def file_to_array(quantity = 'I1', notebook_number = 7, sim_number = 1, version = ''):

    # assemble relative path
    rel_path = './etc/outputs/'
    if version != '':
         rel_path += 'DECSKS-' + str(version) + '_outputs/'

    # assemble filepath
    sim_name = 'out_' + quantity + '_' + 's0' + str(notebook_number) + '-0' + str(sim_number)
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


T = 60.
Nt_03 = 60 # number of time steps, does not include time zero
t_03 = np.linspace(0,T, Nt_03+1)
t = t_03.copy()

Nt_04 = 200
t_04 = np.linspace(0,T, Nt_04+1)

Nt_05 = 60
t_05 = np.linspace(0,T, Nt_05+1)

# I1

# s07-03
I1 = file_to_array(quantity = 'I1', notebook_number = 7, sim_number = 3)

I1 = np.abs((I1 - I1[0]) / I1[0])
print I1
plt.semilogy(t, I1, 'oc', markersize = 12, linewidth = 2, label = 'DECSKS-2.2: s7-03')

# s07-04
I1 = file_to_array(quantity = 'I1', notebook_number = 7, sim_number = 4, version = 2.2)
I1 = np.abs((I1 - I1[0]) / I1[0])
plt.semilogy(t_04, I1, '-r', linewidth = 2, label = 'DECSKS-2.2: s7-04')

# s07-05
I1 = file_to_array(quantity = 'I1', notebook_number = 7, sim_number = 5)
I1 = np.abs((I1 - I1[0]) / I1[0])
plt.semilogy(t, I1, 'Db', linewidth = 2, label = 'DECSKS-2.2: s7-05')


# DECSKS-1.2 data
# s07-03
I1 = file_to_array(quantity = 'I1', notebook_number = 7, sim_number = 3, version = 1.2)

I1 = np.abs((I1 - I1[0]) / I1[0])
plt.semilogy(t, I1, '*m', markersize = 10, linewidth = 2, label = 'DECSKS-1.2: s7-03')

# s07-04
I1 = file_to_array(quantity = 'I1', notebook_number = 7, sim_number = 4, version = 1.2)
I1 = np.abs((I1 - I1[0]) / I1[0])
plt.semilogy(t_04, I1, '-g', linewidth = 2, label = 'DECSKS-1.2: s7-04')

# s07-05
I1 = file_to_array(quantity = 'I1', notebook_number = 7, sim_number = 5, version = 1.2)
I1 = np.abs((I1 - I1[0]) / I1[0])
plt.semilogy(t, I1, 'Hy', linewidth = 2, label = 'DECSKS-1.2: s7-05')



plt.xlabel('time', fontsize = 14)
plt.ylabel('Abs. value of fractional error in the $L^1$ norm', fontsize = 12)
plt.title('Fractional error in the $L^1$ norm invariant (mass conservation)', fontsize = 14)

plt.axis([0,60.,1e-16, 1e-6])
plt.legend(loc = 'best')
plt.grid()
plt.figure()

# I2

# s07-03
I2 = file_to_array(quantity = 'I2', notebook_number = 7, sim_number = 3)
I2 = np.abs((I2 - I2[0]) / I2[0])
plt.semilogy(t, I2, 'oc', markersize = 12, linewidth = 2, label = 'DECSKS-2.2: s7-03')

# s07-04
I2 = file_to_array(quantity = 'I2', notebook_number = 7, sim_number = 4, version = 2.2)
I2 = np.abs((I2 - I2[0]) / I2[0])
plt.semilogy(t_04, I2, '-r', linewidth = 2, label = 'DECSKS-2.2: s7-04')

# s07-05
I2 = file_to_array(quantity = 'I2', notebook_number = 7, sim_number = 5)
I2 = np.abs((I2 - I2[0]) / I2[0])
plt.semilogy(t, I2, 'Db', linewidth = 2, label = 'DECSKS-2.2: s7-05')



# DECSKS-1.2 data
# s07-03
I2 = file_to_array(quantity = 'I2', notebook_number = 7, sim_number = 3, version = 1.2)

I2 = np.abs((I2 - I2[0]) / I2[0])
plt.semilogy(t, I2, '*m', markersize = 10, linewidth = 2, label = 'DECSKS-1.2: s7-03')

# s07-04
I2 = file_to_array(quantity = 'I2', notebook_number = 7, sim_number = 4, version = 1.2)
I2 = np.abs((I2 - I2[0]) / I2[0])
plt.semilogy(t_04, I2, '-g', linewidth = 2, label = 'DECSKS-1.2: s7-04')

# s07-05
I2 = file_to_array(quantity = 'I2', notebook_number = 7, sim_number = 5, version = 1.2)
I2 = np.abs((I2 - I2[0]) / I2[0])
plt.semilogy(t, I2, 'Hy', linewidth = 2, label = 'DECSKS-1.2: s7-05')


plt.xlabel('time', fontsize = 14)
plt.ylabel('Abs. value of fractional error in the $L^2$ norm', fontsize = 12)
plt.title('Fractional error in the $L^2$ norm invariant', fontsize = 14)
plt.axis([0,60.,1e-16, 1e-1])
plt.legend(loc = 'lower right')
plt.grid()
plt.figure()


# IW

# s07-03
IW = file_to_array(quantity = 'IW', notebook_number = 7, sim_number = 3)
IW = np.abs((IW - IW[0]) / IW[0])
plt.semilogy(t, IW, 'oc', markersize = 12, linewidth = 2, label = 'DECSKS-2.2: s7-03')

# s07-04
IW = file_to_array(quantity = 'IW', notebook_number = 7, sim_number = 4, version = 2.2)
IW = np.abs((IW - IW[0]) / IW[0])
plt.semilogy(t_04, IW, '-r', linewidth = 2, label = 'DECSKS-2.2: s7-04')

# s07-05
IW = file_to_array(quantity = 'IW', notebook_number = 7, sim_number = 5)
IW = np.abs((IW - IW[0]) / IW[0])
plt.semilogy(t, IW, 'Db', linewidth = 2, label = 'DECSKS-2.2: s7-05')


# DECSKS-1.2 data
# s07-03
IW = file_to_array(quantity = 'IW', notebook_number = 7, sim_number = 3, version = 1.2)

IW = np.abs((IW - IW[0]) / IW[0])
plt.semilogy(t, IW, '*m', markersize = 10, linewidth = 2, label = 'DECSKS-1.2: s7-03')

# s07-04
IW = file_to_array(quantity = 'IW', notebook_number = 7, sim_number = 4, version = 1.2)
IW = np.abs((IW - IW[0]) / IW[0])
plt.semilogy(t_04, IW, '-g', linewidth = 2, label = 'DECSKS-1.2: s7-04')

# s07-05
IW = file_to_array(quantity = 'IW', notebook_number = 7, sim_number = 5, version = 1.2)
IW = np.abs((IW - IW[0]) / IW[0])
plt.semilogy(t, IW, 'Hy', linewidth = 2, label = 'DECSKS-1.2: s7-05')



plt.axis([0,60.,1e-16, 1e-1])
plt.xlabel('time', fontsize = 14)
plt.ylabel('Abs. value of fractional error in total energy, $W$', fontsize = 12)
plt.title('Fractional error in the total energy invariant $W$', fontsize = 14)
plt.legend(loc = 'best')
plt.grid()
plt.figure()

# S

# s07-03
S = file_to_array(quantity = 'S', notebook_number = 7, sim_number = 3)
S = np.abs((S - S[0]) / S[0])
plt.semilogy(t, S, 'oc', markersize = 12, linewidth = 2, label = 'DECSKS-2.2: s7-03')

# s07-04
S = file_to_array(quantity = 'S', notebook_number = 7, sim_number = 4, version = 2.2)
S = np.abs((S - S[0]) / S[0])
plt.semilogy(t_04, S, '-r', linewidth = 2, label = 'DECSKS2.2: s7-04')

# s07-05
S = file_to_array(quantity = 'S', notebook_number = 7, sim_number = 5)
S = np.abs((S - S[0]) / S[0])
plt.semilogy(t, S, 'Db', linewidth = 2, label = 'DECSKS-2.2: s7-05')


# DECSKS-1.2 data
# s07-03
S = file_to_array(quantity = 'S', notebook_number = 7, sim_number = 3, version = 1.2)

S = np.abs((S - S[0]) / S[0])
plt.semilogy(t, S, '*m', markersize = 10, linewidth = 2, label = 'DECSKS-1.2: s7-03')

# s07-04
S = file_to_array(quantity = 'S', notebook_number = 7, sim_number = 4, version = 1.2)
S = np.abs((S - S[0]) / S[0])
plt.semilogy(t_04, S, '-g', linewidth = 2, label = 'DECSKS-1.2: s7-04')


# s07-05
S = file_to_array(quantity = 'S', notebook_number = 7, sim_number = 5, version = 1.2)
S = np.abs((S - S[0]) / S[0])
plt.semilogy(t, S, 'Hy', linewidth = 2, label = 'DECSKS-1.2: s7-05')

plt.xlabel('time', fontsize = 14)
plt.ylabel('Abs. value of fractional error in the entropy $S$ invariant', fontsize = 12)
plt.title('Fractional error in the entropy invariant $S$', fontsize = 14)

plt.axis([0,60.,1e-16, 1e-1])
plt.legend(loc = 'best')
plt.grid()
plt.figure()

# WE

# s07-03
#WE = file_to_array(quantity = 'WE', notebook_number = 7, sim_number = 3)
#WE = WE / WE[0]
#plt.semilogy(t, WE, 'oc', markersize = 12, linewidth = 2, label = 'DECSKS-2.2: s7-03')


# s07-04
#WE = file_to_array(quantity = 'WE', notebook_number = 7, sim_number = 4, version = 2.2)
#WE = WE / WE[0]
#plt.semilogy(t_04, WE, '-r', linewidth = 2, label = 'DECSKS-2.2: s7-04')


# s07-05
WE = file_to_array(quantity = 'WE', notebook_number = 7, sim_number = 5)
WE = WE / WE[0]
plt.semilogy(t, WE, 'Db', linewidth = 2, label = 'DECSKS-2.2: s7-05')


# Linear Landau damping prediction:

gamma = -0.153359 # from linear theory for this distribution function
WE_linear = np.exp(2*gamma*t) # Linear Landau theory prediction for normalized W_E / W_E0
plt.semilogy(t, WE_linear, '--k', linewidth = 2, label = 'linear theory')


# DECSKS-1.2 data
# s07-03
WE = file_to_array(quantity = 'WE', notebook_number = 7, sim_number = 3, version = 1.2)

WE = WE / WE[0]
plt.semilogy(t, WE, '*m', markersize = 10, linewidth = 2, label = 'DECSKS-1.2: s7-03')

# s07-04
WE = file_to_array(quantity = 'WE', notebook_number = 7, sim_number = 4, version = 1.2)
WE = WE / WE[0]
plt.semilogy(t_04, WE, '-g', linewidth = 2, label = 'DECSKS-1.2: s7-04')

# s07-05
#WE = file_to_array(quantity = 'WE', notebook_number = 7, sim_number = 5, version = 1.2)
#WE = WE / WE[0]
#plt.semilogy(t, WE, 'Hy', linewidth = 2, label = 'DECSKS-1.2: s7-05')


plt.xlabel('time', fontsize = 14)
plt.ylabel('Normalized electrostatic energy, $W_E/W_{E0}$', fontsize = 12)
plt.title('DECSKS-2.2: Normalized electrostatic energy, $W_E / W_{E0}$')

plt.legend(loc = 'lower left')
plt.grid()
