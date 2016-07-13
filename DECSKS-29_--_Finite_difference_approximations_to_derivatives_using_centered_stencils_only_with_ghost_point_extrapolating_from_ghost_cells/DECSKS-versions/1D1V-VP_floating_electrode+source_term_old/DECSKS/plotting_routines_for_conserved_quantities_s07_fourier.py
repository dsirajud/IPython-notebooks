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
Nt = 60 # number of time steps, does not include time zero
t = np.linspace(0,T, Nt+1)

Nt_fine = 600
t_fine = np.linspace(0,T, Nt_fine+1)

##### I1

# v2.3 data
# s07-01
I1 = file_to_array(quantity = 'I1', notebook_number = 7, sim_number = 1, version = 2.3)

I1 = np.abs((I1 - I1[0]) / I1[0])
plt.semilogy(t, I1, 'py', markersize = 12, linewidth = 2, label = 'DECSKS-2.3: s7-01')

# s07-02
I1 = file_to_array(quantity = 'I1', notebook_number = 7, sim_number = 2, version = 2.3)
I1 = np.abs((I1 - I1[0]) / I1[0])
plt.semilogy(t_fine, I1, '-b', linewidth = 2, label = 'DECSKS-2.3: s7-02')


# v2.2 data
# s07-01
I1 = file_to_array(quantity = 'I1', notebook_number = 7, sim_number = 1, version = 2.2)

I1 = np.abs((I1 - I1[0]) / I1[0])
plt.semilogy(t, I1, 'oc', markersize = 12, linewidth = 2, label = 'DECSKS-2.2: s7-01')

# s07-02
I1 = file_to_array(quantity = 'I1', notebook_number = 7, sim_number = 2, version = 2.2)
I1 = np.abs((I1 - I1[0]) / I1[0])
plt.semilogy(t_fine, I1, '-r', linewidth = 2, label = 'DECSKS-2.2: s7-02')


# DECSKS-1.2 data
# s07-01
I1 = file_to_array(quantity = 'I1', notebook_number = 7, sim_number = 1, version = 1.2)

I1 = np.abs((I1 - I1[0]) / I1[0])
plt.semilogy(t, I1, '*m', markersize = 10, linewidth = 2, label = 'DECSKS-1.2: s7-01')

# s07-02
I1 = file_to_array(quantity = 'I1', notebook_number = 7, sim_number = 2, version = 1.2)
I1 = np.abs((I1 - I1[0]) / I1[0])
plt.semilogy(t_fine, I1, '-g', linewidth = 2, label = 'DECSKS-1.2: s7-02')


plt.xlabel('time', fontsize = 14)
plt.ylabel('Abs. value of fractional error in the $L^1$ norm', fontsize = 12)
plt.title('Fractional error in the $L^1$ norm invariant (mass conservation)', fontsize = 14)

plt.axis([0,60.,1e-16, 1e-6])
plt.legend(loc = 'best')
plt.grid()
plt.figure()

##### I2

# v2.3 data
# s07-01
I2 = file_to_array(quantity = 'I2', notebook_number = 7, sim_number = 1, version = 2.3)

I2 = np.abs((I2 - I2[0]) / I2[0])
plt.semilogy(t, I2, 'py', markersize = 12, linewidth = 2, label = 'DECSKS-2.3: s7-01')

# s07-02
I2 = file_to_array(quantity = 'I2', notebook_number = 7, sim_number = 2, version = 2.3)
I2 = np.abs((I2 - I2[0]) / I2[0])
plt.semilogy(t_fine, I2, '-b', linewidth = 2, label = 'DECSKS-2.3: s7-02')



# v2.2 data
# s07-01
I2 = file_to_array(quantity = 'I2', notebook_number = 7, sim_number = 1, version = 2.2)
I2 = np.abs((I2 - I2[0]) / I2[0])
plt.semilogy(t, I2, 'oc', markersize = 12, linewidth = 2, label = 'DECSKS-2.2: s7-01')

# s07-02
I2 = file_to_array(quantity = 'I2', notebook_number = 7, sim_number = 2, version = 2.2)
I2 = np.abs((I2 - I2[0]) / I2[0])
plt.semilogy(t_fine, I2, '-r', linewidth = 2, label = 'DECSKS-2.2: s7-02')

# DECSKS-1.2 data
# s07-01
I2 = file_to_array(quantity = 'I2', notebook_number = 7, sim_number = 1, version = 1.2)

I2 = np.abs((I2 - I2[0]) / I2[0])
plt.semilogy(t, I2, '*m', markersize = 10, linewidth = 2, label = 'DECSKS-1.2: s7-01')

# s07-02
I2 = file_to_array(quantity = 'I2', notebook_number = 7, sim_number = 2, version = 1.2)
I2 = np.abs((I2 - I2[0]) / I2[0])
plt.semilogy(t_fine, I2, '-g', linewidth = 2, label = 'DECSKS-1.2: s7-02')


plt.xlabel('time', fontsize = 14)
plt.ylabel('Abs. value of fractional error in the $L^2$ norm', fontsize = 12)
plt.title('Fractional error in the $L^2$ norm invariant', fontsize = 14)
plt.axis([0,60.,1e-16, 1e-6])
plt.legend()
plt.grid()
plt.figure()



##### IW

# v2.3 data
# s07-01
IW = file_to_array(quantity = 'IW', notebook_number = 7, sim_number = 1, version = 2.3)

IW = np.abs((IW - IW[0]) / IW[0])
plt.semilogy(t, IW, 'py', markersize = 12, linewidth = 2, label = 'DECSKS-2.3: s7-01')

# s07-02
IW = file_to_array(quantity = 'IW', notebook_number = 7, sim_number = 2, version = 2.3)
IW = np.abs((IW - IW[0]) / IW[0])
plt.semilogy(t_fine, IW, '-b', linewidth = 2, label = 'DECSKS-2.3: s7-02')


# v2.2 data
# s07-01
IW = file_to_array(quantity = 'IW', notebook_number = 7, sim_number = 1, version = 2.2)
IW = np.abs((IW - IW[0]) / IW[0])
plt.semilogy(t, IW, 'oc', markersize = 12, linewidth = 2, label = 'DECSKS-2.2: s7-01')

# s07-02
IW = file_to_array(quantity = 'IW', notebook_number = 7, sim_number = 2, version = 2.2)
IW = np.abs((IW - IW[0]) / IW[0])
plt.semilogy(t_fine, IW, '-r', linewidth = 2, label = 'DECSKS-2.2: s7-02')

# DECSKS-1.2 data
# s07-01
IW = file_to_array(quantity = 'IW', notebook_number = 7, sim_number = 1, version = 1.2)

IW = np.abs((IW - IW[0]) / IW[0])
plt.semilogy(t, IW, '*m', markersize = 10, linewidth = 2, label = 'DECSKS-1.2: s7-01')

# s07-02
IW = file_to_array(quantity = 'IW', notebook_number = 7, sim_number = 2, version = 1.2)
IW = np.abs((IW - IW[0]) / IW[0])
plt.semilogy(t_fine, IW, '-g', linewidth = 2, label = 'DECSKS-1.2: s7-02')

plt.axis([0,60.,1e-16, 1e-6])
plt.xlabel('time', fontsize = 14)
plt.ylabel('Abs. value of fractional error in total energy, $W$', fontsize = 12)
plt.title('Fractional error in the total energy invariant $W$', fontsize = 14)
plt.legend(loc = 'best')
plt.grid()
plt.figure()

##### S

# v2.3 data
# s07-01
S = file_to_array(quantity = 'S', notebook_number = 7, sim_number = 1, version = 2.3)

S = np.abs((S - S[0]) / S[0])
plt.semilogy(t, S, 'py', markersize = 12, linewidth = 2, label = 'DECSKS-2.3: s7-01')

# s07-02
S = file_to_array(quantity = 'S', notebook_number = 7, sim_number = 2, version = 2.3)
S = np.abs((S - S[0]) / S[0])
plt.semilogy(t_fine, S, '-b', linewidth = 2, label = 'DECSKS-2.3: s7-02')


# v2.2 data
# s07-01
S = file_to_array(quantity = 'S', notebook_number = 7, sim_number = 1, version = 2.2)
S = np.abs((S - S[0]) / S[0])
plt.semilogy(t, S, 'oc', markersize = 12, linewidth = 2, label = 'DECSKS-2.2: s7-01')

# s07-02
S = file_to_array(quantity = 'S', notebook_number = 7, sim_number = 2, version = 2.2)
S = np.abs((S - S[0]) / S[0])
plt.semilogy(t_fine, S, '-r', linewidth = 2, label = 'DECSKS2.2: s7-02')

# DECSKS-1.2 data
# s07-01
S = file_to_array(quantity = 'S', notebook_number = 7, sim_number = 1, version = 1.2)

S = np.abs((S - S[0]) / S[0])
plt.semilogy(t, S, '*m', markersize = 10, linewidth = 2, label = 'DECSKS-1.2: s7-01')

# s07-02
S = file_to_array(quantity = 'S', notebook_number = 7, sim_number = 2, version = 1.2)
S = np.abs((S - S[0]) / S[0])
plt.semilogy(t_fine, S, '-g', linewidth = 2, label = 'DECSKS-1.2: s7-02')

plt.xlabel('time', fontsize = 14)
plt.ylabel('Abs. value of fractional error in the entropy $S$ invariant', fontsize = 12)
plt.title('Fractional error in the entropy invariant $S$', fontsize = 14)

plt.axis([0,60.,1e-16, 1e-6])
plt.legend(loc = 'best')
plt.grid()
plt.figure()

##### WE


# v2.3 data, s07-01 is plotted last so that the smaller sized markers can be overlay on top, not hidden underneath
# s07-02
WE = file_to_array(quantity = 'WE', notebook_number = 7, sim_number = 2, version = 2.3)
WE = WE / WE[0]
plt.semilogy(t_fine, WE, '-.y', linewidth = 2, label = 'DECSKS-2.3: s7-02')

# v2.2 data
# s07-01
WE = file_to_array(quantity = 'WE', notebook_number = 7, sim_number = 1, version = 2.2)
WE = WE / WE[0]
plt.semilogy(t, WE, 'oc', markersize = 12, linewidth = 2, label = 'DECSKS-2.2: s7-01')


# s07-02
WE = file_to_array(quantity = 'WE', notebook_number = 7, sim_number = 2, version = 2.2)
WE = WE / WE[0]
plt.semilogy(t_fine, WE, '--r', linewidth = 2, label = 'DECSKS-2.2: s7-02')


# Linear Landau damping prediction:

gamma = -0.153359 # from linear theory for this distribution function
WE_linear = np.exp(2*gamma*t) # Linear Landau theory prediction for normalized W_E / W_E0
plt.semilogy(t, WE_linear, '--k', linewidth = 2, label = 'linear theory')


# DECSKS-1.2 data
# s07-01
WE = file_to_array(quantity = 'WE', notebook_number = 7, sim_number = 1, version = 1.2)

WE = WE / WE[0]
plt.semilogy(t, WE, '*m', markersize = 11, linewidth = 2, label = 'DECSKS-1.2: s7-01')

# s07-02
WE = file_to_array(quantity = 'WE', notebook_number = 7, sim_number = 2, version = 1.2)
WE = WE / WE[0]
plt.semilogy(t_fine, WE, ':g', linewidth = 2, label = 'DECSKS-1.2: s7-02')

# v2.3 data
# s07-01
WE = file_to_array(quantity = 'WE', notebook_number = 7, sim_number = 1, version = 2.3)

WE = WE / WE[0]
plt.semilogy(t, WE, 'py', markersize = 5, linewidth = 2, label = 'DECSKS-2.3: s7-01')


plt.xlabel('time', fontsize = 14)
plt.ylabel('Normalized electrostatic energy, $W_E/W_{E0}$', fontsize = 12)
plt.title('DECSKS-2.2: Normalized electrostatic energy, $W_E / W_{E0}$')

plt.legend()
plt.grid()
