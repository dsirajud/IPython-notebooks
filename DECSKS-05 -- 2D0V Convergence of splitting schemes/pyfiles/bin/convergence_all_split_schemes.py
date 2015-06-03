import matplotlib.pyplot as plt
import numpy as np


LF2_L2 = [0.0904846431239,
          0.0567891299228,
          0.0164363443106,
          0.0072981720409,
          0.00181532179268,
          0.000652533032571,
          0.000163022550367,
          4.07485580026e-05]

Y4_L2 = [0.090307529582,
      0.0567891299228,
      0.00259476847737,
      0.000508087566529,
      3.15827066908e-05,
      4.08834154401e-06
      ]

Y4_L1 = [6.93889390391e-18,
         6.41847686111e-17,
         1.19695919842e-16,
         2.15105711021e-16,
         3.81639164715e-16,
         4.7011006199e-16
         ]

O6_4_L2 = [0.00020686775511,
        9.99512075318e-06,
        5.8197029775e-07,
        1.13419682722e-07,
        7.03821650757e-09]

O6_4_L1 = [2.94902990916e-17,
         1.71737624122e-16,
         2.9316826744e-16,
         3.22658566532e-16,
         3.57353036051e-16,
         1.3357370765e-16
         ]

O11_6_L2 = [1.2851594292e-06,
            1.97273069998e-08,
            3.11795117233e-10,
            8.09948306913e-11,
            1.45263248677e-10,
            2.13965048952e-10]

O11_6_L1 = [1.38777878078e-16,
            2.53269627493e-16,
            3.22658566532e-16,
            3.29597460436e-16,
            2.20309881449e-16,
            3.1918911958e-16
            ]


O14_6_L2 = [0.000847297047215, # Nt2
            7.30905831321e-06, # Nt3
            1.25119513264e-06, # Nt4
            3.22327851392e-07, # Nt5
            4.21754379716e-08, # Nt7
            4.92326502571e-09, # Nt10
            1.64626135507e-09, # Nt12
            4.3309629071e-10, #Nt15
            9.84353656488e-11, #Nt20
            9.45205643465e-11, # Nt30
            1.73608980936e-10, #Nt60
            2.74796916074e-10] # Nt100

O14_6_L1 = [6.24500451352e-17,
            7.2858385991e-17,
            1.52655665886e-16,
            1.31838984174e-16,
            2.68882138776e-16,
            2.671474153e-16,
            3.08780778724e-16,
            3.81639164715e-16,
            4.00721122951e-16,
            4.12864187282e-16,
            4.12864187282e-16,
            2.08166817117e-16
            ]


LF2_t = [5, 10, 20, 30, 60, 100, 200, 400]
Y4_t = [5, 10, 20, 30, 60,100]
O6_4_t = [5, 10, 20, 30, 60]
O11_6_t = [5, 10, 20, 30, 60, 100]
O14_6_t = [2, 3, 4, 5, 7, 10, 12, 15, 20, 30, 60, 100]

LF2_t = np.array(LF2_t)
Y4_t = np.array(Y4_t)
O6_4_t = np.array(O6_4_t)
O11_6_t = np.array(O11_6_t)
O14_6_t = np.array(O14_6_t)

S_LF2 = 4.
S_Y4 = 8.
S_O6_4 = 13.
S_O11_6 = 23.
S_O14_6 = 29.

SNdt_LF2 = S_LF2*LF2_t
SNdt_Y4 = S_Y4*Y4_t
SNdt_O6_4 = S_O6_4*O6_4_t
SNdt_O11_6 = S_O11_6*O11_6_t
SNdt_O14_6 = S_O14_6*O14_6_t

LF2_L2 = np.array(LF2_L2)
Y4_L2 = np.array(Y4_L2)
O6_4_L2 = np.array(O6_4_L2)
O11_6_L2 = np.array(O11_6_L2)
O14_6_L2 = np.array(O14_6_L2)

plt.loglog(SNdt_LF2, LF2_L2, '-b', marker = 'o',linewidth = 3,label='LF2')
plt.loglog(SNdt_Y4, Y4_L2, '--g', marker = 'D', linewidth = 3, label = 'Y4')
plt.loglog(SNdt_O6_4, O6_4_L2, '-.r', marker = '^', linewidth = 3, label = 'O6-4')
plt.loglog(SNdt_O11_6, O11_6_L2, ':m', marker = 'v', linewidth = 3, label = 'O11-6')
plt.loglog(SNdt_O14_6, O14_6_L2, '--c', marker = '*', linewidth = 3, label = 'O14-6')

#plt.loglog(LF2_t, LF2_L2, '-b', marker = 'o',linewidth = 3,label='LF2')
#plt.loglog(Y4_t, Y4_L2, '--g', marker = 'D', linewidth = 3, label = 'Y4')
#plt.loglog(O6_4_t, O6_4_L2, '-.r', marker = '^', linewidth = 3, label = 'O6-4')
#plt.loglog(O11_6_t, O11_6_L2, ':m', marker = 'v', linewidth = 3, label = 'O11-6')
#plt.loglog(O14_6_t, O14_6_L2, '-c', marker = '*', linewidth = 3, label = 'O14-6')
plt.title('$L^2$ error, $T = 1.0, N_x = N_y = 256 = 256$', fontsize = 18)
plt.xlabel('Total number of time substeps per rotation, $SN_{\Delta t}$', fontsize = 14)
plt.legend(loc='upper right')
plt.ylabel('$L^2$ error', fontsize = 18)
plt.axis([10, 1e4, 1e-14, 2])
plt.grid()
plt.savefig('L2_error_F10Nx256Ny256T1_various_split_schemes.png')
plt.clf()

log_LF2_L2 = np.log10(LF2_L2)
log_LF2_t = np.log10(SNdt_LF2)

log_Y4_L2 = np.log10(Y4_L2)
log_Y4_t = np.log10(SNdt_Y4)

log_O6_4_L2 = np.log10(O6_4_L2)
log_O6_4_t = np.log10(SNdt_O6_4)

log_O11_6_L2 = np.log10(O11_6_L2)
log_O11_6_t = np.log10(SNdt_O11_6)

log_O14_6_L2 = np.log10(O14_6_L2)
log_O14_6_t = np.log10(SNdt_O14_6)

slope_LF2 = np.zeros(len(log_LF2_t))
slope_Y4 = np.zeros(len(log_LF2_t))
slope_O6_4 = np.zeros(len(log_O6_4_t))
slope_O11_6 = np.zeros(len(log_O11_6_t))
slope_O14_6 = np.zeros(len(log_O14_6_t))

for i in range(1,len(log_LF2_t)):
    slope_LF2[i] = ( log_LF2_L2[i] - log_LF2_L2[i-1]) / (log_LF2_t[i] - log_LF2_t[i-1])

for i in range(1,len(log_Y4_t)):
    slope_Y4[i] = ( log_Y4_L2[i] - log_Y4_L2[i-1]) / (log_Y4_t[i] - log_Y4_t[i-1])

for i in range(1,len(log_O6_4_t)):
    slope_O6_4[i] = ( log_O6_4_L2[i] - log_O6_4_L2[i-1]) / (log_O6_4_t[i] - log_O6_4_t[i-1])

for i in range(1,len(log_O11_6_t)):
    slope_O11_6[i] = ( log_O11_6_L2[i] - log_O11_6_L2[i-1]) / (log_O11_6_t[i] - log_O11_6_t[i-1])

for i in range(1,len(log_O14_6_t)):
    slope_O14_6[i] = ( log_O14_6_L2[i] - log_O14_6_L2[i-1]) / (log_O14_6_t[i] - log_O14_6_t[i-1])

#pretty print table

print "SN_dt        LF2 L2 error (max)        LF2 order (based on L2 error)"
print "---------------------------------------------------------------------"
for i in range(len(LF2_t)):
    if i == 0:
        print "%d           %2.5e                       --" % (SNdt_LF2[i], LF2_L2[i])
    else:
        print "%d         %2.5e                       %2.5f" % (SNdt_LF2[i], LF2_L2[i], slope_LF2[i])

print "\n"
print "SN_dt        Y4 L2 error (max)        Y4 order (based on L2 error)        Y4 L2 error (absmax)"
print "-----------------------------------------------------------------------------------------------"
for i in range(len(Y4_t)):
    if i == 0:
        print "%d            %2.5e                    --   " % (SNdt_Y4[i], Y4_L2[i])
    else:
        print "%d              %2.5e                 %2.5f                          %2.5e" % (SNdt_Y4[i], Y4_L2[i], slope_Y4[i], Y4_L1[i])


print "\n"
print "SN_dt        O6_4 L2 error (max)        O6_4 order (based on L2 error)        O6_4 L2 error (absmax)"
print "-----------------------------------------------------------------------------------------------"
for i in range(len(O6_4_t)):
    if i == 0:
        print "%d            %2.5e                    --   " % (SNdt_O6_4[i], O6_4_L2[i])
    else:
        print "%d              %2.5e                 %2.5f                          %2.5e" % (SNdt_O6_4[i], O6_4_L2[i], slope_O6_4[i], O6_4_L1[i])

print "\n"
print "SN_dt        O11_6 L2 error (max)        O11_6 order (based on L2 error)        O11_6 L2 error (absmax)"
print "-----------------------------------------------------------------------------------------------"
for i in range(len(O11_6_t)):
    if i == 0:
        print "%d            %2.5e                    --   " % (SNdt_O11_6[i], O11_6_L2[i])
    else:
        print "%d              %2.5e                 %2.5f                          %2.5e" % (SNdt_O11_6[i], O11_6_L2[i], slope_O11_6[i], O11_6_L1[i])

print "\n"
print "SN_dt        O14_6 L2 error (max)        O14_6 order (based on L2 error)        O14_6 L2 error (absmax)"
print "-----------------------------------------------------------------------------------------------"
for i in range(len(O14_6_t)):
    if i == 0:
        print "%d            %2.5e                    --   " % (SNdt_O14_6[i], O14_6_L2[i])
    else:
        print "%d              %2.5e                 %2.5f                          %2.5e" % (SNdt_O14_6[i], O14_6_L2[i], slope_O14_6[i], O14_6_L1[i])


methods = ['LF2', 'Y4', 'O6-4', 'O11-6', 'O14-6']
simtimes = [90.65,
            175.88,
            321.06,
            487.84,
            580.11]

simtimes_new = [43.91,
                90.21,
                153.7,
                258.6,
                329.2]

numstages = [4,
             8,
             13,
             23,
             29]


print '\n'
print 'method         proc. time per time step (matrix calc.) [sec]        Number of stages per time step'
print '------------------------------------------------------------------------------------------'

for i in range(len(methods)):
    print "%s                      %g                                           %d" % (methods[i], simtimes_new[i], numstages[i])

