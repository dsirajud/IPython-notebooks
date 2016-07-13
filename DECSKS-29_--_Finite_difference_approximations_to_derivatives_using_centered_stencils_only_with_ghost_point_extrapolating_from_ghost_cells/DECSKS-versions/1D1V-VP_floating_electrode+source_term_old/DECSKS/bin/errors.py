from main import * # numpy imported as np, matplotlib.pyplot as plt
import matplotlib.pyplot as plt
import math as math

 

# REMEMBER TO CHANGE THIS SO L2 IS BIG ENOUGH
NO = 7 # number of runs, O = order

#---------------------------------------------------------------------------#
#outfile and plot strings

pound    =  '#'
eq      =  '='
dash    =  '-'
dot     =  '.'
spc     =  ' '
underscore   = '_'
vertbar = '|'
newline = '\n'
indent  = '    '    # i.e. ctrl + u + spc
hdashline  = spc  + spc + dash*75 + spc + spc
hpounddashline = pound + spc + dash * 100 + spc + pound
hdotline   = dot * 10
hdotlongline = spc * 2 + dot * 75 + spc * 2
domain_title = spc * (len(header)/2) + 'Error reports N_F8' + spc * (len(header) / 2)
header     = pound + spc + eq * 75  + spc + pound
footer = header
hardcopy_stem = 'plot_-_N_F8_' # for eps plot save at final time
eps_format = '.eps'

outfilename = 'error_report_F8_N'
out_simdata_stem = 'F8_N_simdata_'

outfile = open(outfilename,'w')
outfile.write(header + newline)
outfile.write(domain_title + newline)
outfile.write(hdashline + newline)

#===========================================================================#
# error calculations and write operations

L2 = np.zeros(NO)

#L2[0], x_exact, f_exact, x, f_CS, Nx, Nt, sim_time  = main(8,25)
#L2_str = 'L2[0] = %1.9e' % L2[0]
#time_str = 'total sim time = %1.9e seconds,  ' % sim_time
#Nx_str = 'Nx8'
#Nt_str = 'Nt25'
#NxNt_str = 'for ' + Nx_str + Nt_str

#plt.plot(x,f_CS,'ob', label = "F8")
#plt.plot(x_exact,f_exact,'g', linewidth = 2.0, label = "Exact initial density")
#plt.grid()
#plt.axis([-.5,.5,0,1.2])
#plt.xlabel('$x$', fontsize = 18)
#plt.ylabel('$f(x)$', fontsize = 18)
#plt.legend(bbox_to_anchor = (0., 1.02, 1., .102), loc = 3,borderaxespad = 0., fancybox = True, shadow = True, ncol = 5)
#hardcopy = hardcopy_stem + underscore + Nx_str + Nt_str + '_w_f0_itmax' + eps_format #at time t = T (it_max)
#plt.text(0.2,1.14, '$N_x = %d, \, N_t = %d$' % (Nx,Nt), fontsize = 16)
#plt.savefig(hardcopy)
#plt.clf()

#outfile.write(NxNt_str + spc + ':' + spc)
#outfile.write(time_str + L2_str + newline)

# write final vector data to files

#out2name = out_simdata_stem + Nx_str + Nt_str + underscore + 'f_CS'
#out3name = out_simdata_stem + Nx_str + Nt_str + underscore + 'x_CS'
#out4name = out_simdata_stem + Nx_str + Nt_str + underscore + 'f_exact'
#out5name = out_simdata_stem + Nx_str + Nt_str + underscore + 'x_exact'

#out2 = open(out2name,'w')
#out3 = open(out3name,'w')
#out4 = open(out4name,'w')
#out5 = open(out5name,'w')

#write all vectors to separate files in a single column (i.e. with \n chars)

#f_CS.tofile(out2,'\n')
#x.tofile(out3,'\n')
#f_exact.tofile(out4,'\n')
#x_exact.tofile(out5,'\n')

# see "read_files_into_an_array.py" for recipe of extracting arrays interactively for plotting, etc.

#out2.close(), out3.close(), out4.close(), out5.close()


#---------------------------------------------------------------------------#
#---------------------------------------------------------------------------#
#L2[1], x_exact, f_exact, x, f_CS, Nx, Nt, sim_time = main(16,50)
#L2_str = 'L2[1] = %1.9e' % L2[1]
#time_str = 'total sim time = %1.9e seconds,  ' % sim_time
#Nx_str = 'Nx16'
#Nt_str = 'Nt50'
#NxNt_str = 'for ' + Nx_str + Nt_str

#plt.plot(x,f_CS,'ob', label = "F8")
#plt.plot(x_exact,f_exact,'g', linewidth = 2.0, label = "Exact initial density")
#plt.grid()
#plt.axis([-.5,.5,0,1.2])
#plt.xlabel('$x$', fontsize = 18)
#plt.ylabel('$f(x)$', fontsize = 18)
#plt.legend(bbox_to_anchor = (0., 1.02, 1., .102), loc = 3,borderaxespad = 0., fancybox = True, shadow = True, ncol = 5)
#hardcopy = hardcopy_stem + underscore + Nx_str + Nt_str + '_w_f0_itmax' + eps_format #at time t = T (it_max)
#plt.text(0.2,1.14, '$N_x = %d, \, N_t = %d$' % (Nx,Nt), fontsize = 16)
#plt.savefig(hardcopy)
#plt.clf()

#outfile.write(NxNt_str + spc + ':' + spc)
#outfile.write(time_str + L2_str + newline)

# write final vector data to files

#out2name = out_simdata_stem + Nx_str + Nt_str + underscore + 'f_CS'
#out3name = out_simdata_stem + Nx_str + Nt_str + underscore + 'x_CS'
#out4name = out_simdata_stem + Nx_str + Nt_str + underscore + 'f_exact'
#out5name = out_simdata_stem + Nx_str + Nt_str + underscore + 'x_exact'

#out2 = open(out2name,'w')
#out3 = open(out3name,'w')
#out4 = open(out4name,'w')
#out5 = open(out5name,'w')

#write all vectors to separate files in a single column (i.e. with \n chars)

#f_CS.tofile(out2,'\n')
#x.tofile(out3,'\n')
#f_exact.tofile(out4,'\n')
#x_exact.tofile(out5,'\n')

# see "read_files_into_an_array.py" for recipe of extracting arrays interactively for plotting, etc.

#out2.close(), out3.close(), out4.close(), out5.close()

#---------------------------------------------------------------------------#
#---------------------------------------------------------------------------#
#L2[2], x_exact, f_exact, x, f_CS, Nx, Nt, sim_time = main(32,100)
#L2_str = 'L2[2] = %1.9e' % L2[2]
#time_str = 'total sim time = %1.9e seconds,  ' % sim_time
#Nx_str = 'Nx32'
#Nt_str = 'Nt100'
#NxNt_str = 'for ' + Nx_str + Nt_str

#plt.plot(x,f_CS,'ob', label = "F8")
#plt.plot(x_exact,f_exact,'g', linewidth = 2.0, label = "Exact initial density")
#plt.grid()
#plt.axis([-.5,.5,0,1.2])
#plt.xlabel('$x$', fontsize = 18)
#plt.ylabel('$f(x)$', fontsize = 18)
#plt.legend(bbox_to_anchor = (0., 1.02, 1., .102), loc = 3,borderaxespad = 0., fancybox = True, shadow = True, ncol = 5)
#hardcopy = hardcopy_stem + underscore + Nx_str + Nt_str + '_w_f0_itmax' + eps_format #at time t = T (it_max)
#plt.text(0.2,1.14, '$N_x = %d, \, N_t = %d$' % (Nx,Nt), fontsize = 16)
#plt.savefig(hardcopy)
#plt.clf()

#outfile.write(NxNt_str + spc + ':' + spc)
#outfile.write(time_str + L2_str + newline)

# write final vector data to files

#out2name = out_simdata_stem + Nx_str + Nt_str + underscore + 'f_CS'
#out3name = out_simdata_stem + Nx_str + Nt_str + underscore + 'x_CS'
#out4name = out_simdata_stem + Nx_str + Nt_str + underscore + 'f_exact'
#out5name = out_simdata_stem + Nx_str + Nt_str + underscore + 'x_exact'

#out2 = open(out2name,'w')
#out3 = open(out3name,'w')
#out4 = open(out4name,'w')
#out5 = open(out5name,'w')

#write all vectors to separate files in a single column (i.e. with \n chars)

#f_CS.tofile(out2,'\n')
#x.tofile(out3,'\n')
#f_exact.tofile(out4,'\n')
#x_exact.tofile(out5,'\n')

# see "read_files_into_an_array.py" for recipe of extracting arrays interactively for plotting, etc.

#out2.close(), out3.close(), out4.close(), out5.close()

#---------------------------------------------------------------------------#
#---------------------------------------------------------------------------#
L2[3], x_exact, f_exact, x, f_CS, Nx, Nt, sim_time = main(128,128,100)
L2_str = 'L2[3] = %1.9e' % L2[3]
time_str = 'total sim time = %1.9e seconds,  ' % sim_time
Nx_str = 'Nx128'
Ny_str = 'Ny128'
Nt_str = 'Nt100'
NxNt_str = 'for ' + Nx_str + Nt_str

plt.plot(x,f_CS,'ob', label = "F8")
plt.plot(x_exact,f_exact,'g', linewidth = 2.0, label = "Exact initial density")
plt.grid()
plt.axis([-.5,.5,0,1.2])
plt.xlabel('$x$', fontsize = 18)
plt.ylabel('$f(x)$', fontsize = 18)
plt.legend(bbox_to_anchor = (0., 1.02, 1., .102), loc = 3,borderaxespad = 0., fancybox = True, shadow = True, ncol = 5)
hardcopy = hardcopy_stem + underscore + Nx_str + Nt_str + '_w_f0_itmax' + eps_format #at time t = T (it_max)
plt.text(0.2,1.14, '$N_x = %d, \, N_t = %d$' % (Nx,Nt), fontsize = 16)
plt.savefig(hardcopy)
plt.clf()

outfile.write(NxNt_str + spc + ':' + spc)
outfile.write(time_str + L2_str + newline)

# write final vector data to files

out2name = out_simdata_stem + Nx_str + Nt_str + underscore + 'f_CS'
out3name = out_simdata_stem + Nx_str + Nt_str + underscore + 'x_CS'
out4name = out_simdata_stem + Nx_str + Nt_str + underscore + 'f_exact'
out5name = out_simdata_stem + Nx_str + Nt_str + underscore + 'x_exact'

out2 = open(out2name,'w')
out3 = open(out3name,'w')
out4 = open(out4name,'w')
out5 = open(out5name,'w')

#write all vectors to separate files in a single column (i.e. with \n chars)

f_CS.tofile(out2,'\n')
x.tofile(out3,'\n')
f_exact.tofile(out4,'\n')
x_exact.tofile(out5,'\n')

from scitools.easyviz import *

stem = 'plot_-_N_F8_'
suffix = '_w_f0_it0*'

filename = stem + Nx_str + Ny_str + Nt_str + suffix + eps_format

movie(filename,encoder='ppmtompeg',fps=24,output_file='tmpmovie.mpeg')

#x = linspace(-.5,.5,1024)
#vv = v_2pix(x)
#clf()
#plt.plot(x,vv,'g', linewidth = 3.0, label = "$v_{LCS}$ distribution")
#plt.grid()
#plt.axis([-.5,.5,-1.1,1.1])
#plt.xlabel('$x$')
#plt.xlabel(prop = {fontsize: '18'})
#plt.ylabel('$v_{LCS}(x)$')
#plt.ylabel(prop = {fontsize: '18'})
#plt.legend(bbox_to_anchor = (0., 1.02, 1., .102), loc = 3,borderaxespad = 0., fancybox = True, shadow = True, ncol = 5)
#plt.savefig('v_LCS.eps')

# see "read_files_into_an_array.py" for recipe of extracting arrays interactively for plotting, etc.

out2.close(), out3.close(), out4.close(), out5.close()

#---------------------------------------------------------------------------#
#---------------------------------------------------------------------------#
#L2[4], x_exact, f_exact, x, f_CS, Nx, Nt, sim_time = main(128,400)
#L2_str = 'L2[4] = %1.9e' % L2[4]
#time_str = 'total sim time = %1.9e seconds,  ' % sim_time
#Nx_str = 'Nx128'
#Nt_str = 'Nt400'
#NxNt_str = 'for ' + Nx_str + Nt_str

#plt.plot(x,f_CS,'ob', label = "F8")
#plt.plot(x_exact,f_exact,'g', linewidth = 2.0, label = "Exact initial density")
#plt.grid()
#plt.axis([-.5,.5,0,1.2])
#plt.xlabel('$x$', fontsize = 18)
#plt.ylabel('$f(x)$', fontsize = 18)
#plt.legend(bbox_to_anchor = (0., 1.02, 1., .102), loc = 3,borderaxespad = 0., fancybox = True, shadow = True, ncol = 5)
#hardcopy = hardcopy_stem + underscore + Nx_str + Nt_str + '_w_f0_itmax' + eps_format #at time t = T (it_max)
#plt.text(0.2,1.14, '$N_x = %d, \, N_t = %d$' % (Nx,Nt), fontsize = 16)
#plt.savefig(hardcopy)
#plt.clf()

#outfile.write(NxNt_str + spc + ':' + spc)
#outfile.write(time_str + L2_str + newline)

# write final vector data to files

#out2name = out_simdata_stem + Nx_str + Nt_str + underscore + 'f_CS'
#out3name = out_simdata_stem + Nx_str + Nt_str + underscore + 'x_CS'
#out4name = out_simdata_stem + Nx_str + Nt_str + underscore + 'f_exact'
#out5name = out_simdata_stem + Nx_str + Nt_str + underscore + 'x_exact'

#out2 = open(out2name,'w')
#out3 = open(out3name,'w')
#out4 = open(out4name,'w')
#out5 = open(out5name,'w')

#write all vectors to separate files in a single column (i.e. with \n chars)

#f_CS.tofile(out2,'\n')
#x.tofile(out3,'\n')
#f_exact.tofile(out4,'\n')
#x_exact.tofile(out5,'\n')

# see "read_files_into_an_array.py" for recipe of extracting arrays interactively for plotting, etc.

#out2.close(), out3.close(), out4.close(), out5.close()

#---------------------------------------------------------------------------#
#---------------------------------------------------------------------------#
#L2[5], x_exact, f_exact, x, f_CS, Nx, Nt, sim_time = main(256,800)
#L2_str = 'L2[5] = %1.9e' % L2[5]
#time_str = 'total sim time = %1.9e seconds,  ' % sim_time
#Nx_str = 'Nx256'
#Nt_str = 'Nt800'
#NxNt_str = 'for ' + Nx_str + Nt_str

#plt.plot(x,f_CS,'ob', label = "F8")
#plt.plot(x_exact,f_exact,'g', linewidth = 2.0, label = "Exact initial density")
#plt.grid()
#plt.axis([-.5,.5,0,1.2])
#plt.xlabel('$x$', fontsize = 18)
#plt.ylabel('$f(x)$', fontsize = 18)
#plt.legend(bbox_to_anchor = (0., 1.02, 1., .102), loc = 3,borderaxespad = 0., fancybox = True, shadow = True, ncol = 5)
#hardcopy = hardcopy_stem + underscore + Nx_str + Nt_str + '_w_f0_itmax' + eps_format #at time t = T (it_max)
#plt.text(0.2,1.14, '$N_x = %d, \, N_t = %d$' % (Nx,Nt), fontsize = 16)
#plt.savefig(hardcopy)
#plt.clf()

#outfile.write(NxNt_str + spc + ':' + spc)
#outfile.write(time_str + L2_str + newline)

# write final vector data to files

#out2name = out_simdata_stem + Nx_str + Nt_str + underscore + 'f_CS'
#out3name = out_simdata_stem + Nx_str + Nt_str + underscore + 'x_CS'
#out4name = out_simdata_stem + Nx_str + Nt_str + underscore + 'f_exact'
#out5name = out_simdata_stem + Nx_str + Nt_str + underscore + 'x_exact'

#out2 = open(out2name,'w')
#out3 = open(out3name,'w')
#out4 = open(out4name,'w')
#out5 = open(out5name,'w')

#write all vectors to separate files in a single column (i.e. with \n chars)

#f_CS.tofile(out2,'\n')
#x.tofile(out3,'\n')
#f_exact.tofile(out4,'\n')
#x_exact.tofile(out5,'\n')

# see "read_files_into_an_array.py" for recipe of extracting arrays interactively for plotting, etc.

#out2.close(), out3.close(), out4.close(), out5.close()

#---------------------------------------------------------------------------#
#---------------------------------------------------------------------------#
#L2[6], x_exact, f_exact, x, f_CS, Nx, Nt, sim_time = main(512,1600)
#L2_str = 'L2[6] = %1.9e' % L2[6]
#time_str = 'total sim time = %1.9e seconds,  ' % sim_time
#Nx_str = 'Nx512'
#Nt_str = 'Nt1600'
#NxNt_str = 'for ' + Nx_str + Nt_str

#plt.plot(x,f_CS,'ob', label = "F8")
#plt.plot(x_exact,f_exact,'g', linewidth = 2.0, label = "Exact initial density")
#plt.grid()
#plt.axis([-.5,.5,0,1.2])
#plt.xlabel('$x$', fontsize = 18)
#plt.ylabel('$f(x)$', fontsize = 18)
#plt.legend(bbox_to_anchor = (0., 1.02, 1., .102), loc = 3,borderaxespad = 0., fancybox = True, shadow = True, ncol = 5)
#hardcopy = hardcopy_stem + underscore + Nx_str + Nt_str + '_w_f0_itmax' + eps_format #at time t = T (it_max)
#plt.text(0.2,1.14, '$N_x = %d, \, N_t = %d$' % (Nx,Nt), fontsize = 16)
#plt.savefig(hardcopy)
#plt.clf()

#outfile.write(NxNt_str + spc + ':' + spc)
#outfile.write(time_str + L2_str + newline)

# write final vector data to files

#out2name = out_simdata_stem + Nx_str + Nt_str + underscore + 'f_CS'
#out3name = out_simdata_stem + Nx_str + Nt_str + underscore + 'x_CS'
#out4name = out_simdata_stem + Nx_str + Nt_str + underscore + 'f_exact'
#out5name = out_simdata_stem + Nx_str + Nt_str + underscore + 'x_exact'

#out2 = open(out2name,'w')
#out3 = open(out3name,'w')
#out4 = open(out4name,'w')
#out5 = open(out5name,'w')

#write all vectors to separate files in a single column (i.e. with \n chars)

#f_CS.tofile(out2,'\n')
#x.tofile(out3,'\n')
#f_exact.tofile(out4,'\n')
#x_exact.tofile(out5,'\n')

# see "read_files_into_an_array.py" for recipe of extracting arrays interactively for plotting, etc.

#out2.close(), out3.close(), out4.close(), out5.close()


#---------------------------------------------------------------------------#
#---------------------------------------------------------------------------#
#L2[7], x_exact, f_exact, x, f_CS, Nx, Nt, sim_time = main(1024,3200)
#L2_str = 'L2[7] = %1.9e' % L2[7]
#time_str = 'total sim time = %1.9e seconds,  ' % sim_time
#Nx_str = 'Nx1024'
#Nt_str = 'Nt3200'
#NxNt_str = 'for ' + Nx_str + Nt_str

#plt.plot(x,f_CS,'ob', label = "F8")
#plt.plot(x_exact,f_exact,'g', linewidth = 2.0, label = "Exact initial density")
#plt.grid()
#plt.axis([-.5,.5,0,1.2])
#plt.xlabel('$x$', fontsize = 18)
#plt.ylabel('$f(x)$', fontsize = 18)
#plt.legend(bbox_to_anchor = (0., 1.02, 1., .102), loc = 3,borderaxespad = 0., fancybox = True, shadow = True, ncol = 5)
#hardcopy = hardcopy_stem + underscore + Nx_str + Nt_str + '_w_f0_itmax' + eps_format #at time t = T (it_max)
#plt.text(0.2,1.14, '$N_x = %d, \, N_t = %d$' % (Nx,Nt), fontsize = 16)
#plt.savefig(hardcopy)
#plt.clf()

#outfile.write(NxNt_str + spc + ':' + spc)
#outfile.write(time_str + L2_str + newline)

# write final vector data to files

#out2name = out_simdata_stem + Nx_str + Nt_str + underscore + 'f_CS'
#out3name = out_simdata_stem + Nx_str + Nt_str + underscore + 'x_CS'
#out4name = out_simdata_stem + Nx_str + Nt_str + underscore + 'f_exact'
#out5name = out_simdata_stem + Nx_str + Nt_str + underscore + 'x_exact'

#out2 = open(out2name,'w')
#out3 = open(out3name,'w')
#out4 = open(out4name,'w')
#out5 = open(out5name,'w')

#write all vectors to separate files in a single column (i.e. with \n chars)

#f_CS.tofile(out2,'\n')
#x.tofile(out3,'\n')
#f_exact.tofile(out4,'\n')
#x_exact.tofile(out5,'\n')

# see "read_files_into_an_array.py" for recipe of extracting arrays interactively for plotting, etc.

#out2.close(), out3.close(), out4.close(), out5.close()

#---------------------------------------------------------------------------#
#---------------------------------------------------------------------------#
#L2[8], x_exact, f_exact, x, f_CS, Nx, Nt, sim_time = main(2048,6400) stored as L2_6
#---------------------------------------------------------------------------#
#---------------------------------------------------------------------------#
#L2[9], x_exact, f_exact, x, f_CS, Nx, Nt, sim_time = main(4096,100)

outfile.write(newline + newline)
outfile.write(hdashline + newline)
outfile.write(newline)

#for q in range(1,NO):
#    print "for q = %d" % q
#    order =  math.log(L2[q-1] / L2[q] , 2)
#    order_str = 'L2[%d]/L2[%d] order = %1.9e' % ((q-1), q, order)
#    outfile.write(order_str + newline)

outfile.close()
