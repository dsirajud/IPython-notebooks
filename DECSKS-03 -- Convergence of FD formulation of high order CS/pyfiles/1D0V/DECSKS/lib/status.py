import os
import time

other_dir = ['/lib', '/bin', '/plots', '/etc', '/etc/finite_difference_schemes', '/bin/finite_difference_schemes']
parent_dir = './..'
up_one_level = '/..'
current_dir = '.'

# Called in /DECSKS/main.py
def check_and_clean(t, n, tic, rm_plots, cleanup_report):
    """Checks the status of the simulation and broadcasts
    progress reports as percent complete

    inputs:
    t -- (instance) time
    n -- (int) current time step
    tic -- (float) epoch
    png_cleanup -- (int) user specified switch (1 = yes, 0 = no) to indicate
                    if user wants plots to be removed

    outputs:
    None
    """
    if n == t.N/10:

        print "10%% done : \
        n = %d of %d" % (n, t.N)

    if n == t.N/4:

        print "25%% done : \
        n = %d of %d" % (n, t.N)

    if n == t.N/2:

        print "50%% done : \
        n = %d of %d" % (n, t.N)

    if n == 3*t.N/4:

        print "75%% done : \
        n = %d of %d" % (n, t.N)

    if n == 9*t.N/10:

        print "90%% done : \
        n = %d of %d" % (n, t.N)

    if n == .95*t.N:

        print "95%% done : \
        n = %d of %d" % (n, t.N)

    if n == t.N:

        print "100% done : \
        end of simulation"

        toc = time.clock()     # mark ending time
        sim_time = toc - tic   # simulation processing time
        print "total time for simulation = %g seconds \n\n" % sim_time
        cleanup(rm_plots, cleanup_report)

    return None

def cleanup(rm_plots, cleanup_report):
    """Removes files (.pyc, .*~, and .png if png_cleanup == 1)

    rm_plots -- (int) user specified switch (1 = yes, 0 = no) to indicate
                    if user wants plots to be removed
    """
    if cleanup_report == 1:
        print "Cleanup started...\n\n"

    # ./DECSKS/
    pycfilelist = [pycfile for pycfile in os.listdir('.') if pycfile.endswith('.pyc')]
    cwd = os.getcwd()
    if pycfilelist != []:
        for pycfile in pycfilelist:
            os.remove(pycfile)

    if cleanup_report == 1:
        print "removing .pyc files from (%s):" % cwd
        print "-------------------------------------------------------\n"
        for pycfile in pycfilelist:
            print pycfile
        print '\n'

    os.chdir(current_dir + other_dir[0])        # ./DECSKS/lib
    cwd = os.getcwd()
    pycfilelist = [pycfile for pycfile in os.listdir('.') if pycfile.endswith('.pyc')]
    if pycfilelist != []:
        for pycfile in pycfilelist:
            os.remove(pycfile)

    if cleanup_report == 1:
        print "removing .pyc files from (%s):" % cwd
        print "-------------------------------------------------------\n"
        for pycfile in pycfilelist:
            print pycfile
        print '\n'
    os.chdir(parent_dir + other_dir[1])    # ./DECSKS/bin

    cwd = os.getcwd()
    pycfilelist = [pycfile for pycfile in os.listdir('.') if pycfile.endswith('.pyc')]
    if pycfilelist != []:
        for pycfile in pycfilelist:
            os.remove(pycfile)

    if cleanup_report == 1:
        print "removing .pyc files from (%s):" % cwd
        print "-------------------------------------------------------\n"
        for pycfile in pycfilelist:
            print pycfile

    os.chdir(parent_dir)    # ./DECSKS/
    cwd = os.getcwd()
    tmpfilelist = [tmpfile for tmpfile in os.listdir('.') if tmpfile.endswith('~')]
    if tmpfilelist != []:
        for tmpfile in tmpfilelist:
            os.remove(tmpfile)

    if cleanup_report == 1:
        print "\n"
        print "removing temporary files from (%s):" % cwd
        print "-------------------------------------------------------\n"
        for tmpfile in tmpfilelist:
            print tmpfile
        print '\n'

    os.chdir(current_dir + other_dir[0])    # ./DECSKS/lib
    cwd = os.getcwd()
    tmpfilelist = [tmpfile for tmpfile in os.listdir('.') if tmpfile.endswith('~')]
    if tmpfilelist != []:
        for tmpfile in tmpfilelist:
            os.remove(tmpfile)

    if cleanup_report == 1:
        print "\n"
        print "removing temporary files from (%s):" % cwd
        print "-------------------------------------------------------\n"
        for tmpfile in tmpfilelist:
            print tmpfile
        print '\n'

    os.chdir(parent_dir + other_dir[1])    # ./DECSKS/bin
    cwd = os.getcwd()
    tmpfilelist = [tmpfile for tmpfile in os.listdir('.') if tmpfile.endswith('~')]
    if tmpfilelist != []:
        for tmpfile in tmpfilelist:
            os.remove(tmpfile)

    if cleanup_report == 1:
        print "\n"
        print "removing temporary files from (%s):" % cwd
        print "-------------------------------------------------------\n"
        for tmpfile in tmpfilelist:
            print tmpfile
        print '\n'

    os.chdir(parent_dir + other_dir[3])    # ./DECSKS/etc
    cwd = os.getcwd()
    tmpfilelist = [tmpfile for tmpfile in os.listdir('.') if tmpfile.endswith('~')]
    if tmpfilelist != []:
        for tmpfile in tmpfilelist:
            os.remove(tmpfile)

    if cleanup_report == 1:
        print "\n"
        print "removing temporary files from (%s):" % cwd
        print "-------------------------------------------------------\n"
        for tmpfile in tmpfilelist:
            print tmpfile
        print '\n'

    os.chdir(parent_dir + other_dir[4])    # ./DECSKS/etc
    cwd = os.getcwd()
    tmpfilelist = [tmpfile for tmpfile in os.listdir('.') if tmpfile.endswith('~')]
    if tmpfilelist != []:
        for tmpfile in tmpfilelist:
            os.remove(tmpfile)

    if cleanup_report == 1:
        print "\n"
        print "removing temporary files from (%s):" % cwd
        print "-------------------------------------------------------\n"
        for tmpfile in tmpfilelist:
            print tmpfile
        print '\n'

    os.chdir(parent_dir + up_one_level + other_dir[5])    # ./DECSKS/etc
    cwd = os.getcwd()
    tmpfilelist = [tmpfile for tmpfile in os.listdir('.') if tmpfile.endswith('~')]
    if tmpfilelist != []:
        for tmpfile in tmpfilelist:
            os.remove(tmpfile)

    if cleanup_report == 1:
        print "\n"
        print "removing temporary files from (%s):" % cwd
        print "-------------------------------------------------------\n"
        for tmpfile in tmpfilelist:
            print tmpfile
        print '\n'

    if rm_plots == 1: # plot files are always generated, no check empty file list needed
        os.chdir(parent_dir + up_one_level + other_dir[2])    # ./DECSKS/plots
        cwd = os.getcwd()
        pngfilelist = [pngfile for pngfile in os.listdir('.') if pngfile.endswith('.png')]
        for pngfile in pngfilelist:
            os.remove(pngfile)

        if cleanup_report == 1:
            print "\n"
            print "removing plot files from (%s):" % cwd
            print "-------------------------------------------------------\n"
            for pngfile in pngfilelist:
                print pngfile
                os.remove(pngfile)

            print "======================================================="

            print "\ncleanup complete\n"

    return None
