One possible NON edge case convergence prescription:

(1) velocityfields.py -- set v = 0.5 in const velocity profile

(1.5) lib/diagonostics -- in L2 calc, replace x.cells with 
                          (x.cells - 0.5) for the position of 
                          the final density according to the 
                          initial density you specify

                          current: 'triple gaussian bell asymmetric'
                          (can change in params_Nx...dat files)

(2) settle on which order you wish to correct for, create
    an indicative folder, and inside ./params/ change the
    global error to the desired error (= LTE - 1)

(3) DECSKS/main_sh.py -- change the input params directory path according
                  to (2) above, near top of script

(4) DECSKS/bin/finite_differences../ generate the required FD tables needed
    for the chosen method

    python generate_tables_of_...._required_for_given_GE...py GE

(5) run DECSKS/convergence_tests.sh in shell, will output L2 error to terminal
